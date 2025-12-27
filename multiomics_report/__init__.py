import os
import csv
import logging
from multiqc import config
from pkg_resources import get_distribution

log = logging.getLogger("multiqc")

# Save this plugin's version number (defined in setup.py) to the MultiQC config
try:
    config.multiomics_report_version = get_distribution("multiomics_report").version
except Exception:
    config.multiomics_report_version = "unknown"

# MODULE_DEFAULTS should use a primitive value (not a list) for 'use_filename_as_sample_name'
MODULE_DEFAULTS = {
    'fastp': {
        's_name_filenames': True,
    },
}


def before_config():
    """
    Hook that runs before config is finalized (before file search).
    This is where we set up sample renaming rules.
    
    IMPORTANT: This runs BEFORE file search, so patterns are set before
    sample names are cleaned from filenames. This is the correct place to
    set fn_clean_exts.
    """
    # --- 1. Clean Extensions (Crucial step!) ---
    # Strip common extensions so filenames become clean sample names
    # e.g., 'test1_1_fastp.json' -> 'test1_1' (after _fastp pattern matches)
    fn_clean_enabled = getattr(config, 'fn_clean_sample_names', True)
    log.info(f"Plugin: fn_clean_sample_names = {fn_clean_enabled}")
    if not fn_clean_enabled:
        log.warning(
            "Plugin: fn_clean_sample_names is False! "
            "Patterns won't be applied. Enable it in config or set it to True."
        )
        # Optionally enable it
        config.fn_clean_sample_names = True
        log.info("Plugin: Enabled fn_clean_sample_names")
    
    # Initialize fn_clean_exts if it doesn't exist
    if not hasattr(config, 'fn_clean_exts'):
        config.fn_clean_exts = []
    
    # Add RNA-seq specific extensions to clean
    # CRITICAL: Order matters! Patterns are checked in order, and the FIRST match wins.
    # More specific patterns MUST come before less specific ones.
    # Example: '_fastp' must come before '.json' so 'test1_fastp.json' -> 'test1' not 'test1_fastp'
    extra_exts = [
        '_Log.final.out',  # Most specific first
        '.metrics.tsv',    # More specific before '.metrics'
        '.metrics',
        '.isoforms',
        '.genes',
        '_fastp',          # Must come before .json, .html, etc. in defaults
    ]
    
    # Prepend to the beginning of the list (checked first)
    config.fn_clean_exts[0:0] = extra_exts
    
    log.debug(f"Added {len(extra_exts)} patterns to fn_clean_exts")
    
    # --- 2. Dynamic Sample Renaming ---
    # Find and parse the sample annotation CSV
    annotation_file = find_sample_sheet()
    
    if annotation_file:
        new_rules_dict = parse_sample_sheet_dict(annotation_file)
        
        # Generate a TSV file for --replace-names
        # This is the most reliable way to rename samples in MultiQC
        tsv_path = generate_rename_tsv(new_rules_dict, annotation_file)
        
        if tsv_path:
            # Load the rename file using MultiQC's built-in function
            # This is what --replace-names does internally
            try:
                config.load_replace_names(tsv_path)
                print(f"Plugin Info: Loaded sample renaming file with {len(new_rules_dict)} rules from {tsv_path}")
                # Print first few examples
                for i, (old, new) in enumerate(list(new_rules_dict.items())[:3]):
                    print(f"  {old} -> {new}")
                if len(new_rules_dict) > 3:
                    print(f"  ... and {len(new_rules_dict) - 3} more")
            except Exception as e:
                print(f"Plugin Error: Failed to load rename file: {e}")
                import traceback
                traceback.print_exc()
    else:
        print(f"Plugin Warning: Sample annotation file not found. Skipping sample renaming.")
    
    # Set default config for modules
    for module_name, module_config in MODULE_DEFAULTS.items():
        # Only set if not already configured by user
        if not hasattr(config, module_name):
            config.update({module_name: module_config})
            log.debug(f"Set default config for {module_name}: {module_config}")
        else:
            # Merge or overwrite based on type
            existing = getattr(config, module_name, None)
            if isinstance(module_config, dict) and isinstance(existing, dict):
                # Merge: user config overrides defaults
                merged = {**module_config, **existing}
            else:
                # If the user has a value (existing), keep it. If not, use the default.
                merged = existing if existing is not None else module_config
            config.update({module_name: merged})
            log.debug(f"Merged config for {module_name}: {merged}")

    final_patterns = list(config.fn_clean_exts) if hasattr(config, 'fn_clean_exts') else []

    log.info(f"Plugin: Set fn_clean_exts to {len(final_patterns)} patterns total")
    log.debug(f"Plugin: First 10 patterns: {final_patterns[:10]}")
    log.debug(f"Plugin: Your custom patterns are at positions: {[final_patterns.index(p) for p in extra_exts if p in final_patterns]}")


def config_loaded():
    """
    Hook that runs after all configs are loaded.
    Verify that fn_clean_exts was set correctly.
    """
    # Just verify - the actual setting happens in before_config
    if hasattr(config, 'fn_clean_exts'):
        final_patterns = list(config.fn_clean_exts)
        log.debug(f"Plugin: Verified fn_clean_exts has {len(final_patterns)} patterns")
        if config.verbose:
            log.debug(f"Plugin: First 15 patterns: {final_patterns[:15]}")
    else:
        log.warning("Plugin: fn_clean_exts not set! This should not happen.")


def execution_start():
    """ 
    Code to run when MultiQC starts (after config is loaded).
    Register search patterns for our custom modules.
    """
    from multiqc.utils.util_functions import update_dict
    
    # Register search patterns
    search_patterns = {
        'rnaseq/rnaseqqc': {
            'fn': '*metrics.tsv',
        },
        'rnaseq/gene_type_counts': {
            'fn': '*.json',
            'contents': 'gene_type_count'
        },
        "rsem": [
            {
                "fn": "*.json",
                "exclude": None,
                'contents': 'num_genes_detected'
            },
            {
                "fn": "*.cnt",
                "exclude": None,
                'contents': 'num_genes_detected',
            }
        ],
        'rnaseq/mad_rna' : {
            'fn': '*summary_metrics.json', 
            'contents': 'MAD of log ratios'
        }
    }
    # Merge search patterns into config.sp (adds to beginning to supersede defaults)
    config.sp = update_dict(config.sp, search_patterns, add_in_the_beginning=True)

def find_sample_sheet():
    """
    Locate the sample annotation CSV file with priority logic.
    Priority: 1) Environment variable, 2) Analysis directories, 3) Output directory, 4) Current directory
    """
    # Check environment variable first (set by Snakemake)
    if 'MULTIQC_SAMPLE_SHEET' in os.environ:
        env_path = os.environ['MULTIQC_SAMPLE_SHEET']
        if os.path.exists(env_path):
            return env_path
    
    # Check in analysis directories (where MultiQC searches for files)
    analysis_dirs = getattr(config, 'analysis_dir', None)
    if analysis_dirs:
        # Handle both list and single string
        if isinstance(analysis_dirs, str):
            analysis_dirs = [analysis_dirs]
        elif not isinstance(analysis_dirs, list):
            analysis_dirs = []
        
        # Check each analysis directory for sample_annotation.csv
        for analysis_dir in analysis_dirs:
            if isinstance(analysis_dir, str) and os.path.exists(analysis_dir):
                potential_path = os.path.join(analysis_dir, 'sample_annotation.csv')
                if os.path.exists(potential_path):
                    return potential_path
    
    # Check in the output directory (where we copy the file)
    output_dir = getattr(config, 'output_dir', None)
    if output_dir and isinstance(output_dir, str) and os.path.exists(output_dir):
        potential_path = os.path.join(output_dir, 'sample_annotation.csv')
        if os.path.exists(potential_path):
            return potential_path
    
    # Check in current working directory as last resort
    potential_path = os.path.abspath('sample_annotation.csv')
    if os.path.exists(potential_path):
        return potential_path
    
    return None

def parse_sample_sheet_dict(filepath):
    """
    Reads the sample annotation CSV and returns a dictionary of renaming rules.
    
    Returns: {'test1_1': 'test1', 'test1_2': 'test1', ...}
    
    This format is what MultiQC's config.sample_names_replace expects.
    """
    rules = {}
    try:
        with open(filepath, 'r') as f:
            reader = csv.DictReader(f)
            for row in reader:
                # Get the target name (the final merged name)
                target_name = row.get('sample_name', '').strip()
                
                # Get the run number
                run_id = row.get('run', '').strip()
                
                # Skip if either is missing
                if not target_name or not run_id:
                    continue
                
                # Construct the source name (what MultiQC will see after cleaning)
                # e.g., 'test1_1_fastp.json' becomes 'test1_1' after cleaning, then we rename to 'test1'
                source_name = f"{target_name}_{run_id}"
                
                # Add the renaming rule: source -> target
                rules[source_name] = target_name
                
    except Exception as e:
        print(f"Plugin Error: Failed to parse sample sheet {filepath}: {e}")
        import traceback
        traceback.print_exc()
    
    return rules

def generate_rename_tsv(rules_dict, annotation_file):
    """
    Generates a TSV file for MultiQC's --replace-names option.
    Returns the path to the generated TSV file.
    """
    try:
        # Create TSV file in the same directory as the annotation file
        tsv_path = annotation_file.replace('.csv', '_multiqc_rename.tsv')
        
        with open(tsv_path, 'w') as f:
            for source, target in rules_dict.items():
                f.write(f"{source}\t{target}\n")
        
        return tsv_path
    except Exception as e:
        print(f"Plugin Error: Failed to generate rename TSV: {e}")
        import traceback
        traceback.print_exc()
        return None

def multiomics_report_after_modules():
    """
    Hook that runs after all modules are initialized.
    This allows us to access data from other modules like multiomics.
    """
    # This hook runs after modules, but we'll handle data access in the module itself
    # using a different approach
    pass