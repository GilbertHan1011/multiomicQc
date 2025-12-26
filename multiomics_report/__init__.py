import os
import csv
from multiqc.utils import config

def before_config():
    """
    Hook that runs before config is finalized (before file search).
    This is where we set up sample renaming rules.
    """
    # --- 1. Clean Extensions (Crucial step!) ---
    # Strip common extensions so filenames become clean sample names
    # e.g., 'test1_1_fastp.json' -> 'test1_1'
    if not hasattr(config, 'fn_clean_exts'):
        config.fn_clean_exts = []
    
    # Add RNA-seq specific extensions to clean
    extra_exts = [
        '_fastp',
        '.metrics',
        '.metrics.tsv',
        '_Log.final.out',
        '.genes',
        '.isoforms',
    ]
    
    for ext in extra_exts:
        if ext not in config.fn_clean_exts:
            config.fn_clean_exts.append(ext)
    
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

def execution_start():
    """ 
    Code to run when MultiQC starts (after config is loaded).
    Register search patterns for our custom modules.
    """
    # Register search patterns
    search_patterns = {
        'rnaseq/rnaseqqc': {
            'fn': '*metrics.tsv',
        },
        'rnaseq/gene_type_counts': {
            'fn': '*.json',
            'contents': 'gene_type_count'
        },
        'rnaseq/rsem': {
            'fn': ['*.json', '*.cnt'],
            'contents': 'num_genes_detected'
        },
    }
    config.update_dict(config.sp, search_patterns)

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