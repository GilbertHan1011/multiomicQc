import os
import csv
import logging
import re
from typing import Dict, List, Set, Optional, Any
from collections import OrderedDict
from multiqc import config, report
from multiqc.types import ColumnKey, SectionKey
from pkg_resources import get_distribution

log = logging.getLogger("multiqc")

# --- CONSTANTS & CONFIGURATION ---

# 1. Metrics to move/merge (The "MAD" headers)
MAD_METRIC_KEYS = [
    'MAD of log ratios', 
    'Pearson correlation', 
    'Spearman correlation', 
    'SD of log ratios', 
    'num_pairs_evaluated'
]

# 2. Regex to identify Parent Samples (Group B) vs Child Samples (Group A)
# Parents do NOT have a hyphen (e.g., 'test', 'apple')
IS_PARENT_REGEX = re.compile(r'^[^-]+$')

# 3. Defaults for Config
MODULE_DEFAULTS = {
    'fastp': {'s_name_filenames': True},
}

# Save this plugin's version number (defined in setup.py) to the MultiQC config
try:
    config.multiomics_report_version = get_distribution("multiomics_report").version
except Exception:
    config.multiomics_report_version = "unknown"


def before_config():
    """Set up configuration defaults and cleaning rules."""
    # 1. Clean Extensions
    if not getattr(config, 'fn_clean_sample_names', True):
        log.warning("Plugin: fn_clean_sample_names was False. Enabling it.")
        config.fn_clean_sample_names = True

    if not hasattr(config, 'fn_clean_exts'):
        config.fn_clean_exts = []

    extra_exts = [
        '_Log.final.out', '.summary_metrics.json', '.metrics.tsv', 
        '.metrics', '.isoforms', '.genes', '_fastp'
    ]
    # Prepend to ensure higher priority
    config.fn_clean_exts[0:0] = extra_exts
    log.debug(f"Plugin: Added {len(extra_exts)} patterns to fn_clean_exts")

    # 2. Module Defaults
    for mod, defaults in MODULE_DEFAULTS.items():
        if not hasattr(config, mod):
            config.update({mod: defaults})
        else:
            # Shallow merge: user config wins
            existing = getattr(config, mod)
            if isinstance(existing, dict):
                config.update({mod: {**defaults, **existing}})

    # 3. Table Sample Merge Rules
    if not hasattr(config, 'table_sample_merge') or config.table_sample_merge is None:
        config.table_sample_merge = {}
    
    merge_rule = {
        '': [{'type': 'regex', 'pattern': r'(_run|_)\d+$'}]
    }
    config.table_sample_merge.update(merge_rule)


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
    """Register search patterns."""
    from multiqc.utils.util_functions import update_dict
    
    search_patterns = {
        'multiomics_report/rnaseqqc': {'fn': '*metrics.tsv'},
        'multiomics_report/gene_type_counts': {'fn': '*.json', 'contents': 'gene_type_count'},
        "multiomics_report/rsem": [
            {'fn': "*.json", 'contents': 'num_genes_detected'},
            {'fn': "*.cnt", 'contents': 'num_genes_detected'},
        ],
        'multiomics_report/mad_rna': {
            'fn_re': r'.*summary_metrics\.json$',
            'contents': 'MAD of log ratios',
            'num_lines': 10
        }
    }
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
                source_name = f"{target_name}_run{run_id}"
                
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
    Main Logic Hook:
    1. Harvests 'Parent' (B) data and removes it from the table.
    2. Harvests 'MAD' headers and removes them from the source section.
    3. Merges Parent data into 'Child' (A) samples.
    """
    
    # Phase 1: extract B samples (test, apple) and remove them from display
    parent_data = _extract_and_hide_parents()
    
    if not parent_data:
        log.info("Plugin: No parent samples found. Skipping merge.")
        return

    # Phase 2: Get headers for the metrics we want to transfer
    mad_headers = _extract_mad_headers()

    # Phase 3: The Left Join
    _merge_parent_data_into_children(parent_data, mad_headers)


# ------------------------------------------------------------------------------
# HELPER FUNCTIONS (The "Elegance" part)
# ------------------------------------------------------------------------------

def _extract_and_hide_parents() -> Dict[str, Dict]:
    """Finds samples matching the Parent Regex, extracts data, and deletes them from report."""
    parent_data = {}
    removed_count = 0

    # Iterate over a list(keys) because we will modify the dictionary size during iteration
    for section_key, samples_dict in list(report.general_stats_data.items()):
        for sample_group in list(samples_dict.keys()):
            sample_name = str(sample_group)

            # Check if this is a parent (e.g., no hyphens)
            if IS_PARENT_REGEX.match(sample_name):
                
                # Init storage
                if sample_name not in parent_data:
                    parent_data[sample_name] = {}

                # Harvest data from all rows
                for row in samples_dict[sample_group]:
                    if hasattr(row, 'data') and row.data:
                        parent_data[sample_name].update(row.data)
                
                # Remove from report
                del samples_dict[sample_group]
                removed_count += 1

        # Clean up empty sections
        if not samples_dict:
            del report.general_stats_data[section_key]
            if section_key in report.general_stats_headers:
                del report.general_stats_headers[section_key]

    if removed_count > 0:
        log.info(f"Plugin: Extracted and hid {removed_count} parent sample entries.")
    
    return parent_data


def _extract_mad_headers() -> OrderedDict:
    """Finds MAD headers in multiomics section, removes them there, returns them."""
    mad_headers = OrderedDict()
    
    # Find the source section
    multiomics_section_key = next(
        (k for k in report.general_stats_headers if 'multiomic' in str(k).lower()), 
        None
    )

    if multiomics_section_key:
        headers_dict = report.general_stats_headers[multiomics_section_key]
        
        # Copy matching headers
        for key_str in MAD_METRIC_KEYS:
            # Check for Key object or String matches
            col_key = ColumnKey(key_str)
            if col_key in headers_dict:
                mad_headers[col_key] = headers_dict[col_key].copy()
                del headers_dict[col_key] # Remove from source
            elif key_str in headers_dict:
                mad_headers[ColumnKey(key_str)] = headers_dict[key_str].copy()
                del headers_dict[key_str] # Remove from source

        # Also scrub the data from the multiomics section to be clean
        if multiomics_section_key in report.general_stats_data:
            for _, rows in report.general_stats_data[multiomics_section_key].items():
                for row in rows:
                    if hasattr(row, 'data') and row.data:
                        for k in list(row.data.keys()):
                            if str(k) in MAD_METRIC_KEYS:
                                del row.data[k]

    # Fallback: Create defaults if missing
    if not mad_headers:
        log.warning("Plugin: MAD headers not found. Using defaults.")
        mad_headers[ColumnKey('MAD of log ratios')] = {
            'title': 'MAD',
            'description': 'MAD QC: Median Absolute Deviation of log ratios',
            'format': '{:.4f}',
            'scale': 'RdYlGn',
            'namespace': 'multiomics'
        }
        mad_headers[ColumnKey('Pearson correlation')] = {
            'title': 'Pearson',
            'description': 'MAD QC: Pearson correlation coefficient',
            'min': 0, 'max': 1, 'format': '{:.4f}', 'scale': 'RdYlGn', 'namespace': 'multiomics'
        }
        mad_headers[ColumnKey('Spearman correlation')] = {
            'title': 'Spearman',
            'description': 'MAD QC: Spearman correlation coefficient',
            'min': 0, 'max': 1, 'format': '{:.4f}', 'scale': 'RdYlGn', 'namespace': 'multiomics'
        }
        mad_headers[ColumnKey('SD of log ratios')] = {
            'title': 'SD',
            'description': 'MAD QC: Standard Deviation of log ratios',
            'format': '{:.4f}', 'scale': 'RdYlGn', 'namespace': 'multiomics'
        }
        mad_headers[ColumnKey('num_pairs_evaluated')] = {
            'title': 'Pairs',
            'description': 'MAD QC: Number of pairs evaluated',
            'format': '{:,.0f}', 'scale': 'Blues', 'namespace': 'multiomics'
        }
            
    return mad_headers


def _merge_parent_data_into_children(parent_data: Dict, headers: OrderedDict):
    """Iterates all child samples and merges parent data if a match is found."""
    
    match_count = 0
    # Track where we've added headers so we only do it once per section
    sections_with_headers_added = set()

    for section_key, samples_dict in report.general_stats_data.items():
        
        # Skip multiomics section (source)
        if 'multiomic' in str(section_key).lower():
            continue

        section_modified = False

        for sample_group, rows in samples_dict.items():
            sample_name = str(sample_group)
            
            # Logic: Child 'test1-1' -> Parent 'test1'
            # Adjust this split logic if your naming convention changes!
            base_name = sample_name.split('-')[0]

            if base_name in parent_data:
                source_stats = parent_data[base_name]
                
                # Check for merge collision (Child already has this data?)
                # We check the first row for any MAD metric key
                has_mad_data = False
                for row in rows:
                    if hasattr(row, 'data') and row.data:
                        for mad_key in MAD_METRIC_KEYS:
                            if ColumnKey(mad_key) in row.data or mad_key in row.data:
                                has_mad_data = True
                                break
                        if has_mad_data:
                            break
                
                if has_mad_data:
                    continue # Skip, already merged

                # Perform Merge
                for row in rows:
                    if hasattr(row, 'data') and row.data:
                        # Update child row with parent data (only new keys)
                        for k, v in source_stats.items():
                            if k not in row.data:
                                row.data[k] = v
                
                match_count += 1
                section_modified = True

        # If we merged data in this section, we must add the headers
        if section_modified and section_key not in sections_with_headers_added:
            if section_key in report.general_stats_headers:
                # Merge dictionaries
                report.general_stats_headers[section_key].update(headers)
                sections_with_headers_added.add(section_key)

    log.info(f"Plugin: Left-joined stats into {match_count} child samples.")