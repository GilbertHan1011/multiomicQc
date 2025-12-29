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
# Parents do NOT have a hyphen
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
        '.metrics', '.isoforms', 
        '.genes', '_fastp',"_gene_type_count",
        "_coverage.tsv", "_peakcount.txt", "_jaccard.txt",
        {'type': 'remove', 'pattern': 'sambamba_markdup_'},  # Remove prefix
        {'type': 'remove', 'pattern': 'bowtie2_'},            # Remove prefix
        {'type': 'remove', 'pattern': '.err'},                 # Remove suffix
        {'type': 'remove', 'pattern': 'frip_'}, 
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
        'multiomics_report/rnaseqqc': {
            'fn': '*metrics.tsv', 
            'contents': 'Duplicate Rate of Mapped'
            },
        'multiomics_report/gene_type_counts': {
            'fn': '*gene_type_count.json',  # More specific filename pattern
            'contents': 'gene_type_count',
            'num_lines': 50  # Check more lines for JSON files that might be formatted
        },
        "multiomics_report/rsem": [
            {'fn': "*.json", 'contents': 'num_genes_detected'},
            {'fn': "*.cnt", 'contents': 'num_genes_detected'},
        ],
        'multiomics_report/mad_rna': {
            'fn_re': r'.*summary_metrics\.json$',
            'contents': 'MAD of log ratios',
            'num_lines': 10
        },
        'multiomics_report/peak_coverage': {'fn': "*_coverage.tsv"},
        'multiomics_report/peak_count' : {'fn': "*_peakcount.txt"},
        'multiomics_report/jacarrd' : {'fn' : "*_jaccard.txt"},
        'multiomics_report/fraglen' : {'fn': "*_fraglen.txt"},
        'multiomics_report/frip' : {'fn': "frip*.tsv"},
    }
    config.sp = update_dict(config.sp, search_patterns, add_in_the_beginning=True)


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

    _add_sambamba_calculated_metrics()


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



def _add_sambamba_calculated_metrics():
    """
    Add calculated metrics from sambamba markdup data:
    - nrf_sm: Non-Redundant Fraction = 100 - duplicate_rate
    - duplicate_read: Number of duplicate reads
    - unique_read: (sorted_end_pairs * 2 - duplicate_reads) / 2
    - total_reads: Number of paired-end fragments seen by Sambamba
    """
    # Find Sambamba module section in general_stats_data
    sambamba_section_key = None
    for section_key in report.general_stats_data.keys():
        section_str = str(section_key).lower()
        if 'sambamba' in section_str or 'markdup' in section_str:
            sambamba_section_key = section_key
            break
    
    if not sambamba_section_key:
        log.debug("Plugin: Sambamba module not found in general stats")
        return
    
    log.info(f"Plugin: Found Sambamba section: {sambamba_section_key}")
    
    # Process each sample in the sambamba section
    samples_dict = report.general_stats_data.get(sambamba_section_key, {})
    modified_count = 0
    
    for sample_group, rows in samples_dict.items():
        sample_name = str(sample_group)
        
        # Extract data from the first row (usually there's only one row per sample)
        sample_data = {}
        for row in rows:
            if hasattr(row, 'data') and row.data:
                sample_data.update(row.data)
        
        if not sample_data:
            continue
        
        # Calculate nrf_sm = 100 - duplicate_rate
        if "duplicate_rate" in sample_data:
            try:
                duplicate_rate = float(sample_data["duplicate_rate"])
                nrf_sm = 100.0 - duplicate_rate
                # Add to all rows for this sample
                for row in rows:
                    if hasattr(row, 'data') and row.data:
                        row.data["nrf_sm"] = nrf_sm
                log.debug(f"Plugin: Added nrf_sm for {sample_name}: {nrf_sm:.2f}%")
            except (ValueError, TypeError) as e:
                log.debug(f"Plugin: Could not calculate nrf_sm for {sample_name}: {e}")
        
        # Add duplicate_reads as duplicate_read (if not already present)
        if "duplicate_reads" in sample_data and "duplicate_read" not in sample_data:
            for row in rows:
                if hasattr(row, 'data') and row.data:
                    row.data["duplicate_read"] = sample_data["duplicate_reads"]
            log.debug(f"Plugin: Added duplicate_read for {sample_name}: {sample_data['duplicate_reads']}")
        
        # Calculate unique_read
        if "sorted_end_pairs" in sample_data and "duplicate_reads" in sample_data:
            try:
                sorted_end_pairs = int(sample_data["sorted_end_pairs"])
                duplicate_reads = int(sample_data["duplicate_reads"])
                unique_read = (sorted_end_pairs * 2 - duplicate_reads) / 2
                for row in rows:
                    if hasattr(row, 'data') and row.data:
                        row.data["unique_read"] = unique_read
                log.debug(f"Plugin: Added unique_read for {sample_name}: {unique_read:.0f}")
            except (ValueError, TypeError) as e:
                log.debug(f"Plugin: Could not calculate unique_read for {sample_name}: {e}")
        
        # Add total_reads (from sorted_end_pairs)
        if "sorted_end_pairs" in sample_data:
            try:
                total_reads = int(sample_data["sorted_end_pairs"])
                for row in rows:
                    if hasattr(row, 'data') and row.data:
                        row.data["total_reads"] = total_reads
                log.debug(f"Plugin: Added total_reads for {sample_name}: {total_reads}")
            except (ValueError, TypeError):
                pass
        
        modified_count += 1
    
    # Add headers for new metrics
    if modified_count > 0:
        _add_sambamba_headers(sambamba_section_key)
        log.info(f"Plugin: Successfully added calculated sambamba metrics for {modified_count} samples")


def _add_sambamba_headers(section_key: SectionKey):
    """Add headers for calculated sambamba metrics."""
    if section_key not in report.general_stats_headers:
        report.general_stats_headers[section_key] = {}
    
    headers = report.general_stats_headers[section_key]
    
    # nrf_sm header
    if ColumnKey('nrf_sm') not in headers:
        headers[ColumnKey('nrf_sm')] = {
            'title': 'NRF',
            'description': 'Sambamba: Non-Redundant Fraction (100 - duplicate_rate)',
            'min': 0,
            'max': 100,
            'suffix': '%',
            'format': '{:,.2f}',
            'scale': 'RdYlGn',
            'namespace': 'sambamba'
        }
    
    # duplicate_read header
    if ColumnKey('duplicate_read') not in headers:
        headers[ColumnKey('duplicate_read')] = {
            'title': 'Dup Reads',
            'description': 'Sambamba: Number of duplicate reads',
            'format': '{:,.0f}',
            'scale': 'OrRd',
            'namespace': 'sambamba'
        }
    
    # unique_read header
    if ColumnKey('unique_read') not in headers:
        headers[ColumnKey('unique_read')] = {
            'title': 'Unique Reads',
            'description': 'Sambamba: Unique reads ((sorted_end_pairs * 2 - duplicate_reads) / 2)',
            'format': '{:,.0f}',
            'scale': 'GnBu',
            'namespace': 'sambamba'
        }
    
    # total_reads header
    if ColumnKey('total_reads') not in headers:
        headers[ColumnKey('total_reads')] = {
            'title': 'Total Reads',
            'description': 'Sambamba: Total paired-end fragments (sorted_end_pairs)',
            'format': '{:,.0f}',
            'scale': 'Blues',
            'namespace': 'sambamba'
        }