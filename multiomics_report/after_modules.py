"""
After modules hook and helper functions for MultiQC plugin.
Handles multiomics_report_after_modules hook and all merge logic.
"""
import logging
from typing import Dict
from collections import OrderedDict
from multiqc import report
from multiqc.types import ColumnKey, SectionKey
from .config import MAD_METRIC_KEYS, IS_PARENT_REGEX

log = logging.getLogger("multiqc")


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
    mad_headers = _extract_reproducibility_headers()

    # Phase 3: The Left Join
    _merge_parent_data_into_children(parent_data, mad_headers)

    _add_sambamba_calculated_metrics()


# ------------------------------------------------------------------------------
# HELPER FUNCTIONS
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
                
                log.debug(f"Plugin: Extracted parent '{sample_name}' from section '{section_key}' with keys: {list(parent_data[sample_name].keys())}")
                
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
        # Debug: show all parent data keys
        for parent_name, parent_stats in parent_data.items():
            log.debug(f"Plugin: Parent '{parent_name}' has metrics: {list(parent_stats.keys())}")
    
    return parent_data


def _extract_reproducibility_headers() -> OrderedDict:
    """Finds reproducibility headers in multiomics section, removes them there, returns them."""
    mad_headers = OrderedDict()
    
    # Find ALL multiomics sections (not just the first one)
    multiomics_section_keys = [
        k for k in report.general_stats_headers 
        if 'multiomic' in str(k).lower()
    ]
    
    log.debug(f"Plugin: Found {len(multiomics_section_keys)} multiomics sections: {multiomics_section_keys}")

    # Search through all multiomics sections
    for multiomics_section_key in multiomics_section_keys:
        headers_dict = report.general_stats_headers[multiomics_section_key]
        log.debug(f"Plugin: Checking section '{multiomics_section_key}' with headers: {list(headers_dict.keys())}")
        
        # Copy matching headers
        for key_str in MAD_METRIC_KEYS:
            # Skip if we already have this header
            col_key = ColumnKey(key_str)
            if col_key in mad_headers:
                continue
                
            # Check for Key object or String matches
            if col_key in headers_dict:
                mad_headers[col_key] = headers_dict[col_key].copy()
                del headers_dict[col_key] # Remove from source
                log.debug(f"Plugin: Extracted header '{key_str}' (as ColumnKey) from '{multiomics_section_key}'")
            elif key_str in headers_dict:
                mad_headers[ColumnKey(key_str)] = headers_dict[key_str].copy()
                del headers_dict[key_str] # Remove from source
                log.debug(f"Plugin: Extracted header '{key_str}' (as string) from '{multiomics_section_key}'")

        # Also scrub the data from the multiomics section to be clean
        if multiomics_section_key in report.general_stats_data:
            for sample_name, rows in report.general_stats_data[multiomics_section_key].items():
                for row in rows:
                    if hasattr(row, 'data') and row.data:
                        for k in list(row.data.keys()):
                            if str(k) in MAD_METRIC_KEYS:
                                log.debug(f"Plugin: Removing '{k}' from multiomics data for sample '{sample_name}'")
                                del row.data[k]

    # Fallback: Create defaults for missing headers
    log.debug(f"Plugin: Looking for MAD_METRIC_KEYS: {MAD_METRIC_KEYS}")
    log.debug(f"Plugin: Currently have headers: {list(mad_headers.keys())}")
    
    # Add missing MAD headers
    if ColumnKey('MAD of log ratios') not in mad_headers:
        mad_headers[ColumnKey('MAD of log ratios')] = {
            'title': 'MAD',
            'description': 'MAD QC: Median Absolute Deviation of log ratios',
            'format': '{:.4f}',
            'scale': 'RdYlGn',
            'namespace': 'multiomics'
        }
    if ColumnKey('Pearson correlation') not in mad_headers:
        mad_headers[ColumnKey('Pearson correlation')] = {
            'title': 'Pearson',
            'description': 'MAD QC: Pearson correlation coefficient',
            'min': 0, 'max': 1, 'format': '{:.4f}', 'scale': 'RdYlGn', 'namespace': 'multiomics'
        }
    if ColumnKey('Spearman correlation') not in mad_headers:
        mad_headers[ColumnKey('Spearman correlation')] = {
            'title': 'Spearman',
            'description': 'MAD QC: Spearman correlation coefficient',
            'min': 0, 'max': 1, 'format': '{:.4f}', 'scale': 'RdYlGn', 'namespace': 'multiomics'
        }
    
    # Add jaccard header if missing (this is the key fix!)
    if ColumnKey('jaccard') not in mad_headers:
        log.debug("Plugin: Creating default jaccard header")
        mad_headers[ColumnKey('jaccard')] = {
            'title': 'Jaccard',
            'description': 'Jaccard similarity coefficient',
            'min': 0,
            'max': 1,
            'format': '{:.4f}',
            'scale': 'RdYlGn',
            'namespace': 'multiomics'
        }
    
    # Optional: Add other jaccard-related headers if needed
    # (intersection, union, n_intersections) - but you removed these from MAD_METRIC_KEYS
    
    log.debug(f"Plugin: Final headers to merge: {list(mad_headers.keys())}")
    return mad_headers


def _merge_parent_data_into_children(parent_data: Dict, headers: OrderedDict):
    """Iterates all child samples and merges parent data if a match is found."""
    
    match_count = 0
    # Track where we've added headers so we only do it once per section
    sections_with_headers_added = set()
    
    log.debug(f"Plugin: Starting merge with {len(parent_data)} parent samples")
    log.debug(f"Plugin: Available headers to merge: {list(headers.keys())}")

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
                log.debug(f"Plugin: Found match: child '{sample_name}' -> parent '{base_name}' in section '{section_key}'")
                log.debug(f"Plugin: Parent '{base_name}' has keys: {list(source_stats.keys())}")
                
                # Check for merge collision (Child already has this data?)
                # We check the first row for any MAD metric key
                has_mad_data = False
                existing_keys = set()
                for row in rows:
                    if hasattr(row, 'data') and row.data:
                        existing_keys.update(row.data.keys())
                        for mad_key in MAD_METRIC_KEYS:
                            if ColumnKey(mad_key) in row.data or mad_key in row.data:
                                has_mad_data = True
                                break
                        if has_mad_data:
                            break
                
                if has_mad_data:
                    log.debug(f"Plugin: Child '{sample_name}' already has MAD data, skipping merge")
                    continue # Skip, already merged

                # Perform Merge
                merged_keys = []
                for row in rows:
                    if hasattr(row, 'data') and row.data:
                        # Update child row with parent data (only new keys)
                        for k, v in source_stats.items():
                            if k not in row.data:
                                row.data[k] = v
                                merged_keys.append(str(k))
                
                log.debug(f"Plugin: Merged keys into '{sample_name}': {merged_keys}")
                match_count += 1
                section_modified = True

        # If we merged data in this section, we must add the headers
        if section_modified and section_key not in sections_with_headers_added:
            if section_key in report.general_stats_headers:
                # Merge dictionaries
                report.general_stats_headers[section_key].update(headers)
                sections_with_headers_added.add(section_key)
                log.debug(f"Plugin: Added headers to section '{section_key}'")

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

