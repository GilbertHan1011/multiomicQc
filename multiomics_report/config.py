"""
Configuration and setup hooks for MultiQC plugin.
Handles before_config and config_loaded hooks.
"""
import logging
import re
from multiqc import config
from pkg_resources import get_distribution

log = logging.getLogger("multiqc")

# --- CONSTANTS & CONFIGURATION ---

# 1. Metrics to move/merge (The "MAD" headers and related QC metrics)
# These keys are propagated from parent samples to child samples in
# `after_modules._merge_parent_data_into_children()`.
MAD_METRIC_KEYS = [
    # MAD QC metrics
    'MAD of log ratios',
    'Pearson correlation',
    'Spearman correlation',
    # Peak correlation p-value (from replicate_correlations.tsv)
    'pearson_p',
    # Reproducibility QC metrics (from *reproducibility.qc.json)
    'rescue_ratio',
    'self_consistency_ratio',
    'N_optimal',
    # Jaccard metrics (from *_jaccard.txt)
    'jaccard',
    # BAM Correlation metrics (from *_bam_correlation_stats_mqc.tsv)
    'BAM_Correlation',
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
        "_coverage.tsv", "_peakcount.txt",
        "_jaccard.txt", "_bam_correlation_stats_mqc.tsv",
        ".prealign.stats.tsv",".stats.tsv",
        ".tss_histogram.csv", ".replicate_correlations.tsv",
        ".mapstat", ".pairstat", "hicpro.RSstat", 
        '.dedup.stats', '_loglog_fits.csv', 
        '.loop_counts.tsv', 'complexity',
        {'type': 'remove', 'pattern': 'sambamba_markdup_'},  # Remove prefix
        {'type': 'remove', 'pattern': 'bowtie2_'},            # Remove prefix
        {'type': 'remove', 'pattern': '.err'},                 # Remove suffix
        {'type': 'remove', 'pattern': 'frip_'}, 
        {'type': 'remove', 'pattern': 'lcextrap_'},
        {'type': 'remove', 'pattern': r'\.mapq.*'},  # Remove .mapq and anything after it (e.g., .mapq_30.1000) 
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
        '': [
            {'type': 'regex', 'pattern': r'(_run|_)\d+$'},
            {'type': 'regex', 'pattern': r'(_R|_)\d+$'}
        ]
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

