"""
Search pattern registration for MultiQC plugin.
Handles execution_start hook.
"""
import logging
from multiqc import config
from multiqc.utils.util_functions import update_dict

log = logging.getLogger("multiqc")


def execution_start():
    """Register search patterns."""
    search_patterns = {
        'multiomics_report/rnaseqqc': {
            'fn': '*metrics.tsv', 
            'contents': 'Duplicate Rate of Mapped'
        },
        'multiomics_report/gene_type_counts': {
            'fn': '*gene_type_count.json',
            'contents': 'gene_type_count',
            'num_lines': 50
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
        'multiomics_report/peak_count': {'fn': "*_peakcount.txt"},
        'multiomics_report/jaccard': {'fn': "*_jaccard.txt"},
        'multiomics_report/fraglen': {'fn': "*_fraglen.txt"},
        'multiomics_report/frip': {'fn': "frip*.tsv"},
        'multiomics_report/preseq' : {'fn': "lcextrap*.txt"},
        'multiomics_report/bam_correlation': {'fn': "*_bam_correlation_stats_mqc.tsv"},
        'multiomics_report/frip_atac': {'fn': "*stats.tsv", 'contents': 'frip'},
        'multiomics_report/tss': {'fn': "*tss_histogram.csv"},
        'multiomics_report/prealign': {'fn': "*prealign.stats.tsv", 'contents': 'percent_filtered'},
        'multiomics_report/samblaster': {'fn': "*.samblaster.log", 'contents': 'Marked'},  
        'multiomics_report/replicate_correlations': {'fn': "*replicate_correlations.tsv"},
        'multiomics_report/reproducibility_qc': {'fn': "*reproducibility.qc.json", 'contents': 'reproducibility'},
        'multiomics_report/hic_mapstat' :{'fn': "*.mapstat", 'contents': 'total'},
        'multiomics_report/hic_pairstat' :{'fn': "*.pairstat"},
        'multiomics_report/RSstat' :{'fn': "*.RSstat"},
        'multiomics_report/dedup_stats' :{'fn': "*.dedup.stats"},
    }
    
    # Modify custom_content search pattern to exclude correlation files
    if 'custom_content' in config.sp:
        custom_content_pattern = config.sp['custom_content']
        
        # Ensure it's a dict (not a list)
        if isinstance(custom_content_pattern, dict):
            # Add exclude pattern
            if 'exclude_fn' not in custom_content_pattern:
                custom_content_pattern['exclude_fn'] = []
            elif not isinstance(custom_content_pattern['exclude_fn'], list):
                custom_content_pattern['exclude_fn'] = [custom_content_pattern['exclude_fn']]
            
            # Add our exclude patterns
            exclude_patterns = [
                "*bam_correlation*.yaml",
                "*bam_correlation*.yml",
            ]
            for excl in exclude_patterns:
                if excl not in custom_content_pattern['exclude_fn']:
                    custom_content_pattern['exclude_fn'].append(excl)
            
            log.info(f"Plugin: Added exclude patterns to custom_content: {exclude_patterns}")
    
    config.sp = update_dict(config.sp, search_patterns, add_in_the_beginning=True)

