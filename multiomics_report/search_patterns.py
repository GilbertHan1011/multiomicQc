"""
Search pattern registration for MultiQC plugin.
Handles execution_start hook.
"""
from multiqc import config
from multiqc.utils.util_functions import update_dict


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
    }
    config.sp = update_dict(config.sp, search_patterns, add_in_the_beginning=True)

