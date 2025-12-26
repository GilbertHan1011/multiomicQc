#!/usr/bin/env python
""" MultiQC BSF plugin functions
We can add any custom Python functions here and call them
using the setuptools plugin hooks.
"""


from __future__ import print_function
from pkg_resources import get_distribution
import logging
import os, csv

from multiqc.utils import report, util_functions, config

# Initialise the main MultiQC logger
log = logging.getLogger('multiqc')

# Save this plugin's version number (defined in setup.py) to the MultiQC config
config.rnaseq_report_version = get_distribution("rnaseq_report").version


# Add default config options for the things that are used in atacseq_report
def rnaseq_report_execution_start():
    """
    Code to execute after the config files and
    command line flags have been parsed self.
    this setuptools hook is the earliest that will be able
    to use custom command line flags.
    """
    # Halt execution if we've disabled the plugin
    if config.kwargs.get('disable_rnaseq_report', True):
        return None

    log.info("Running rnaseq_report MultiQC Plugin v{}, use --disable-rnaseq-report to disable".format(config.rnaseq_report_version))

    # Add to the search patterns used by atacseq module
    if 'rnaseq/rnaseqqc' not in config.sp:
        config.update_dict(config.sp, {'rnaseq/rnaseqqc': {'fn': '*.metrics.tsv', 'contents': 'Mapping Rate'}})
        log.info("updated config.sp for rnaseq/rnaseqqc")
    if 'rnaseq/gene_type_counts' not in config.sp:
        config.update_dict(config.sp, {'rnaseq/gene_type_counts': {'fn': '*.json', 'contents': 'gene_type_count'}})
        log.info("updated config.sp for rnaseq/gene_type_counts")
    if 'rnaseq/rsem' not in config.sp:
        config.update_dict(config.sp, {'rnaseq/rsem': {'fn': '*.json', 'contents': 'num_genes_detected'}})
        log.info("updated config.sp for rnaseq/rsem")


def rnaseq_report_after_modules():
    """
    Hook that runs after all modules are initialized.
    This allows us to access data from other modules like rnaseq.
    """
    # This hook runs after modules, but we'll handle data access in the module itself
    # using a different approach
    pass