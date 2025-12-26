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
config.multiomics_report_version = get_distribution("multiomics_report").version


def multiomics_report_after_modules():
    """
    Hook that runs after all modules are initialized.
    This allows us to access data from other modules like multiomics.
    """
    # This hook runs after modules, but we'll handle data access in the module itself
    # using a different approach
    pass