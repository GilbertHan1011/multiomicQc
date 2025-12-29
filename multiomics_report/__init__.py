"""
MultiQC plugin for multiomics reporting.
Main entry point that imports and exposes hooks from separate modules.
"""
from .config import before_config, config_loaded
from .search_patterns import execution_start
from .after_modules import multiomics_report_after_modules

# Re-export hooks for MultiQC to discover
__all__ = ['before_config', 'config_loaded', 'execution_start', 'multiomics_report_after_modules']