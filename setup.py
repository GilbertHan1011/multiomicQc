#!/usr/bin/env python3
"""
MultiQC plugin for reporting results of the RNA-seq pipeline.
"""

from setuptools import setup, find_packages

setup(
    name='multiomics_report',
    version='0.1.0',
    description='Custom multiomics plugin for MultiQC',
    packages=find_packages(),
    include_package_data=True,
    url = 'https://github.com/GilbertHan1011/multiomicQc',
    download_url = 'https://github.com/GilbertHan1011/multiomicQc',
    install_requires=[
        'multiqc',
        'click',
        'tables'
    ],
    entry_points={
        'multiqc.modules.v1': [
            # 'name = package.module:Class'
            'multiomics_report = multiomics_report.modules.multiomics:MultiqcModule',
        ],
        'multiqc.hooks.v1': [
            # Hook before_config to configure sample renaming (runs before file search)
            'before_config = multiomics_report:before_config',
            'config_loaded = multiomics_report:config_loaded',
            # Hook execution_start to load search patterns
            'execution_start = multiomics_report:execution_start',
            'after_modules = multiomics_report:multiomics_report_after_modules',
        ]
    },
    classifiers=[
        'Development Status :: 4 - Beta',
        'Environment :: Console',
        'Environment :: Web Environment',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: MIT License',
        'Natural Language :: English',
        'Operating System :: MacOS :: MacOS X',
        'Operating System :: POSIX',
        'Operating System :: Unix',
        'Programming Language :: Python',
        'Programming Language :: JavaScript',
        'Topic :: Scientific/Engineering',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Topic :: Scientific/Engineering :: Visualization',
    ]
)