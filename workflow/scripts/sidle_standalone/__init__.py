"""
Standalone SIDLE - Multi-region 16S reconstruction without QIIME2

From q2-sidle: https://github.com/jwdebelius/q2-sidle
Copyright (c) 2020, Justine Debelius. BSD-3-Clause License

Modified to remove QIIME2 dependencies while preserving core functionality.
"""

from .extract import prepare_extracted_region, trim_dada2_posthoc
from .align import align_regional_kmers
from .reconstruct import reconstruct_counts
from .build_database import reconstruct_database
from .taxonomy import reconstruct_taxonomy

__version__ = '1.0.0-standalone'
__all__ = [
    'prepare_extracted_region',
    'trim_dada2_posthoc',
    'align_regional_kmers',
    'reconstruct_database',
    'reconstruct_counts',
    'reconstruct_taxonomy'
]
