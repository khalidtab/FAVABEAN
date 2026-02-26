"""
Standalone SIDLE utilities
From q2-sidle: https://github.com/jwdebelius/q2-sidle
Copyright (c) 2020, Justine Debelius. BSD-3-Clause License
"""
import numpy as np
import pandas as pd
from Bio import SeqIO

# Degenerate nucleotide mappings
degenerate_map = {
    "R": ['A', 'G'], 'Y': ['C', 'T'], 'S': ['G', 'C'], 'W': ['A', 'T'],
    'K': ['G', 'T'], 'M': ['A', 'C'], 'B': ['C', 'G', 'T'], 'D': ['A', 'G', 'T'],
    'H': ['A', 'C', 'T'], 'V': ['A', 'C', 'G'], 'N': ['A', 'C', 'G', 'T'],
}

degen_sub2 = {
    'R': '[AG]', 'Y': '[CT]', 'S': '[CG]', 'W': '[AT]', 'K': '[GT]', 'M': '[AC]',
    'B': '[CGT]', 'D': '[AGT]', 'H': '[ACT]', 'V': '[ACG]', 'N': '[ACGT]'
}

database_params = {
    # All delimiters use ';' (no space) because FASTA headers are normalized
    # by read_fasta() and the extraction rules ('; ' → ';').
    'greengenes': {'delim': ';', 'defined': lambda x: len(x) > 3},
    'silva': {'delim': ';', 'defined': lambda x: not (('uncul' in x) or ('metagenome' in x))},
    'ehomd': {
        'delim': ';',
        'defined': lambda x: not (('unclassified' in x.lower()) or
                                  ('uncultured' in x.lower()) or
                                  len(x.strip()) <= 3)
    },
    'none': {'delim': ';', 'defined': lambda x: True},
}

def read_fasta(filepath):
    """Read FASTA into a {seq_id: sequence} Series.

    Uses rec.description (the full header line) and normalizes spaces
    after semicolons so that taxonomy-as-ID formats like eHOMD
    ('; species') or Greengenes ('; p__') are not truncated by
    BioPython's rec.id whitespace split.
    """
    result = {}
    for rec in SeqIO.parse(filepath, "fasta"):
        seq_id = rec.description.replace('; ', ';')
        result[seq_id] = str(rec.seq).upper()
    return pd.Series(result)

def write_fasta(sequences, filepath):
    """Write sequences to FASTA, normalizing headers ('; ' → ';')."""
    with open(filepath, 'w') as f:
        items = sequences.items() if isinstance(sequences, (pd.Series, dict)) else sequences
        for seq_id, seq in items:
            norm_id = str(seq_id).replace('; ', ';')
            f.write(f">{norm_id}\n{seq}\n")

def _check_regions(region):
    region, region_idx = np.unique(region, return_index=True)
    region_order = {r: i for i, r in zip(region_idx, region)}
    region_names = {i: r for r, i in region_order.items()}
    return region_order, region_names, len(region_order)
