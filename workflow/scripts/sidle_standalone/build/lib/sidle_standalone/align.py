"""
Standalone SIDLE k-mer alignment
From q2-sidle: https://github.com/jwdebelius/q2-sidle
Copyright (c) 2020, Justine Debelius. BSD-3-Clause License
"""
import warnings
warnings.filterwarnings('ignore', category=RuntimeWarning)

import numpy as np
import pandas as pd

from .utils import read_fasta


def align_regional_kmers(kmers_file, rep_seq_file, region, output_file,
                         max_mismatch=2, **kwargs):
    """
    Align ASV sequences to regional k-mer database using vectorized numpy.

    Parameters
    ----------
    kmers_file : str
        Path to FASTA file with reference k-mers
    rep_seq_file : str
        Path to FASTA file with ASV representative sequences
    region : str
        Region identifier
    output_file : str
        Path for output alignment TSV
    max_mismatch : int
        Maximum mismatches allowed (default: 2)

    Returns
    -------
    pd.DataFrame
        Alignment results with columns: kmer, asv, length, mismatch, region, max-mismatch
    """
    # Load sequences
    rep_seq = read_fasta(rep_seq_file)
    kmers = read_fasta(kmers_file)

    # Check consistent lengths
    asv_lengths = rep_seq.apply(len).unique()
    kmer_lengths = kmers.apply(len).unique()

    if len(asv_lengths) > 1:
        raise ValueError(f'ASV sequences have inconsistent lengths: {asv_lengths}. '
                         f'Run trim_dada2_posthoc() first to trim ASVs to uniform length.')
    if len(kmer_lengths) > 1:
        raise ValueError(f'K-mer sequences have inconsistent lengths: {kmer_lengths}')
    if asv_lengths[0] != kmer_lengths[0]:
        raise ValueError(f'K-mer length ({kmer_lengths[0]}) != ASV length ({asv_lengths[0]}). '
                         f'Trim ASVs to match k-mer length using trim_dada2_posthoc().')

    seq_length = asv_lengths[0]

    # Vectorized alignment via numpy broadcasting
    # Convert sequences to 2D character arrays: (n_seqs, seq_length)
    kmer_ids = kmers.index.values
    asv_ids = rep_seq.index.values

    kmer_arr = np.array([list(s) for s in kmers.values])   # (K, L)
    asv_arr = np.array([list(s) for s in rep_seq.values])   # (A, L)

    # Broadcast comparison: (K, 1, L) != (1, A, L) -> (K, A, L)
    # Then sum mismatches along L -> (K, A)
    mismatch_matrix = (kmer_arr[:, np.newaxis, :] != asv_arr[np.newaxis, :, :]).sum(axis=2)

    # Find pairs within mismatch threshold
    kmer_idx, asv_idx = np.where(mismatch_matrix <= max_mismatch)

    result = pd.DataFrame({
        'kmer': kmer_ids[kmer_idx],
        'asv': asv_ids[asv_idx],
        'length': seq_length,
        'mismatch': mismatch_matrix[kmer_idx, asv_idx],
        'region': region,
        'max-mismatch': max_mismatch,
    })

    # Save and return
    result.to_csv(output_file, sep='\t', index=False)
    print(f"Aligned {len(kmers)} k-mers x {len(rep_seq)} ASVs -> {len(result)} hits "
          f"(max {max_mismatch} mismatches)")
    return result
