"""
Standalone SIDLE database extraction
From q2-sidle: https://github.com/jwdebelius/q2-sidle
Copyright (c) 2020, Justine Debelius. BSD-3-Clause License
"""
import numpy as np
import pandas as pd
from skbio import DNA

from .utils import read_fasta, write_fasta


def trim_dada2_posthoc(rep_seqs_file, count_file, output_fasta, output_counts, trim_length=0):
    """
    Trim DADA2 representative sequences to uniform length for SIDLE alignment.

    DADA2 merged paired-end reads produce variable-length ASVs, but SIDLE's
    k-mer alignment requires all ASVs to be the same length. This function
    trims sequences, discards those shorter than the threshold, and collapses
    any duplicates created by trimming (summing their counts).

    Mirrors q2-sidle's trim_dada2_posthoc but operates on FASTA + TSV files
    instead of QIIME2 Artifacts.

    Parameters
    ----------
    rep_seqs_file : str
        Path to FASTA file of DADA2 representative sequences
    count_file : str
        Path to count table TSV (rows = ASV IDs, columns = samples)
    output_fasta : str
        Path for output trimmed FASTA
    output_counts : str
        Path for output trimmed/collapsed count table TSV
    trim_length : int
        Length to trim ASVs to. If 0, uses the minimum sequence length
        in the input. If positive, trims from the left (5' end). If
        negative, trims from the right (3' end).

    Returns
    -------
    tuple
        (trimmed_seqs as pd.Series, trimmed_counts as pd.DataFrame)
    """
    import hashlib
    import warnings

    # Read representative sequences
    rep_seqs = read_fasta(rep_seqs_file)

    # Read count table
    counts = pd.read_csv(count_file, sep='\t', index_col=0)
    counts = counts.apply(pd.to_numeric, errors='coerce').fillna(0)

    # Determine trim length
    seq_lengths = rep_seqs.apply(len)
    if trim_length == 0:
        trim_length = seq_lengths.min()
        print(f"Auto-detected trim length: {trim_length} bp (minimum sequence length)")

    abs_trim = abs(trim_length)

    # Warn about and remove short sequences
    too_short = seq_lengths < abs_trim
    if too_short.any():
        n_short = too_short.sum()
        warnings.warn(
            f"Discarding {n_short} ASV(s) shorter than {abs_trim} bp.",
            UserWarning
        )
    rep_seqs = rep_seqs[~too_short].copy()

    # Filter count table to matching ASVs
    shared_ids = rep_seqs.index.intersection(counts.index)
    if len(shared_ids) == 0:
        raise ValueError(
            "No ASV IDs match between FASTA and count table. "
            f"FASTA IDs sample: {list(rep_seqs.index[:5])}, "
            f"Count table IDs sample: {list(counts.index[:5])}"
        )
    rep_seqs = rep_seqs.loc[shared_ids]
    counts = counts.loc[shared_ids]

    # Trim sequences
    if trim_length > 0:
        trimmed = rep_seqs.apply(lambda x: x[:trim_length])
    else:
        trimmed = rep_seqs.apply(lambda x: x[trim_length:])

    # Collapse duplicates: group ASVs that became identical after trimming
    # and sum their counts
    before = len(trimmed)
    seq_to_ids = trimmed.reset_index()
    seq_to_ids.columns = ['asv_id', 'sequence']

    # Group by the trimmed sequence
    groups = seq_to_ids.groupby('sequence')['asv_id'].apply(list)

    collapsed_seqs = {}
    collapsed_counts = pd.DataFrame(0, index=[], columns=counts.columns)

    for seq, asv_ids in groups.items():
        # Use first ASV ID as the representative
        rep_id = asv_ids[0]
        collapsed_seqs[rep_id] = seq
        collapsed_counts.loc[rep_id] = counts.loc[asv_ids].sum(axis=0)

    after = len(collapsed_seqs)
    if before != after:
        print(f"Trimming collapsed {before - after} duplicate sequence(s) "
              f"({before} -> {after} ASVs)")

    # Convert back to Series
    trimmed_seqs = pd.Series(collapsed_seqs)
    trimmed_seqs.index.name = 'feature-id'
    trimmed_seqs.name = 'sequence'

    # Save outputs
    write_fasta(trimmed_seqs, output_fasta)
    collapsed_counts.to_csv(output_counts, sep='\t')

    print(f"Trimmed {len(trimmed_seqs)} sequences to {abs_trim} bp")

    return trimmed_seqs, collapsed_counts


def prepare_extracted_region(sequences_file, region, trim_length, fwd_primer, rev_primer,
                             kmer_output, map_output, reverse_complement_rev=True,
                             reverse_complement_result=False, chunk_size=1000):
    """
    Prepare reference database for a specific region
    
    Parameters
    ----------
    sequences_file : str
        Path to reference FASTA (full 16S sequences)
    region : str
        Region name identifier
    trim_length : int
        Length to trim sequences to
    fwd_primer : str
        Forward primer sequence
    rev_primer : str
        Reverse primer sequence
    kmer_output : str
        Output path for k-mer FASTA
    map_output : str
        Output path for k-mer map TSV
    reverse_complement_rev : bool
        Whether to RC the reverse primer
    reverse_complement_result : bool
        Whether to RC the final extracted sequences
        
    Returns
    -------
    tuple
        (kmer sequences as Series, kmer map as DataFrame)
    """
    # RC reverse primer if needed
    if reverse_complement_rev:
        rev_primer = str(DNA(rev_primer).reverse_complement())

    # Read sequences
    sequences = read_fasta(sequences_file)

    # Expand degenerate bases, then trim to uniform length
    all_seqs = _block_seqs(sequences.to_dict())
    trimmed = _artifical_trim(all_seqs, trim_length)
    condensed = _condense_seqs(trimmed)
    
    # Collapse and optionally RC
    kmers, group2 = _collapse_all_sequences(condensed, reverse_complement_result)
    
    # Build k-mer map
    kmer_map = _expand_ids(group2, fwd_primer, rev_primer, region, trim_length, chunk_size)
    
    # Save outputs
    write_fasta(kmers, kmer_output)
    kmer_map.to_csv(map_output, sep='\t')
    
    return kmers, kmer_map


def _block_seqs(seqs, degen_thresh=3):
    """Expand degenerate sequences"""
    expanded = []
    for seq_id, seq_str in seqs.items():
        seq = DNA(seq_str)
        if seq.has_degenerates():
            exp = [str(s) for s in seq.expand_degenerates()]
            exp_series = pd.Series(exp)
            exp_series.index = [f'{seq_id}@{str(i+1).zfill(4)}' for i in range(len(exp))]
        else:
            exp_series = pd.Series([str(seq)], index=[seq_id])
        expanded.append(exp_series)
    
    result = pd.concat(expanded)
    result.name = 'sequence'
    result.index.name = 'seq-name'
    df = result.reset_index()
    df['db-seq'] = df['seq-name'].apply(lambda x: x.split("@")[0])
    return df


def _artifical_trim(seqs, trim_length):
    """Trim sequences to specified length"""
    seqs['extract-length'] = seqs['sequence'].apply(len)
    seqs = seqs[seqs['extract-length'].abs() >= abs(trim_length)]
    
    if trim_length > 0:
        seqs['amplicon'] = seqs['sequence'].apply(lambda x: x[:trim_length])
    else:
        seqs['amplicon'] = seqs['sequence'].apply(lambda x: x[trim_length:])
    
    seqs['db-seq'] = seqs['seq-name'].apply(lambda x: x.split('@')[0])
    seqs.drop_duplicates(['db-seq', 'amplicon'], inplace=True)
    
    return seqs[['seq-name', 'amplicon']]


def _condense_seqs(seqs):
    """Collapse duplicate sequences"""
    seqs.sort_values(['amplicon', 'seq-name'], inplace=True)
    fragment = seqs.groupby('amplicon')['seq-name'].apply(lambda x: "|".join(x))
    return fragment.sort_index().reset_index()


def _collapse_all_sequences(condensed, reverse_complement_result):
    """Collapse and optionally RC sequences"""
    grouped = condensed.groupby('amplicon')['seq-name'].apply(
        lambda x: '>%s' % "|".join(np.sort(x.values))
    )
    group2 = grouped.reset_index()
    group2.sort_values('seq-name', inplace=True)
    
    if reverse_complement_result:
        group2['seq'] = group2['amplicon'].apply(lambda x: str(DNA(x).reverse_complement()))
    else:
        group2['seq'] = group2['amplicon']
    
    # Create k-mer sequences as Series
    kmers = group2.set_index('seq-name')['seq']
    kmers.index = kmers.index.str.replace('>', '')
    
    return kmers, group2


def _expand_ids(group2, fwd_primer, rev_primer, region, trim_length, chunk_size):
    """Build k-mer map from collapsed sequences"""
    ids = group2['seq-name'].apply(lambda x: x.strip('>'))
    
    # Split pipe-delimited IDs
    expanded = []
    for kmer_id in ids:
        for seq_name in kmer_id.split("|"):
            expanded.append({
                'kmer': kmer_id,
                'seq-name': seq_name,
                'db-seq': seq_name.split("@")[0],
                'fwd-primer': fwd_primer,
                'rev-primer': rev_primer,
                'region': region,
                'kmer-length': trim_length
            })
    
    kmer_map = pd.DataFrame(expanded)
    kmer_map.set_index('db-seq', inplace=True)
    kmer_map.sort_index(inplace=True)
    
    return kmer_map
