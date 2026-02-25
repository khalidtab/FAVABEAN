"""
Standalone SIDLE reconstruction with EM algorithm
From q2-sidle: https://github.com/jwdebelius/q2-sidle
Copyright (c) 2020, Justine Debelius. BSD-3-Clause License
"""
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

import numpy as np
import pandas as pd

from .utils import _check_regions


def reconstruct_counts(region, alignment_files, count_files, database_map_file, database_summary_file,
                       output_file, per_nucleotide_error=0.005, min_abund=1e-5, min_counts=1000,
                       region_normalize='average'):
    """
    Reconstruct counts from multiple regions using EM algorithm
    
    Parameters
    ----------
    region : list of str
        Region names in order
    alignment_files : list of str
        Paths to alignment TSV files (from align_regional_kmers)
    count_files : list of str
        Paths to count TSV files (ASV tables from DADA2)
    database_map_file : str
        Path to database map TSV
    database_summary_file : str
        Path to database summary TSV
    output_file : str
        Path for output reconstructed counts
    per_nucleotide_error : float
        Per-nucleotide sequencing error rate
    min_abund : float
        Minimum relative abundance threshold
    min_counts : int
        Minimum total counts per sample
    region_normalize : str
        'average', 'weighted', or 'unweighted'
        
    Returns
    -------
    pd.DataFrame
        Reconstructed counts (features × samples)
    """
    region_order, region_names, num_regions = _check_regions(region)
    
    # Load data
    database_map = pd.read_csv(database_map_file, sep='\t', index_col=0)
    database_summary = pd.read_csv(database_summary_file, sep='\t', index_col=0)
    database_summary.index.name = 'clean_name'
    
    if region_normalize == 'unweighted':
        database_summary['num-regions'] = 1
    
    # Load alignments
    regional_alignment = [pd.read_csv(f, sep='\t') for f in alignment_files]
    
    # Construct alignment matrix
    align_mat = _construct_align_mat(
        matches=regional_alignment,
        sequence_map=database_map['clean_name'] if 'clean_name' in database_map.columns else database_map.iloc[:, 0],
        seq_summary=database_summary,
        nucleotide_error=per_nucleotide_error
    )
    
    # Load count tables per region and concatenate
    # ASV IDs are region-specific (from seqtab_to_fasta.R), so we concatenate
    # rather than sum — preserving per-region identity
    count_tables = []
    for count_file, reg in zip(count_files, region):
        table = pd.read_csv(count_file, sep='\t', index_col=0)
        table = table.apply(pd.to_numeric, errors='coerce').fillna(0)
        count_tables.append(table)

    counts = pd.concat(count_tables, axis=0)
    # If any ASV ID appears in multiple regions (shouldn't normally), sum them
    if counts.index.duplicated().any():
        counts = counts.groupby(counts.index).sum()

    # Filter to aligned ASVs only
    aligned_asvs = set(align_mat['asv'].values)
    keep_asvs = sorted(aligned_asvs & set(counts.index))
    if len(keep_asvs) == 0:
        raise ValueError(
            "No ASVs from count tables match alignment results. "
            "Check that ASV IDs are consistent between FASTA and count files.\n"
            f"  Count table IDs (first 5): {list(counts.index[:5])}\n"
            f"  Alignment ASV IDs (first 5): {sorted(aligned_asvs)[:5]}"
        )
    counts = counts.loc[keep_asvs]
    print(f"Matched {len(keep_asvs)}/{len(aligned_asvs)} aligned ASVs to count table")

    # Filter low-count samples
    sample_totals = counts.sum(axis=0)
    keep_samples = sample_totals > min_counts
    if not keep_samples.all():
        n_removed = (~keep_samples).sum()
        print(f"Warning: Removing {n_removed} sample(s) with <{min_counts} total counts")
    counts = counts.loc[:, keep_samples]
    counts = counts.loc[counts.sum(axis=1) > 0]

    align_mat = align_mat[align_mat['asv'].isin(counts.index)]
    
    # Normalize counts
    n_table = counts / counts.sum(axis=0)
    
    # Solve EM for relative abundance
    rel_abund = _solve_iterative_noisy(
        align_mat=align_mat,
        table=n_table,
        min_abund=min_abund,
        seq_summary=database_summary
    )
    
    # Scale to counts
    count_table = _scale_relative_abundance(
        align_mat=align_mat,
        relative=rel_abund,
        counts=counts,
        region_normalize=region_normalize,
        num_regions=num_regions,
        seq_summary=database_summary
    )
    
    # Save and return
    count_table.to_csv(output_file, sep='\t')
    return count_table


def _construct_align_mat(matches, sequence_map, seq_summary, nucleotide_error=0.005):
    """Build alignment probability matrix.

    Uses merge() instead of .loc indexing to robustly handle cases where
    alignment db-seq values don't perfectly match the database map keys.
    """
    # Build a clean lookup DataFrame from sequence_map
    if isinstance(sequence_map, pd.Series):
        map_df = sequence_map.reset_index()
        map_df.columns = ['db-seq', 'clean_name']
    else:
        map_df = sequence_map.reset_index()
        if 'clean_name' not in map_df.columns:
            map_df.columns = ['db-seq', 'clean_name']

    align_list = []

    for match in matches:
        # Expand pipe-delimited kmer IDs into separate rows
        match = _expand_duplicate_sequences(match, 'kmer')

        # Extract the base reference ID (before @)
        match['db-seq'] = match['kmer'].apply(lambda x: str(x).split("@")[0])

        # Map to clean names via merge (robust to missing keys)
        match = match.merge(map_df, on='db-seq', how='inner')

        if len(match) == 0:
            warnings.warn("No matches found for an alignment batch — "
                          "check that kmer IDs match database map entries.")
            continue

        align_list.append(match)

    if len(align_list) == 0:
        raise ValueError(
            "No alignments could be mapped to database sequences. "
            "Check that the database map was built from the same "
            "reference and alignment files."
        )

    align_mat = pd.concat(align_list, axis=0)
    align_mat.drop_duplicates(['clean_name', 'asv', 'region'], inplace=True)

    # Get best match per ASV-seq-region combo
    asv_ref = align_mat.groupby(
        ['clean_name', 'asv', 'region', 'max-mismatch']
    )[['mismatch', 'length']].min().reset_index()

    # Calculate match probabilities (SMURF algorithm)
    asv_ref['match_prob'] = (
        np.power(1 - nucleotide_error, asv_ref['length'] - asv_ref['mismatch']) *
        np.power(nucleotide_error / 3, asv_ref['mismatch'])
    ).astype(float)

    asv_ref['perfect_prob'] = np.power(1 - nucleotide_error, asv_ref['length'])
    asv_ref['max_error_prob'] = (
        np.power(1 - nucleotide_error, asv_ref['length'] - asv_ref['max-mismatch']) *
        np.power(nucleotide_error / 3, asv_ref['max-mismatch'])
    )

    asv_ref['error_thresh'] = (
        asv_ref['perfect_prob'] - 0.1 * asv_ref['max_error_prob'] - np.spacing(1)
    )
    asv_ref['perf_match'] = asv_ref['match_prob'] > asv_ref['error_thresh']

    # Normalize by region coverage — use merge instead of .loc
    region_norm = seq_summary[['num-regions']].copy()
    region_norm.columns = ['region_norm']
    asv_ref = asv_ref.merge(
        region_norm, left_on='clean_name', right_index=True, how='left'
    )
    asv_ref['region_norm'] = asv_ref['region_norm'].fillna(1).astype(int)
    asv_ref['norm'] = asv_ref['match_prob'] / asv_ref['region_norm']
    asv_ref.drop(columns='region_norm', inplace=True)

    return asv_ref


def _solve_iterative_noisy(align_mat, table, seq_summary, tolerance=1e-7, min_abund=1e-10, num_iter=10000):
    """Solve EM for all samples"""
    align_mat = align_mat[align_mat['asv'].isin(table.index)]
    align = align_mat.pivot_table(values='norm', index='asv', columns='clean_name', fill_value=0)
    # Only keep ASVs present in both the alignment and the count table
    shared_idx = align.index.intersection(table.index)
    align = align.loc[shared_idx]
    table = table.loc[shared_idx]
    
    align_seqs = align.columns.values
    recon_list = []
    
    for sample in table.columns:
        col = table[sample]
        filt_align = align.loc[col > 0].values
        abund = col[col > 0].values
        sub_seqs = align_seqs[(filt_align > 0).any(axis=0)]
        filt_align = filt_align[:, (filt_align > 0).sum(axis=0) > 0]
        
        freq = _solve_ml_em_iterative_1_sample(
            align=filt_align,
            abund=abund,
            align_kmers=sub_seqs,
            sample=sample,
            num_iter=num_iter,
            tolerance=tolerance,
            min_abund=min_abund
        )
        recon_list.append(freq)
    
    # Combine samples
    recon = pd.concat(recon_list, axis=1).fillna(0)
    
    # Normalize by regions
    recon_norm = recon.copy()
    for seq_id in recon.index:
        if seq_id in seq_summary.index:
            recon_norm.loc[seq_id] /= seq_summary.loc[seq_id, 'num-regions']
    
    # Normalize to sum=1 per sample
    recon_norm = recon_norm / recon_norm.sum(axis=0)
    
    return recon_norm


def _solve_ml_em_iterative_1_sample(align, abund, align_kmers, sample, num_iter=10000, tolerance=1e-7, min_abund=1e-10):
    """
    Core EM algorithm from SMURF
    
    Expectation-Maximization to reconstruct bacterial frequencies from ASV counts
    """
    # Initial frequency estimate
    bact_freq = np.dot(align.T, abund) / np.sum(np.dot(align.T, abund))
    
    # EM iterations
    for i in range(int(num_iter)):
        # E-step: assign theta estimate
        theta_i = np.dot(align, bact_freq)
        
        # M-step: update frequencies
        r_weighted = abund / (theta_i + np.spacing(1))
        bact_factor = np.dot(r_weighted.T, align)
        
        # Check convergence
        error = np.dot(np.absolute(1 - bact_factor), bact_freq)
        bact_freq = bact_freq * bact_factor
        
        if error < tolerance:
            break
        
        # Filter low-abundance
        high_enough = bact_freq > min_abund
        align = align[:, high_enough]
        bact_freq = bact_freq[high_enough]
        align_kmers = align_kmers[high_enough]
    
    # Hard threshold and normalize
    bact_freq[bact_freq <= min_abund] = 0
    bact_freq = bact_freq / bact_freq.sum()
    
    return pd.Series(bact_freq, index=align_kmers, name=sample)


def _scale_relative_abundance(align_mat, relative, counts, seq_summary, num_regions, region_normalize='average'):
    """Scale relative abundance to absolute counts"""
    align_mat = align_mat[align_mat['clean_name'].isin(relative.index) & align_mat['asv'].isin(counts.index)]
    align = align_mat.pivot_table(values='norm', index='asv', columns='clean_name', fill_value=0)
    align = align[relative.index]

    # Intersect indices so boolean masks match DataFrame dimensions
    # (mirrors the shared_idx logic in _solve_iterative_noisy)
    shared_asvs = align.index.intersection(counts.index)
    align = align.loc[shared_asvs]
    counts = counts.loc[shared_asvs]

    scaled_list = []
    for sample in counts.columns:
        non_zero_asvs = (counts[sample] > 0).values
        non_zero_seqs = (relative[sample] > 0).values
        
        if non_zero_asvs.sum() == 0 or non_zero_seqs.sum() == 0:
            continue
        
        sub_align = align.iloc[non_zero_asvs, non_zero_seqs].values
        sub_freq = relative.iloc[non_zero_seqs][sample].values
        sub_counts = counts.loc[non_zero_asvs, sample].values
        sub_seqs = relative.index[non_zero_seqs]
        
        # Probability calculations
        p_r_and_j = sub_freq * sub_align
        p_r_given_j = p_r_and_j / (p_r_and_j.sum(axis=1, keepdims=True) + np.spacing(1))
        count_j_given_r = sub_counts[:, np.newaxis] * p_r_given_j
        
        scaled_list.append(pd.Series(count_j_given_r.sum(axis=0), index=sub_seqs, name=sample))
    
    scaled_counts = pd.DataFrame(scaled_list).T.fillna(0)
    
    # Apply region normalization
    for seq_id in scaled_counts.index:
        if seq_id in seq_summary.index:
            n_reg = seq_summary.loc[seq_id, 'num-regions']
            if region_normalize == 'average':
                scaled_counts.loc[seq_id] = np.round(scaled_counts.loc[seq_id] / n_reg, 0)
            elif region_normalize == 'weighted':
                scaled_counts.loc[seq_id] = np.round(scaled_counts.loc[seq_id] * num_regions / n_reg, 0)
    
    return scaled_counts


def _expand_duplicate_sequences(df, id_col, delim='|'):
    """Expand pipe-delimited IDs into separate rows"""
    data_cols = [c for c in df.columns if c != id_col]
    wide_exp = pd.concat([
        df[id_col].apply(lambda x: pd.Series(str(x).split(delim))),
        df[data_cols]
    ], axis=1)
    
    long_exp = wide_exp.melt(id_vars=data_cols, value_name=id_col)
    long_exp = long_exp.drop(columns=['variable']).dropna(subset=[id_col])
    long_exp = long_exp[[id_col] + data_cols]
    long_exp.sort_values(by=list(long_exp.columns), inplace=True)
    long_exp.reset_index(drop=True, inplace=True)
    
    return long_exp
