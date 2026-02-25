"""
Standalone SIDLE database reconstruction
From q2-sidle: https://github.com/jwdebelius/q2-sidle
Copyright (c) 2020, Justine Debelius. BSD-3-Clause License

Reconstructs a kmer database based on aligned regions, resolving
ambiguous mappings where multiple reference sequences share identical
k-mers across regions.

The algorithm works in three rounds:
  Round 1: Group k-mers by reference sequence, intersect across regions,
           identify "tidy" groups (1-to-1 kmer-to-clean_name mapping)
  Round 2: Re-process untidy sequences excluding already-assigned ones
  Round 3: Graph-based clustering for remaining ambiguous sequences
"""
import warnings
import numpy as np
import pandas as pd

from .utils import _check_regions


def reconstruct_database(region, alignment_files, kmer_map_files,
                         output_map_file, output_summary_file,
                         count_degenerates=True):
    """
    Reconstruct the k-mer database map from aligned regions.

    This is the critical step that resolves which reference sequences
    are distinguishable from the observed data. Without it, the EM
    algorithm receives many-to-many mappings that degrade accuracy.

    Parameters
    ----------
    region : list of str
        Region identifiers
    alignment_files : list of str
        Paths to per-region alignment TSV files (output of align_regional_kmers)
    kmer_map_files : list of str
        Paths to per-region k-mer map TSV files (output of prepare_extracted_region)
    output_map_file : str
        Path for output database map TSV
    output_summary_file : str
        Path for output database summary TSV
    count_degenerates : bool
        Whether degenerate sequences count as unique k-mers in the
        summary statistics

    Returns
    -------
    reconstruction : pd.DataFrame
        Database map: index=db-seq, columns include clean_name,
        first/last region info
    summary : pd.DataFrame
        Database summary: index=feature-id (clean_name),
        columns: num-regions, total-kmers-mapped,
        mean-kmer-per-region, stdv-kmer-per-region
    """
    region_order, region_names, num_regions = _check_regions(region)

    # Load alignment and kmer map files
    regional_alignment = [pd.read_csv(f, sep='\t') for f in alignment_files]
    kmer_map = [pd.read_csv(f, sep='\t', index_col=0) for f in kmer_map_files]

    # ----------------------------------------------------------------
    # Step 1: Get the set of k-mers that actually aligned to ASVs
    # ----------------------------------------------------------------
    align_map = pd.concat(regional_alignment, axis=0, sort=False)
    align_map.drop_duplicates(['asv', 'kmer'], inplace=True)
    align_map['region'] = align_map['region'].replace(region_order)
    aligned_kmers = _get_unique_kmers(align_map['kmer'])
    print('Alignment map loaded')

    # ----------------------------------------------------------------
    # Step 2: Filter k-mer maps to only aligned k-mers
    # ----------------------------------------------------------------
    kmer_alignments = [_filter_to_aligned(km, aligned_kmers) for km in kmer_map]

    # Get the full set of database sequences that were aligned
    db_seqs_list = [_get_unique_kmers(df['kmer']) for df in kmer_alignments]
    db_seqs = np.unique(np.hstack(db_seqs_list))
    print(f'Database k-mers identified: {len(db_seqs)} sequences')

    # Build the region database: for each db-seq, collect all k-mers per region
    mapped_kmer = pd.concat(
        [_build_region_db(df, db_seqs) for df in kmer_alignments],
        axis=0
    )
    mapped_kmer['region'] = mapped_kmer['region'].replace(region_order)

    # Clean k-mer lists: collapse degenerates (remove @0001 suffixes)
    # and filter to only sequences in db_seqs
    db_seq_set = set(db_seqs)
    mapped_kmer['kmer'] = mapped_kmer['kmer'].apply(_clean_kmer_list)
    mapped_kmer['kmer'] = mapped_kmer['kmer'].apply(
        _check_db_list, ref_seqs=db_seq_set
    )
    db_seqs = mapped_kmer['db-seq'].unique()
    print('K-mer map built')

    # ----------------------------------------------------------------
    # Step 3: Combine sequences into regional groups
    # ----------------------------------------------------------------
    # For each db-seq + region combo, get the union of all k-mers.
    # Use '@@' as separator instead of '-' because db-seq (taxonomy strings)
    # and region names (e.g. "V3-V4") can both contain hyphens.
    _COMPOSITE_SEP = '@@'
    mapped_kmer['composite'] = mapped_kmer.apply(
        lambda x: f"{x['db-seq']}{_COMPOSITE_SEP}{x['region']}", axis=1
    )
    mapped_kmer['kmer'] = mapped_kmer['kmer'].apply(_clean_kmer_list)
    mapped_kmer['kmer'] = mapped_kmer['kmer'].apply(
        _check_db_list, ref_seqs=db_seq_set
    )

    region_db_series = mapped_kmer.groupby('composite')['kmer'].apply(
        _get_regional_seqs
    )
    region_db = region_db_series.reset_index()
    region_db = region_db[['composite', region_db.columns[-1]]]
    region_db.columns = ['composite', 'kmer']
    region_db['region'] = region_db['composite'].apply(
        lambda x: x.split(_COMPOSITE_SEP)[-1]
    )
    region_db['db-seq'] = region_db['composite'].apply(
        lambda x: _COMPOSITE_SEP.join(x.split(_COMPOSITE_SEP)[:-1])
    )
    print('Regions grouped')

    # ----------------------------------------------------------------
    # Step 4: Three rounds of "detangling"
    # ----------------------------------------------------------------
    # Round 1: Find tidy mappings (k-mer -> single clean_name)
    print('Tidying round 1')
    shared_kmers_1, shared_check_1 = _check_intersection(region_db, {'X'})

    tidy_kmer_1 = shared_check_1.loc[shared_check_1['tidy'], 'kmer'].unique()
    tidy_seqs_1 = shared_kmers_1.loc[
        shared_kmers_1['kmer'].isin(tidy_kmer_1), 'db-seq'
    ].unique()
    tidy_maps = [shared_kmers_1.loc[shared_kmers_1['db-seq'].isin(tidy_seqs_1)]]
    print(f'Round 1: {len(tidy_seqs_1)} tidy sequences')

    # Round 2: Re-process excluding round 1 results
    print('Tidying round 2')
    region_db_2 = region_db.loc[~region_db['db-seq'].isin(tidy_seqs_1)].copy()
    shared_kmers_2, shared_check_2 = _check_intersection(
        region_db_2, {'X'} | set(tidy_seqs_1)
    )

    tidy_kmer_2 = shared_check_2.loc[shared_check_2['tidy'], 'kmer'].unique()
    tidy_seqs_2 = shared_kmers_2.loc[
        shared_kmers_2['kmer'].isin(tidy_kmer_2), 'db-seq'
    ].unique()
    tidy_maps.append(
        shared_kmers_2.loc[shared_kmers_2['db-seq'].isin(tidy_seqs_2)]
    )
    print(f'Round 2: {len(tidy_seqs_2)} additional tidy sequences')

    # Round 3: Graph-based clustering for remaining ambiguous sequences
    print('Tidying round 3')
    unkempt = region_db_2.loc[~region_db_2['db-seq'].isin(tidy_seqs_2)]
    matched_so_far = {'X'} | set(tidy_seqs_1) | set(tidy_seqs_2)

    unkempt_grouped = unkempt.groupby('db-seq').apply(
        _define_shared, matched=matched_so_far
    )
    # Ensure unkempt_grouped is a clean Series (drop extra index levels from newer pandas)
    if isinstance(unkempt_grouped, pd.DataFrame):
        unkempt_grouped = unkempt_grouped.iloc[:, -1]

    if len(unkempt_grouped) > 0 and unkempt_grouped.str.len().sum() > 0:
        # Expand pipe-delimited mappings into long form
        to_map = unkempt_grouped.apply(
            lambda x: pd.Series(x.split("|")) if x else pd.Series(dtype=str)
        ).unstack().dropna()
        to_map = to_map.reset_index()
        # Take only the last 3 columns (handles extra index cols from newer pandas)
        to_map = to_map.iloc[:, -3:]
        to_map.columns = ['counter', 'db-seq', 'clean_name']
        detangled = _detangle_names(to_map.copy())
        detangled = detangled.reset_index()
        detangled.rename(columns={'clean_name': 'kmer'}, inplace=True)
    else:
        detangled = pd.DataFrame(columns=['db-seq', 'kmer'])

    print(f'Round 3: {len(detangled)} detangled sequences')

    # ----------------------------------------------------------------
    # Step 5: Combine all rounds into final database map
    # ----------------------------------------------------------------
    print('Combining solutions')
    combined = pd.concat(tidy_maps + [detangled], axis=0)
    combined.rename(columns={'kmer': 'clean_name'}, inplace=True)
    combined = combined[['db-seq', 'clean_name']].drop_duplicates()
    combined.set_index('db-seq', inplace=True)

    # Build the reconstruction map with first/last region info
    all_kmer_data = pd.concat(kmer_map, axis=0)
    if 'db-seq' not in all_kmer_data.columns and all_kmer_data.index.name == 'db-seq':
        all_kmer_data = all_kmer_data.reset_index()

    tidy_db = all_kmer_data.loc[all_kmer_data['db-seq'].isin(combined.index)].copy()
    tidy_db.set_index('db-seq', inplace=True)
    tidy_db['clean_name'] = combined['clean_name']
    tidy_db['region'] = tidy_db['region'].replace(region_order)
    tidy_db.sort_values(['db-seq', 'region'], ascending=True, inplace=True)

    # Build reconstruction: first and last region per db-seq
    first_ = tidy_db.groupby(level=0, sort=False).first()
    last_ = tidy_db.groupby(level=0, sort=False).last()

    reconstruction = pd.concat([
        first_['clean_name'],
        first_[['region', 'fwd-primer']].add_prefix('first-'),
        last_[['region', 'fwd-primer', 'kmer-length']].add_prefix('last-')
    ], axis=1)

    print('Database map complete')

    # ----------------------------------------------------------------
    # Step 6: Build summary statistics
    # ----------------------------------------------------------------
    db_summary = _count_mapping(
        tidy_db.reset_index(),
        count_degenerates,
        kmer='seq-name' if 'seq-name' in tidy_db.reset_index().columns else 'kmer'
    )

    # Map aligned ASVs to clean names
    aligned_asvs = _map_aligned_asvs(align_map, combined)
    if len(aligned_asvs) > 0:
        db_summary['mapped-asvs'] = aligned_asvs.groupby('clean_name').apply(
            lambda x: '|'.join(x)
        )

    db_summary.index.set_names('feature-id', inplace=True)
    print(f'Database summarized: {len(db_summary)} features')

    # Save outputs
    reconstruction.to_csv(output_map_file, sep='\t')
    db_summary.to_csv(output_summary_file, sep='\t')

    return reconstruction, db_summary


# ==================================================================
# Helper functions
# ==================================================================

def _get_unique_kmers(series):
    """Extract unique base reference IDs from pipe-delimited k-mer strings"""
    kmers = np.hstack([
        [a.split("@")[0] for a in kmer.split('|')]
        for kmer in series
    ])
    return np.sort(np.unique(kmers))


def _filter_to_aligned(df, aligned_kmers):
    """Filter k-mer map to only k-mers that had alignments"""
    if df.index.name == 'db-seq':
        df = df.reset_index()
    keep = df.index[df.index.isin(aligned_kmers)] if df.index.name else None

    # Try filtering by index (if index is k-mer IDs)
    if 'kmer' in df.columns:
        # Check which rows have k-mers overlapping with aligned set
        mask = df['kmer'].apply(
            lambda x: any(k.split("@")[0] in aligned_kmers
                          for k in str(x).split("|"))
        )
        return df.loc[mask].copy()
    else:
        mask = df.index.isin(aligned_kmers)
        return df.loc[mask].reset_index()[['db-seq', 'region', 'kmer']]


def _build_region_db(df, kmers):
    """Filter to sequences in the database and return core columns"""
    if 'db-seq' not in df.columns and df.index.name == 'db-seq':
        df = df.reset_index()
    filtered = df.loc[df['db-seq'].isin(kmers)].copy()
    return filtered[['db-seq', 'region', 'kmer']]


def _clean_kmer_list(x):
    """Remove degenerate suffixes (@0001) from pipe-delimited k-mer IDs"""
    return '|'.join(np.unique([y.split('@')[0] for y in str(x).split("|")]))


def _check_db_list(x, ref_seqs):
    """Filter out sequence IDs not in the reference database"""
    f_ = lambda y: y if y in ref_seqs else 'X'
    return '|'.join(np.unique([f_(y) for y in str(x).split('|')]))


def _get_regional_seqs(x):
    """Get union of all k-mer sequences within a region group"""
    degenerates = np.unique(np.hstack([y.split("|") for y in x]))
    return '|'.join(degenerates)


def _define_shared(g, matched):
    """
    For a group of k-mers belonging to one db-seq, find the
    intersection of k-mer sets across all regions, excluding
    already-matched sequences.
    """
    regional_sets = g['kmer'].apply(lambda x: set(str(x).split("|"))).values
    if len(regional_sets) == 0:
        return ''
    intersect = set.intersection(*regional_sets) - matched
    return '|'.join(sorted(intersect))


def _check_intersection(df, matched=frozenset({'X'})):
    """
    Identify tidy vs untidy k-mer mappings.

    A mapping is "tidy" when the k-mer group name equals the
    pipe-joined sorted list of db-seqs that share it — meaning
    there's a 1-to-1 relationship.
    """
    matched = set(matched)

    shared_series = df.groupby('db-seq').apply(
        _define_shared, matched=matched
    )
    shared_kmers = shared_series.reset_index()
    # Keep only the group key and the result value (handles varying pandas versions)
    shared_kmers = shared_kmers[['db-seq', shared_kmers.columns[-1]]]
    shared_kmers.columns = ['db-seq', 'kmer']

    # Drop empty strings (no shared k-mers after exclusion)
    shared_kmers = shared_kmers[shared_kmers['kmer'].str.len() > 0]

    if len(shared_kmers) == 0:
        shared_check = pd.DataFrame(columns=['kmer', 'shared-name', 'tidy'])
        return shared_kmers, shared_check

    # Check tidiness: group by k-mer, join all db-seqs that share it
    shared_check_series = shared_kmers.groupby('kmer')['db-seq'].apply(
        lambda x: '|'.join(sorted(x))
    )
    shared_check = shared_check_series.reset_index()
    shared_check = shared_check[['kmer', shared_check.columns[-1]]]
    shared_check.columns = ['kmer', 'shared-name']

    # Tidy = the k-mer group IS the sorted db-seq list
    # (meaning each k-mer maps to exactly one group of sequences)
    shared_check['tidy'] = shared_check['kmer'] == shared_check['shared-name']

    return shared_kmers, shared_check


def _detangle_names(long_):
    """
    Graph-based clustering for ambiguous sequence mappings.

    When multiple reference sequences share overlapping k-mers across
    regions, this function clusters them by building a co-occurrence
    matrix and finding connected components with matching connection
    patterns.

    Parameters
    ----------
    long_ : DataFrame
        Columns: counter, db-seq, clean_name
        A long-form mapping between database sequences and their
        candidate clean names.

    Returns
    -------
    Series
        Mapping from db-seq to a pipe-delimited clean_name group.
    """
    long_.sort_values(['db-seq', 'clean_name'], ascending=True, inplace=True)
    long_['counter'] = 1
    long_['counter'] = long_.groupby('db-seq')['counter'].cumsum() - 1

    # Build a square co-occurrence matrix
    map_ = long_.set_index(['db-seq', 'clean_name'])['counter'].unstack().notna()
    diagonal = (np.tril(np.ones(map_.shape)) == np.triu(np.ones(map_.shape)))
    map_ = (map_ | diagonal)

    # Sort to cluster similar patterns together
    map_.sort_values(list(map_.columns), inplace=True)
    map_ = map_[map_.index]

    # Calculate connection counts for each sequence
    num_connections = long_.groupby('db-seq')['counter'].max() + 1
    long_.set_index('clean_name', inplace=True)
    long_['co-connection'] = num_connections[long_.index]
    co_min = long_.groupby('db-seq')['co-connection'].min()
    connection_penalty = num_connections - co_min + 1

    # Build penalty mask: only connect sequences with matching
    # connection counts
    penalty_v = (
        np.ones((len(connection_penalty), 1)) *
        connection_penalty.loc[map_.index].values
    )
    penalty_mask = (penalty_v == penalty_v.T)

    overlap = penalty_mask & map_
    overlap = overlap.mask(~overlap, np.nan)
    overlap = overlap.unstack().dropna().reset_index()

    if len(overlap) == 0:
        return pd.Series(dtype=str, name='clean_name')

    # Verify that paired sequences actually share the same connections
    shared = overlap.groupby('db-seq')['clean_name'].apply(
        lambda x: set(x.values)
    )

    def check_shared(x):
        seqs0 = shared.loc[x['clean_name']]
        seqs1 = shared.loc[x['db-seq']]
        if len(seqs0) != len(seqs1):
            return 0
        return len(seqs0 & seqs1) - 1

    check_diff = overlap['clean_name'] != overlap['db-seq']

    if check_diff.any():
        overlap.loc[check_diff, 'seq_count'] = overlap.loc[check_diff].apply(
            check_shared, axis=1
        )
    else:
        overlap['seq_count'] = np.nan

    overlap['matches'] = (overlap['seq_count'] == overlap[0]) | ~check_diff

    overlap.sort_values(['db-seq', 'clean_name'], inplace=True, ascending=True)
    overlap.drop_duplicates(['db-seq', 'clean_name'], inplace=True)

    new_name = overlap.groupby('db-seq')['clean_name'].apply(
        lambda x: "|".join(x)
    )

    return new_name


def _count_mapping(long_, count_degen, kmer='kmer', region='region',
                   clean_name='clean_name', db_seq='db-seq'):
    """
    Count the number of k-mers mapped per region per clean_name.

    Returns summary statistics needed by the EM algorithm for
    region normalization.
    """
    if kmer not in long_.columns:
        kmer = 'kmer' if 'kmer' in long_.columns else long_.columns[0]
    if clean_name not in long_.columns:
        return pd.DataFrame(
            columns=['num-regions', 'total-kmers-mapped',
                     'mean-kmer-per-region', 'stdv-kmer-per-region']
        )

    long_ = long_.copy()

    if count_degen:
        long_.drop_duplicates([kmer, region, clean_name], inplace=True)
    else:
        if db_seq in long_.columns:
            long_.drop_duplicates([db_seq, region, clean_name], inplace=True)

    regional_seqs = long_.groupby([clean_name, region]).count()[[kmer]]
    regional_seqs.reset_index(inplace=True)

    counts = pd.DataFrame({
        'num-regions': regional_seqs.groupby(clean_name)[region].count(),
        'total-kmers-mapped': regional_seqs.groupby(clean_name)[kmer].sum(),
        'mean-kmer-per-region': regional_seqs.groupby(clean_name)[kmer].mean(),
        'stdv-kmer-per-region': regional_seqs.groupby(clean_name)[kmer].std().fillna(0),
    })

    return counts


def _map_aligned_asvs(align_map, seq_map):
    """Map aligned ASVs back to clean names via the sequence map"""
    align_map = align_map.copy()
    align_map.sort_values(['asv', 'region'], inplace=True, ascending=True)

    # Clean k-mer IDs
    aligned_asvs = align_map.set_index(['asv', 'region'])['kmer'].apply(
        _clean_kmer_list
    )
    aligned_asvs = aligned_asvs.apply(
        lambda x: pd.Series(str(x).split("|"))
    )
    aligned_asvs.reset_index(inplace=True)
    aligned_asvs = aligned_asvs.melt(
        id_vars=['asv', 'region'],
        var_name='counter',
        value_name='db-seq'
    )
    aligned_asvs.dropna(inplace=True)
    aligned_asvs.set_index('db-seq', inplace=True)

    # Map to clean names
    keep = aligned_asvs.index.isin(seq_map.index)
    aligned_asvs = aligned_asvs.loc[keep]

    if len(aligned_asvs) == 0:
        return pd.Series(dtype=str, name='asv')

    aligned_asvs['clean_name'] = seq_map.loc[aligned_asvs.index, 'clean_name']
    aligned_asvs.sort_values(['clean_name', 'region', 'asv'], inplace=True)
    aligned_asvs.drop_duplicates(['clean_name', 'asv'], inplace=True)

    return aligned_asvs.set_index('clean_name')['asv']
