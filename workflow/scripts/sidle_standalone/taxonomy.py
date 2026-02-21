"""
Standalone SIDLE taxonomy reconstruction
From q2-sidle: https://github.com/jwdebelius/q2-sidle
Copyright (c) 2020, Justine Debelius. BSD-3-Clause License
"""
import warnings
import numpy as np
import pandas as pd

from .utils import database_params


def reconstruct_taxonomy(reconstruction_map_file, taxonomy_file, output_file,
                         database='none', define_missing='merge', ambiguity_handling='ignore'):
    """
    Reconstruct taxonomy for reconstructed sequences
    
    Parameters
    ----------
    reconstruction_map_file : str
        Path to reconstruction map TSV (db-seq → clean_name)
    taxonomy_file : str
        Path to taxonomy TSV (seq_id → taxonomy string)
    output_file : str
        Output path for reconstructed taxonomy
    database : str
        Database type: 'greengenes', 'silva', or 'none'
    define_missing : str
        How to handle missing taxa: 'merge', 'inherit', or 'ignore'
    ambiguity_handling : str
        How to handle ambiguous taxa: 'missing' or 'ignore'
        
    Returns
    -------
    pd.Series
        Reconstructed taxonomy strings
    """
    # Validate parameters
    if database == 'none' and define_missing != 'ignore':
        warnings.warn('When no database specified, missing values ignored by default')
    if database == 'greengenes' and ambiguity_handling != 'ignore':
        warnings.warn('Greengenes does not include ambiguous taxa')
    
    # Load data
    recon_map = pd.read_csv(reconstruction_map_file, sep='\t', index_col=0)
    taxonomy = pd.read_csv(taxonomy_file, sep='\t', index_col=0, header=None, names=['taxonomy'])
    
    # Get clean_name column
    clean_name_col = 'clean_name' if 'clean_name' in recon_map.columns else recon_map.columns[0]
    
    # Get database parameters
    db_lookup = database_params.get(database, database_params['none'])
    delim = db_lookup['delim']
    defined_f = db_lookup['defined']
    
    # Split taxonomy into levels
    taxonomy = taxonomy.loc[recon_map.index]
    tax_split = taxonomy['taxonomy'].apply(lambda x: pd.Series([s.strip() for s in str(x).split(delim)]))
    
    if len(tax_split.columns) == 1:
        raise ValueError('Only one taxonomic level found. Check database and delimiter.')
    
    # Find undefined levels
    undefined_levels = pd.concat([
        tax_split[c].apply(lambda x: not defined_f(x)) for c in tax_split.columns
    ], axis=1)
    
    ambiguous_levels = pd.concat([
        tax_split[c].apply(lambda x: 'ambig' in str(x).lower()) for c in tax_split.columns
    ], axis=1)
    ambiguous_levels = ambiguous_levels.cummax(axis=1)
    
    undefined = (undefined_levels | (ambiguous_levels & (ambiguity_handling == 'missing'))).astype(bool)
    
    # Handle missing taxa
    if define_missing != 'ignore':
        tax_split.mask(undefined, np.nan, inplace=True)
    
    if define_missing == 'inherit':
        tax_split = tax_split.ffill(axis=1)
    
    # Combine taxonomy for same clean_name
    def _combine_f(x):
        if pd.isnull(x).all():
            return np.nan
        return '|'.join(np.sort(x.dropna().unique()))
    
    def _combine_taxa(g):
        if len(g) == 1:
            return g.iloc[0]
        return g.apply(_combine_f)
    
    tax_split['clean_name'] = recon_map[clean_name_col]
    collapsed = tax_split.groupby('clean_name').apply(_combine_taxa)
    # In newer pandas, 'clean_name' becomes the index after groupby; drop it
    # from columns only if it's still there
    if 'clean_name' in collapsed.columns:
        collapsed.drop(columns=['clean_name'], inplace=True)
    
    # Find splits
    disjoint = pd.concat([
        collapsed[c].apply(lambda x: True if pd.isnull(x) else '|' in str(x))
        for c in collapsed.columns
    ], axis=1).cummax(axis=1)
    
    disjoint_inherit = (disjoint.cummax(axis=1) & ~((disjoint.cumsum(axis=1) == 1) & disjoint))
    collapsed.mask(disjoint_inherit, np.nan, inplace=True)
    
    # Forward fill NaNs
    collapsed = collapsed.ffill(axis=1)
    
    # Combine into taxonomy strings
    new_taxa = collapsed.apply(lambda x: delim.join([str(v) for v in x.values]), axis=1)
    new_taxa.name = 'Taxon'
    new_taxa.index.name = 'Feature ID'
    
    # Save
    new_taxa.to_csv(output_file, sep='\t', header=True)
    
    return new_taxa
