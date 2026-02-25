# Standalone SIDLE

Multi-region 16S rRNA reconstruction without QIIME2 dependencies.

Modified from [q2-sidle](https://github.com/jwdebelius/q2-sidle)  
Copyright (c) 2020, Justine Debelius. BSD-3-Clause License

## Installation

```bash
pip install numpy pandas biopython scikit-bio
```

## Usage

### Step 1: Prepare Regional Databases

Extract primer-specific regions from reference database:

```python
import sidle_standalone as sidle

sidle.prepare_extracted_region(
    sequences_file='silva_138_99.fasta',
    region='V3-V4',
    trim_length=400,
    fwd_primer='CCTACGGGNGGCWGCAG',  # 341F
    rev_primer='GACTACHVGGGTATCTAATCC',  # 805R
    kmer_output='V3V4_kmers.fasta',
    map_output='V3V4_map.tsv'
)
```

### Step 2: Align ASVs to K-mers

Align DADA2 ASVs to regional k-mer database:

```python
sidle.align_regional_kmers(
    kmers_file='V3V4_kmers.fasta',
    rep_seq_file='V3V4_asvs.fasta',  # from DADA2
    region='V3-V4',
    output_file='V3V4_alignment.tsv',
    max_mismatch=2
)
```

### Step 3: Reconstruct Counts

Combine multiple regions using EM algorithm:

```python
sidle.reconstruct_counts(
    region=['V3-V4', 'V4-V5'],
    alignment_files=['V3V4_alignment.tsv', 'V4V5_alignment.tsv'],
    count_files=['V3V4_counts.tsv', 'V4V5_counts.tsv'],
    database_map_file='database_map.tsv',
    database_summary_file='database_summary.tsv',
    output_file='reconstructed_counts.tsv',
    region_normalize='average'
)
```

### Step 4: Reconstruct Taxonomy

```python
sidle.reconstruct_taxonomy(
    reconstruction_map_file='database_map.tsv',
    taxonomy_file='silva_taxonomy.tsv',
    output_file='reconstructed_taxonomy.tsv',
    database='silva'
)
```

## Snakemake Integration

```python
rule sidle_prepare_db:
    input: "reference.fasta"
    output: 
        kmers="kmers/{region}_kmers.fasta",
        map="kmers/{region}_map.tsv"
    params:
        trim_len=lambda w: config['regions'][w.region]['trim_length'],
        fwd=lambda w: config['regions'][w.region]['fwd_primer'],
        rev=lambda w: config['regions'][w.region]['rev_primer']
    run:
        import sidle_standalone as sidle
        sidle.prepare_extracted_region(
            input[0], wildcards.region, params.trim_len,
            params.fwd, params.rev, output.kmers, output.map
        )

rule sidle_align:
    input:
        kmers="kmers/{region}_kmers.fasta",
        asvs="dada2/{region}_asvs.fasta"
    output: "alignments/{region}_alignment.tsv"
    run:
        import sidle_standalone as sidle
        sidle.align_regional_kmers(
            input.kmers, input.asvs, wildcards.region, output[0]
        )

rule sidle_reconstruct:
    input:
        alignments=expand("alignments/{region}_alignment.tsv", region=REGIONS),
        counts=expand("dada2/{region}_counts.tsv", region=REGIONS)
    output: "sidle_reconstructed_counts.tsv"
    params:
        regions=REGIONS
    run:
        import sidle_standalone as sidle
        sidle.reconstruct_counts(
            params.regions, input.alignments, input.counts,
            "database_map.tsv", "database_summary.tsv", output[0]
        )
```

## Requirements

- Python 3.7+
- numpy, pandas, biopython, scikit-bio

## Key Parameters

**max_mismatch**: Maximum nucleotide mismatches in alignment (default: 2)  
**per_nucleotide_error**: Sequencing error rate (default: 0.005)  
**min_abund**: Minimum relative abundance threshold (default: 1e-5)  
**region_normalize**: 'average' (default), 'weighted', or 'unweighted'

## Citation

If using this code, cite:
- Debelius et al. (2021) bioRxiv 2021.03.23.436606
- Fuks et al. (2018) Microbiome 6:17
