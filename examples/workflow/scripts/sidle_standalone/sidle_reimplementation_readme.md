# SIDLE Reimplementation: Changes from q2-sidle

This document describes the changes made when porting the SIDLE (SMURF) multi-region 16S reconstruction algorithm from the original QIIME2 plugin ([q2-sidle](https://github.com/jwdebelius/q2-sidle), Debelius et al. 2021) into a standalone Python package for integration with the FAVABEAN Snakemake pipeline.

Original code copyright (c) 2020, Justine Debelius. Licensed under BSD-3-Clause (see `LICENSE`).

---

## Motivation

The original q2-sidle is tightly coupled to the QIIME2 ecosystem: it uses `biom.Table` for count data, `qiime2.Artifact` for I/O, `qiime2.Metadata` for summaries, and `DNAFASTAFormat` / `DNAIterator` for sequence handling. FAVABEAN uses Snakemake + DADA2 and produces TSV count tables and FASTA files directly. Rather than require users to install the full QIIME2 framework as an intermediate dependency, we refactored the core SIDLE algorithms to operate on standard file formats (TSV, FASTA) with no QIIME2 dependency.

---

## Architectural Changes

### 1. Complete removal of dask

**Original:** q2-sidle uses `dask.delayed` and `dask.distributed.Client` throughout for lazy evaluation and parallel processing. Every function (`prepare_extracted_region`, `align_regional_kmers`, `reconstruct_counts`, `reconstruct_database`) spawns a distributed dask client.

**Reimplementation:** All dask dependencies have been removed entirely. The package has zero dask imports. Instead:

- **Database extraction** (`extract.py`): Degenerate expansion and trimming are applied directly as single pandas/numpy operations on the full dataset, rather than chunked dask.delayed calls.
- **K-mer alignment** (`align.py`): Replaced the `dask.delayed` + `itertools.product` loop with a fully vectorized numpy broadcast. The entire K×A mismatch matrix is computed in one operation: `(kmer_arr[:, np.newaxis, :] != asv_arr[np.newaxis, :, :]).sum(axis=2)`. This is both simpler and faster than the original approach for typical dataset sizes.
- **EM reconstruction** (`reconstruct.py`): The per-sample EM solver runs sequentially within a single Snakemake rule. For parallelism across samples, Snakemake checkpoints can split reconstruction into per-sample rules (handled at the workflow level, not inside this package).
- **Database building** (`build_database.py`): The three-round detangling algorithm uses eager pandas DataFrames throughout, replacing `dask.delayed` and `dd.from_delayed`.

**Rationale:** Dask's distributed scheduler conflicts with Snakemake's own parallelism management. Running both simultaneously causes resource contention and makes debugging difficult. Snakemake already provides DAG-based parallelism, so per-sample and per-region parallelism is better handled at the workflow level.

### 2. Removal of QIIME2 type system

**Original:** Functions accept and return QIIME2 types (`biom.Table`, `pd.Series` wrapped as `DNAFASTAFormat`, `qiime2.Metadata`).

**Reimplementation:** All functions accept and return file paths (strings) or standard pandas objects. Input/output is handled via TSV and FASTA files read with `pd.read_csv()` and BioPython's `SeqIO.parse()`.

### 3. Parallelism delegated to Snakemake

**Original:** Internal parallelism via dask workers.

**Reimplementation:** Each processing step (extract, trim, align, build database, reconstruct, taxonomy) is a separate Snakemake rule. Per-region work is parallelized via Snakemake wildcards (e.g., `{region}`). The EM solver can be further parallelized per-sample using Snakemake checkpoints if needed.

---

## Module-by-Module Changes

### `extract.py`

| Aspect | q2-sidle | Standalone |
|--------|----------|------------|
| Input | `DNAFASTAFormat` artifact | FASTA file path |
| Chunking | `dask.delayed` per chunk | Single pass over all sequences |
| Output | `DNAFASTAFormat` + `pd.DataFrame` | FASTA file + TSV file |

**Added function: `trim_dada2_posthoc()`** — Ports the trimming logic from `q2_sidle/_trim.py`. The original operates on `biom.Table` + `pd.Series` and uses `biom.Table.collapse()` to merge duplicate sequences after trimming. The standalone version reads FASTA + count TSV, trims, collapses duplicates via pandas `groupby` + `sum`, and writes output files. This is necessary because DADA2 merged paired-end reads produce variable-length ASVs that must be trimmed to uniform length before SIDLE alignment.

### `align.py`

| Aspect | q2-sidle | Standalone |
|--------|----------|------------|
| Pairwise comparison | `itertools.product` + `dask.delayed` per chunk pair | Vectorized numpy broadcast |
| Data structures | `dask.dataframe` with `.to_delayed()` | Plain numpy arrays |
| Memory model | Lazy evaluation via dask graph | Eager: full mismatch matrix in memory |

The vectorized approach computes the complete mismatch matrix `(K × A)` in a single numpy operation. For a typical dataset (e.g., 50,000 k-mers × 500 ASVs × 400 bp), this requires approximately 10 GB of memory for the boolean intermediate array — well within the capacity of a modern workstation. For very large datasets, Snakemake rules can split k-mers across multiple alignment jobs.

### `reconstruct.py`

| Aspect | q2-sidle | Standalone |
|--------|----------|------------|
| Count input | `biom.Table` | TSV file(s) via `pd.read_csv()` |
| Database map | `qiime2.Metadata` | TSV file via `pd.read_csv()` |
| Sample iteration | `dask.delayed` per sample | Sequential for loop (Snakemake parallelizes at rule level) |
| Output | `biom.Table` | TSV file |

The core EM algorithm (`_solve_ml_em_iterative_1_sample`) is a faithful port of the original SMURF algorithm. The E-step, M-step, convergence check, and low-abundance filtering are identical.

**Known issues still being addressed (Tasks 6–7):**
- Count loading currently sums all region count tables before reconstruction, losing per-region structure
- `_construct_align_mat()` uses `.loc` indexing that can fail on mismatched keys — needs to be replaced with `merge()` operations

### `build_database.py`

| Aspect | q2-sidle | Standalone |
|--------|----------|------------|
| Detangling | `dask.delayed` + `dd.from_delayed` for lazy pandas | Eager pandas DataFrames |
| Three-round algorithm | Identical logic | Identical logic |
| Output | `pd.DataFrame` + `qiime2.Metadata` | TSV files |

This is a complete port of `q2_sidle/_build_database.py` (~380 lines). The three-round database detangling algorithm is preserved exactly:

1. **Round 1:** Identify unique k-mers and their reference mappings
2. **Round 2:** Resolve shared k-mers using regional co-occurrence
3. **Round 3:** Final cleanup and name assignment

### `taxonomy.py`

No changes from the initial standalone port. Uses `database_params` from `utils.py` for taxonomy string parsing.

### `utils.py`

| Removed | Added |
|---------|-------|
| `import dask` | eHOMD entry in `database_params` |
| `from dask.distributed import Client` | — |
| `_setup_dask_client()` | — |
| `_chunks()` helper | — |

The `database_params` dictionary now includes `'ehomd'` alongside `'greengenes'`, `'silva'`, and `'none'`, supporting FAVABEAN's oral microbiome case studies that use the expanded Human Oral Microbiome Database.

### `cli.py`

No changes. Command-line interface preserved for standalone usage outside Snakemake.

---

## Dependency Changes

| Dependency | q2-sidle | Standalone |
|------------|----------|------------|
| QIIME2 | Required | Not needed |
| biom-format | Required | Not needed |
| dask | Required | Not needed |
| dask.distributed | Required | Not needed |
| numpy | Required | Required |
| pandas | Required | Required |
| biopython | Not used (uses skbio) | Required (FASTA I/O) |
| scikit-bio | Required | Required (DNA degeneracy) |

The total dependency footprint is reduced from ~50+ packages (full QIIME2 environment) to 4 packages.

---

## Integration with FAVABEAN

The standalone SIDLE package is integrated into FAVABEAN via a separate Snakemake rules file (`workflow/rules/sidle.smk`) that chains the modules together:

```
DADA2 condense output
    │
    ├── seqtab_to_fasta.R  (convert condensed TSV → FASTA + counts)
    │
    ├── prepare_extracted_region  (per region: extract k-mers from reference)
    │
    ├── trim_dada2_posthoc  (per region: trim ASVs to uniform length)
    │
    ├── align_regional_kmers  (per region: vectorized k-mer alignment)
    │
    ├── reconstruct_database  (all regions: build database map + summary)
    │
    ├── reconstruct_counts  (all regions: EM reconstruction)
    │
    ├── reconstruct_taxonomy  (assign taxonomy to reconstructed features)
    │
    └── BIOM conversion  (for downstream FALAPhyl analysis)
```

Users enable SIDLE by setting `use_sidle: true` in `favabean.yaml`. The original primer averaging approach (`primer_average.R`) is retained for backward compatibility but is no longer the recommended method for multi-primer datasets.

---

## Citation

If using this reimplementation, please cite:

- Debelius, J. et al. (2021). "SMURF: a Bayesian approach to reconstruction of 16S rRNA sequences from multiple regions." bioRxiv 2021.03.23.436606.
- Fuks, G. et al. (2018). "Combining 16S rRNA gene variable regions enables high-resolution microbial community profiling." Microbiome 6:17.
