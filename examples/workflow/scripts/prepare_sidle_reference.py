#!/usr/bin/env python3
"""
Prepare a SIDLE-ready reference database from DADA2 training FASTAs.

DADA2 taxonomy assignment uses two reference files:
  - ref FASTA:     headers = taxonomy string  (Kingdom;Phylum;...;Genus;)
  - species FASTA: headers = "AccessionID Genus species"

These contain the SAME sequences with different header formats.  This script
merges them into a single FASTA with unique, species-enriched headers and
creates a taxonomy lookup TSV for SIDLE's taxonomy reconstruction.

For databases where the ref FASTA already includes species (e.g., SILVA
wSpecies), or where no species file is provided, the script still ensures
unique headers so that no reference sequences are silently dropped by
read_fasta()'s dict-based storage.

Usage:
  python prepare_sidle_reference.py <ref.fa.gz> <species.fa.gz|none> \
      <out_prepared.fasta> <out_taxonomy.tsv>
"""

import gzip
import os
import sys
from collections import Counter


def parse_fasta(path):
    """Parse FASTA (optionally gzipped) into list of (header, sequence) tuples.

    Normalizes spaces after semicolons in headers so that taxonomy-as-ID
    formats like eHOMD ('; species') or Greengenes ('; p__') are consistent.
    """
    opener = gzip.open if path.endswith('.gz') else open
    entries = []
    with opener(path, 'rt') as f:
        header = None
        seq_parts = []
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith('>'):
                if header is not None:
                    entries.append((header, ''.join(seq_parts)))
                header = line[1:].replace('; ', ';')
                seq_parts = []
            else:
                seq_parts.append(line)
        if header is not None:
            entries.append((header, ''.join(seq_parts)))
    return entries


def prepare_reference(ref_path, species_path, out_fasta, out_taxonomy):
    """Merge ref + species FASTAs into a SIDLE-ready reference.

    Strategy
    --------
    When a species file is available and its sequences match the ref file
    (positionally), we use accession numbers from the species file as unique
    IDs and build full taxonomy strings (Kingdom;...;Genus;Species).

    When no species file is available (or sequences don't match), we
    deduplicate ref headers by appending a counter suffix (#N) to repeats,
    preserving every reference sequence.

    Parameters
    ----------
    ref_path : str
        Path to DADA2 training-set FASTA (taxonomy-as-header, optionally gzipped)
    species_path : str
        Path to DADA2 species-assignment FASTA, or "none"/empty to skip
    out_fasta : str
        Output path for prepared FASTA (unique headers, species-enriched)
    out_taxonomy : str
        Output path for taxonomy TSV (Feature ID → Taxon)
    """
    ref_entries = parse_fasta(ref_path)
    n_ref = len(ref_entries)
    n_unique_orig = len(set(h for h, _ in ref_entries))

    print(f"Reference FASTA: {n_ref} sequences, {n_unique_orig} unique headers")

    # ------------------------------------------------------------------
    # Try to load and match species file
    # ------------------------------------------------------------------
    use_species = False
    spec_accessions = []
    spec_species = []

    if species_path and species_path.lower() != 'none' and os.path.exists(species_path):
        spec_entries = parse_fasta(species_path)

        if len(spec_entries) == n_ref:
            # Verify sequences match positionally
            match = all(r[1] == s[1] for r, s in zip(ref_entries, spec_entries))
            if match:
                use_species = True
                for header, _ in spec_entries:
                    parts = header.split(' ', 1)
                    spec_accessions.append(parts[0])
                    spec_species.append(parts[1] if len(parts) > 1 else '')
                print(f"Species FASTA: {len(spec_entries)} sequences matched positionally")
            else:
                print("WARNING: Species FASTA sequences do not match ref — "
                      "falling back to ref-only mode")
        else:
            print(f"WARNING: Species FASTA has {len(spec_entries)} sequences "
                  f"vs {n_ref} in ref — falling back to ref-only mode")

    # ------------------------------------------------------------------
    # Build entries with unique IDs and enriched taxonomy
    # ------------------------------------------------------------------
    final = []  # list of (unique_id, taxonomy_string, sequence)

    if use_species:
        # Use accession numbers as unique IDs; enrich taxonomy with species
        # Guard against non-unique accessions (rare, but possible in custom DBs)
        accession_counts = Counter(spec_accessions)
        accession_seen = {}

        for i, (ref_header, seq) in enumerate(ref_entries):
            accession = spec_accessions[i]
            species = spec_species[i]
            taxonomy = ref_header.rstrip(';')

            # Count ranks to decide if species should be appended
            ranks = [r.strip() for r in taxonomy.split(';') if r.strip()]
            if len(ranks) < 7 and species:
                taxonomy = taxonomy + ';' + species

            # Ensure unique ID even if accessions repeat
            if accession_counts[accession] > 1:
                accession_seen.setdefault(accession, 0)
                accession_seen[accession] += 1
                unique_id = f"{accession}#seq{accession_seen[accession]}"
            else:
                unique_id = accession

            final.append((unique_id, taxonomy, seq))

        n_with_species = sum(1 for _, t, _ in final
                             if len([r for r in t.split(';') if r.strip()]) >= 7)
        n_deduped_acc = sum(1 for a, c in accession_counts.items() if c > 1)
        print(f"Species-enriched entries: {n_with_species}/{n_ref}")
        if n_deduped_acc > 0:
            print(f"WARNING: {n_deduped_acc} non-unique accession(s) disambiguated with #seqN suffix")

    else:
        # No species file — deduplicate headers with counter suffix
        header_counts = Counter(h for h, _ in ref_entries)
        header_seen = {}

        for header, seq in ref_entries:
            taxonomy = header.rstrip(';')

            if header_counts[header] > 1:
                header_seen.setdefault(header, 0)
                header_seen[header] += 1
                unique_id = f"{taxonomy}#seq{header_seen[header]}"
            else:
                unique_id = taxonomy

            final.append((unique_id, taxonomy, seq))

        n_deduped = sum(1 for h, c in header_counts.items() if c > 1)
        if n_deduped > 0:
            print(f"Deduplicated {n_deduped} non-unique headers "
                  f"({n_ref - n_unique_orig} sequences had duplicate headers)")

    # ------------------------------------------------------------------
    # Write outputs
    # ------------------------------------------------------------------
    os.makedirs(os.path.dirname(out_fasta) or '.', exist_ok=True)
    os.makedirs(os.path.dirname(out_taxonomy) or '.', exist_ok=True)

    with open(out_fasta, 'w') as f:
        for uid, _, seq in final:
            f.write(f'>{uid}\n{seq}\n')

    with open(out_taxonomy, 'w') as f:
        for uid, tax, _ in final:
            f.write(f'{uid}\t{tax}\n')

    n_unique_final = len(set(uid for uid, _, _ in final))
    print(f"Wrote {len(final)} sequences ({n_unique_final} unique IDs) "
          f"to {out_fasta}")
    print(f"Wrote taxonomy TSV to {out_taxonomy}")


if __name__ == '__main__':
    if len(sys.argv) < 5:
        print("Usage: prepare_sidle_reference.py <ref.fa.gz> "
              "<species.fa.gz|none> <out.fasta> <out_taxonomy.tsv>")
        sys.exit(1)

    prepare_reference(
        ref_path=sys.argv[1],
        species_path=sys.argv[2],
        out_fasta=sys.argv[3],
        out_taxonomy=sys.argv[4]
    )
