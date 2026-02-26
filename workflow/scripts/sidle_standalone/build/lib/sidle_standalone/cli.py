#!/usr/bin/env python3
"""
SIDLE CLI - Command-line interface for standalone SIDLE

Usage:
  sidle extract <reference.fasta> <region> <trim_length> <fwd_primer> <rev_primer> <out_kmers.fasta> <out_map.tsv>
  sidle trim <rep_seqs.fasta> <counts.tsv> <out.fasta> <out_counts.tsv> [--trim-length N]
  sidle align <kmers.fasta> <asvs.fasta> <region> <out_alignment.tsv>
  sidle build-db <regions.txt> <alignments.txt> <kmer_maps.txt> <out_map.tsv> <out_summary.tsv>
  sidle reconstruct <regions.txt> <alignments.txt> <counts.txt> <db_map.tsv> <db_summary.tsv> <out_counts.tsv>
  sidle taxonomy <recon_map.tsv> <taxonomy.tsv> <out_taxonomy.tsv> [--database silva]
"""
import sys
import sidle_standalone as sidle


def main():
    if len(sys.argv) < 2:
        print(__doc__)
        sys.exit(1)
    
    cmd = sys.argv[1]
    
    if cmd == 'extract':
        if len(sys.argv) < 9:
            print("Usage: sidle extract <ref.fasta> <region> <trim_len> <fwd> <rev> <out_kmers.fa> <out_map.tsv>")
            sys.exit(1)

        sidle.prepare_extracted_region(
            sequences_file=sys.argv[2],
            region=sys.argv[3],
            trim_length=int(sys.argv[4]),
            fwd_primer=sys.argv[5],
            rev_primer=sys.argv[6],
            kmer_output=sys.argv[7],
            map_output=sys.argv[8]
        )
        print(f"✓ Prepared region {sys.argv[3]}")
    
    elif cmd == 'trim':
        if len(sys.argv) < 6:
            print("Usage: sidle trim <rep_seqs.fa> <counts.tsv> <out.fa> <out_counts.tsv> [--trim-length N]")
            sys.exit(1)

        trim_length = 0
        if '--trim-length' in sys.argv:
            trim_length = int(sys.argv[sys.argv.index('--trim-length') + 1])

        sidle.trim_dada2_posthoc(
            rep_seqs_file=sys.argv[2],
            count_file=sys.argv[3],
            output_fasta=sys.argv[4],
            output_counts=sys.argv[5],
            trim_length=trim_length
        )
        print("✓ Trimmed ASVs")

    elif cmd == 'align':
        if len(sys.argv) < 6:
            print("Usage: sidle align <kmers.fa> <asvs.fa> <region> <out_align.tsv>")
            sys.exit(1)

        sidle.align_regional_kmers(
            kmers_file=sys.argv[2],
            rep_seq_file=sys.argv[3],
            region=sys.argv[4],
            output_file=sys.argv[5]
        )
        print(f"✓ Aligned region {sys.argv[4]}")

    elif cmd == 'build-db':
        if len(sys.argv) < 7:
            print("Usage: sidle build-db <regions.txt> <alignments.txt> <kmer_maps.txt> <out_map.tsv> <out_summary.tsv>")
            sys.exit(1)

        regions = [line.strip() for line in open(sys.argv[2])]
        aligns = [line.strip() for line in open(sys.argv[3])]
        kmer_maps = [line.strip() for line in open(sys.argv[4])]

        sidle.reconstruct_database(
            region=regions,
            alignment_files=aligns,
            kmer_map_files=kmer_maps,
            output_map_file=sys.argv[5],
            output_summary_file=sys.argv[6]
        )
        print(f"✓ Built database map from {len(regions)} regions")

    elif cmd == 'reconstruct':
        if len(sys.argv) < 8:
            print("Usage: sidle reconstruct <regions.txt> <aligns.txt> <counts.txt> <map.tsv> <summary.tsv> <out.tsv>")
            sys.exit(1)
        
        regions = [line.strip() for line in open(sys.argv[2])]
        aligns = [line.strip() for line in open(sys.argv[3])]
        counts = [line.strip() for line in open(sys.argv[4])]
        
        sidle.reconstruct_counts(
            region=regions,
            alignment_files=aligns,
            count_files=counts,
            database_map_file=sys.argv[5],
            database_summary_file=sys.argv[6],
            output_file=sys.argv[7]
        )
        print(f"✓ Reconstructed counts from {len(regions)} regions")
    
    elif cmd == 'taxonomy':
        if len(sys.argv) < 5:
            print("Usage: sidle taxonomy <recon_map.tsv> <taxonomy.tsv> <out_taxonomy.tsv> [--database silva]")
            sys.exit(1)
        
        database = 'none'
        if '--database' in sys.argv:
            database = sys.argv[sys.argv.index('--database') + 1]
        
        sidle.reconstruct_taxonomy(
            reconstruction_map_file=sys.argv[2],
            taxonomy_file=sys.argv[3],
            output_file=sys.argv[4],
            database=database
        )
        print("✓ Reconstructed taxonomy")
    
    else:
        print(f"Unknown command: {cmd}")
        print(__doc__)
        sys.exit(1)


if __name__ == '__main__':
    main()
