#!/usr/bin/env python

import argparse
from Bio import SeqIO
from collections import defaultdict

# written with assistance from ChatGTP 3.5


def read_input_file(input_file):
    with open(input_file, 'r') as f:
        fasta_files = [line.strip() for line in f if line.strip()]
    return fasta_files

def concatenate_alignments(fasta_files, output_file):
    # Dictionary to hold concatenated sequences for each sample
    concatenated_sequences = defaultdict(str)
    
    # Process each fasta file
    for fasta_file in fasta_files:
        for record in SeqIO.parse(fasta_file, "fasta"):
            concatenated_sequences[record.id] += str(record.seq)
    
    # Write the concatenated sequences to the output file
    with open(output_file, "w") as output_handle:
        for sample_id, sequence in concatenated_sequences.items():
            output_handle.write(f">{sample_id}\n{sequence}\n")

def main():
    # Set up argument parsing
    parser = argparse.ArgumentParser(description="Concatenate multiple sequence alignments from FASTA files.")
    parser.add_argument("-i", "--input", required=True, help="File containing list of FASTA files to concatenate")
    parser.add_argument("-o", "--output", required=True, help="Output file for concatenated sequences")

    args = parser.parse_args()

    # Read the list of FASTA files from the input file
    fasta_files = read_input_file(args.input)
    
    # Concatenate the alignments
    concatenate_alignments(fasta_files, args.output)

if __name__ == "__main__":
    main()


