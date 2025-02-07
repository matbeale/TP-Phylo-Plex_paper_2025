#!/usr/bin/env python

import random
import argparse
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def mutate_sequence(sequence, mutation_percentage):
    """
    Introduces random mutations to a specified percentage of sites in a sequence.
    
    Parameters:
        sequence (str): The original sequence.
        mutation_percentage (float): The percentage of sites to mutate (0-100).
        
    Returns:
        str: The mutated sequence.
    """
    sequence = list(sequence)
    num_mutations = int(len(sequence) * mutation_percentage / 100)
    positions = random.sample(range(len(sequence)), num_mutations)
    
    for pos in positions:
        original_base = sequence[pos]
        new_base = random.choice([b for b in "ACGT" if b != original_base])
        sequence[pos] = new_base
    
    return "".join(sequence)

def main():
    # Set up argument parser
    parser = argparse.ArgumentParser(description="Generate multiple mutated versions of a genome sequence.")
    parser.add_argument("input_file", type=str, help="Path to the input FASTA file.")
    parser.add_argument("output_file", type=str, help="Path to the output multi-FASTA file.")
    parser.add_argument("mutation_percentage", type=float, help="Percentage of sites to mutate (0-100).")
    parser.add_argument("num_mutations", type=int, help="Number of mutated versions to generate.")
    args = parser.parse_args()

    # Validate mutation percentage
    if not (0 <= args.mutation_percentage <= 100):
        parser.error("The mutation percentage must be between 0 and 100.")
    
    # Read the input genome sequence
    record = SeqIO.read(args.input_file, "fasta")
    original_sequence = str(record.seq)

    # Generate multiple mutated versions
    mutated_records = []
    for i in range(args.num_mutations):
        mutated_sequence = mutate_sequence(original_sequence, args.mutation_percentage)
        mutated_record = SeqRecord(
            Seq(mutated_sequence),
            id=f"{record.id}_mutated_{i+1}",
            description=f"Mutated version {i+1}"
        )
        mutated_records.append(mutated_record)

    # Save the mutated sequences to a multi-FASTA file
    SeqIO.write(mutated_records, args.output_file, "fasta")

    print(f"{args.num_mutations} mutated sequences saved to {args.output_file}")

if __name__ == "__main__":
    main()

