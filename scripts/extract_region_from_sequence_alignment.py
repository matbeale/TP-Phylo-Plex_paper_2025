#!/usr/bin/env python
# -*- coding: utf-8 -*-

from Bio import SeqIO
from Bio import AlignIO
from Bio.Align import AlignInfo
from Bio.Align import MultipleSeqAlignment
from Bio.Seq import Seq
#from Bio.Alphabet import IUPAC
import pandas as pd
import argparse

description= "Takes a multiple sequence alignment and a csv/bed file, and extracts a subset of positions, creating one new alignment fasta per line. \n"


parser = argparse.ArgumentParser(description=description)
parser.add_argument("seq_alignment_file", help="Fasta file containing alignment")
parser.add_argument("-c", type=str, dest="csv", help="Comma separated region positions (<start>,<end>,<name>)" )
parser.add_argument("-g", type=str, dest="gene", help="Comma separated region positions with orientation (+/-) to allow reverse complementing (<start>,<end>,<name>,<+/->)" )
parser.add_argument("-b", type=str, dest="bed", help="BED file - tab separated region positions (<chrom>\t<start>\t<end>\t<name>)" )



args = parser.parse_args()
align = AlignIO.read(args.seq_alignment_file, 'fasta')


# Extract forward strand sequence and relabel header
def get_subaln(full_aln, mystart, myend, gene_name):
    subalign = full_aln[:, mystart:myend]
    fw_align = MultipleSeqAlignment([])
    for record in subalign:
        myseq = record.seq
        fw_align.add_sequence(record.id + "::" + gene_name + "::+", str(myseq))
    return(fw_align)

# Extract rev strand sequence, reverse complement and relabel header
def get_subaln_rc(full_aln, mystart, myend, gene_name):
    subalign = full_aln[:, mystart:myend]
    rc_align = MultipleSeqAlignment([])
    for record in subalign:
        rev_seq = record.seq.reverse_complement()
        rc_align.add_sequence(record.id + "::" + gene_name + "::-", str(rev_seq))
    return(rc_align)


if args.csv :
    csv_file = pd.read_csv(args.csv, header=None)
    for myline in range(0, len(csv_file.index)):
        subalign = get_subaln(align, int(csv_file.iloc[myline,0])-1, int(csv_file.iloc[myline,1]), csv_file.iloc[myline,2])
        SeqIO.write(subalign, csv_file.iloc[myline,2] + ".fa", 'fasta')


if args.gene :
    gene_file = pd.read_csv(args.gene, header=None)
    for myline in range(0, len(gene_file.index)):
        if  gene_file.iloc[myline,3] == '+':
            subalign = get_subaln(align, int(gene_file.iloc[myline,0])-1, int(gene_file.iloc[myline,1]), gene_file.iloc[myline,2])
            SeqIO.write(subalign, gene_file.iloc[myline,2] + ".fa", 'fasta')
        elif gene_file.iloc[myline,3] == '-':
            subalign = get_subaln_rc(align, int(gene_file.iloc[myline,0])-1, int(gene_file.iloc[myline,1]), gene_file.iloc[myline,2])
            SeqIO.write(subalign, gene_file.iloc[myline,2] + ".fa", 'fasta')
        else:
            print("Error in " + gene_file.iloc[myline,2] + ". Need to specify orientation as + or -")
        #print(subalign)
#    SeqIO.write(subalign, gene_file.iloc[myline,2] + ".fa", 'fasta')



if args.bed :
    bed_file = pd.read_csv(args.bed, sep='\t', header=None)
    for myline in range(0, len(bed_file.index)):
        subalign = get_subaln(align, int(bed_file.iloc[myline,1]), int(bed_file.iloc[myline,2])-1, bed_file.iloc[myline,3])
        SeqIO.write(subalign, bed_file.iloc[myline,2] + ".fa", 'fasta')

