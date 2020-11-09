#!/usr/bin/env python3

import sys, time, argparse, subprocess, os, signal, random, string
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.Seq import MutableSeq
from Bio.Alphabet import generic_dna

Description = """
Testing script
"""

#------------------------------------------------------------
# Get n sequences from the input fasta file
def filter(input_file):
    filename, file_extension = os.path.splitext(input_file)
    input_file_iterator = SeqIO.parse(input_file, "fasta", generic_dna)
    output_path = filename + "." + "filtered" + file_extension
    output_file = open(output_path, "w")

    with open(input_file, "r") as handle:
        for record in SeqIO.parse(handle, "fasta", generic_dna):
            seq = "".join(list([c for c in str(record.seq) if (c in string.ascii_letters)]))
            seq = seq + "$"
            out_record = SeqRecord(Seq(seq, generic_dna), record.id, description="filtered")
            SeqIO.write(out_record, output_file, "fasta")

    output_file.close()

#------------------------------------------------------------
def main():
    parser = argparse.ArgumentParser(description=Description, formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('-i', '--input', help='input file name', type=str, dest="input")

    args = parser.parse_args()

    filter(args.input)


if __name__ == '__main__':
    main()