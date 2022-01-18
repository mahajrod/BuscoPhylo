#!/usr/bin/env python3

from collections import defaultdict
from Bio import SeqIO
from sys import stdin
import argparse


def main():
    sequence_map = defaultdict(str)
    if args.input is not stdin:
        for i in args.input:
            seq_length = None
            for sequence in SeqIO.parse(i, "fasta"):
                sequence_map[sequence.name] += str(sequence.seq)
                if seq_length != len(str(sequence.seq)) and seq_length is not None:
                    print(i, sequence)
                seq_length = len(str(sequence.seq))
    else:
        for sequence in SeqIO.parse(args.input, "fasta"):
            sequence_map[sequence.name] += str(sequence.seq)
            if seq_length != len(str(sequence.seq)) and seq_length is not None:
                print("!!!!!!!!!!!!!!!!!!!!!!!!!", sequence.name)
            seq_length = len(str(sequence.seq))

    outfile = open(args.output, "a")
    for key, value in sequence_map.items():
        outfile.write(f">{key}\n{value}")
    outfile.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="script for concatenate FASTA files into one big FASTA file "
                                                 "by concatenating sequences with the same identifier")
    group_required = parser.add_argument_group('Required options')
    group_required.add_argument('-i', '--input', type=str, default=stdin,
                                nargs="+", help="input concat FASTA file with the same headers or stdin")
    group_required.add_argument('-o', '--output', type=str, help="output FASTA file name")
    args = parser.parse_args()
    main()
