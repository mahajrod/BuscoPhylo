#!/usr/bin/env python3

from collections import defaultdict
from Bio import SeqIO
from sys import stdin
import argparse


def main():
    sequence_map = defaultdict(str)
    seq_length = None
    if args.input is not stdin:
        for i in args.input:
            for sequence in SeqIO.parse(i, "fasta"):
                sequence_map[sequence.name] += str(sequence.seq).upper()
    else:
        for sequence in SeqIO.parse(args.input, "fasta"):
            sequence_map[sequence.name] += str(sequence.seq).upper()
            if seq_length != len(str(sequence.seq)) and seq_length is not None:
                print("!!!!!!!!!!!!!!!!!!!!!!!!!", i)
            seq_length = len(str(sequence.seq))
    # for k, v in sequence_map.items():
    #     print(k, len(v), sep="\t")

    outfile = open(args.output, "w")
    for key in sequence_map.keys():
        outfile.write('>' + key + '\n')
        outfile.write(sequence_map[key] + '\n')
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
