#!/usr/bin/env python3

from Bio import AlignIO
from sys import stdin
import argparse


def main():
    infile = args.input
    outfile = open(args.output, "w")
    AlignIO.convert(infile, "fasta", outfile, "nexus", args.type)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="script for converting FASTA format to NEXUS format")
    group_required = parser.add_argument_group('Required options')
    group_required.add_argument('-i', '--input', type=str, default=stdin, help="input concat FASTA file or stdin")
    group_required.add_argument('-t', '--type', type=str, help="molecular type (DNA, RNA or protein)")
    group_required.add_argument('-o', '--output', type=str, help="output NEXUS file name")
    args = parser.parse_args()
    main()
