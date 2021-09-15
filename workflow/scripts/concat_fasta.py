#!/usr/bin/env python3

from Bio import SeqIO
from collections import defaultdict


def main():
    sequence_map = defaultdict(str)

    for sequence in SeqIO.parse('all.fasta', "fasta"):
        sequence_map[sequence.name] += str(sequence.seq)

    for key in sorted(sequence_map.iterkeys()):
        print('>' + key)
        print(sequence_map[key])

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="script for concatenate FASTA files into one big FASTA file "
                                                 "by concatenating sequences with the same identifier")
    group_required = parser.add_argument_group('Required options')
    group_required.add_argument('-i', '--input', type=str, help="FASTA files directory")
    group_required.add_argument('-s', '--single_copy_files', type=str, nargs='+',
                                help="list of single_copy_busco_sequences files")
    group_required.add_argument('-o', '--output', type=str, help="output FASTA file name")
    args = parser.parse_args()
    main()