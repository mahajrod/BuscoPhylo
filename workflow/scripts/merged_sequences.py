#!/usr/bin/env python3
__author__ = 'tomarovsky'

from pathlib import Path
import argparse


def main():
    common_ids_file = open(args.input, 'r')
    outdir = Path(args.outdir)
    outdir.mkdir()
    for idname in common_ids_file:
        outfile = open(outdir / "merged_" + idname + ".faa", 'a')
        for directory in args.sequence_directories:
            dirpath = Path(directory)
            species_name = dirpath.parents[1].split("/")[-1]
            header = ">" + species_name
            with open(dirpath / idname + ".faa", 'r') as f:
                seq = f.readlines()[1]
            outline = "\n".join([header, seq])
            outfile.write(outline)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="script for merging files with sequences of different species")
    group_required = parser.add_argument_group('Required options')
    group_required.add_argument('-i', '--input', type=str, help="single_copy_busco_sequences.uniq.ids file")
    group_required.add_argument('-d', '--species-directories', type=float, nargs='+',
                                help="list of directories with sequence files")
    group_required.add_argument('-o', '--outdir', type=str, help="output directory name")
    args = parser.parse_args()
    main()