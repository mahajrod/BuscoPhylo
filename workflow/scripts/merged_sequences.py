#!/usr/bin/env python3
from pathlib import Path
import argparse


def main():
    common_ids_file = open(args.common_ids, 'r')
    outdir = Path(args.outdir)
    outdir.mkdir()
    for idname in common_ids_file:
        idname = idname.rstrip()
        file_faa, file_fna = f"{idname}.faa", f"{idname}.fna"
        merged_file_faa, merged_file_fna = f"{idname}.merged.faa", f"{idname}.merged.fna"
        out_faa, out_fna = open(outdir / merged_file_faa, 'a'), open(outdir / merged_file_fna, 'a')
        for directory in args.single_copy_files:
            dirpath = Path(directory)
            header = ">%s" % dirpath.parents[0].stem
            with open(dirpath / file_faa, 'r') as f:
                seq_faa = "%s\n" % f.readlines()[1].strip()
                outline_faa = "\n".join([header, seq_faa])
                out_faa.write(outline_faa)
            with open(dirpath / file_fna, 'r') as f:
                seq_fna = "%s\n" % f.readlines()[1].strip()
                outline_fna = "\n".join([header, seq_fna])
                out_fna.write(outline_fna)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="script for merging files with sequences of different species")
    group_required = parser.add_argument_group('Required options')
    group_required.add_argument('-c', '--common_ids', type=str,
                                help="common_ids file")
    group_required.add_argument('-s', '--single_copy_files', type=str, nargs='+',
                                help="single copy busco sequences directories")
    group_required.add_argument('-o', '--outdir', type=str, help="output directory name")
    args = parser.parse_args()
    main()