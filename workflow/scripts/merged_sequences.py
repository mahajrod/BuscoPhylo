#!/usr/bin/env python3
from pathlib import Path
import argparse


def main():
    common_ids_file = open(args.input, 'r')
    outdir = Path(args.outdir)
    outdir.mkdir()
    for idname in common_ids_file:
        idname = idname.strip()
        file_faa, file_fna = idname + ".faa", idname + ".fna"
        merged_file_faa, merged_file_fna = idname + ".merged.faa", idname + ".merged.fna"
        out_faa, out_fna = open(outdir / merged_file_faa, 'a'), open(outdir / merged_file_fna, 'a')
        for directory in args.single_copy_files:
            dirpath = Path(directory)
            header = ">" + str(dirpath.parents[0].stem)
            with open(dirpath / file_faa, 'r') as f:
                seq_faa = f.readlines()[1].strip() + "\n"
            with open(dirpath / file_fna, 'r') as f:
                seq_fna = f.readlines()[1].strip() + "\n"
            outline_faa, outline_fna = "\n".join([header, seq_faa]), "\n".join([header, seq_fna])
            out_faa.write(outline_faa)
            out_fna.write(outline_fna)



if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="script for merging files with sequences of different species")
    group_required = parser.add_argument_group('Required options')
    group_required.add_argument('-i', '--input', type=str,
                                help="single_copy_busco_sequences.uniq.ids file")
    group_required.add_argument('-s', '--single_copy_files', type=str, nargs='+',
                                help="list of single_copy_busco_sequences files")
    group_required.add_argument('-o', '--outdir', type=str, help="output directory name")
    args = parser.parse_args()
    main()