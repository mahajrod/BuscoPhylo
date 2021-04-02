#!/usr/bin/env python3
__author__ = 'tomarovsky'

from pathlib import Path
import argparse


def main():
    outdir = Path(args.outdir)
    outdir.mkdir()
    outfile_path = outdir / "mafft.tasks.0.sh"

    if args.file_extension == "faa":
        merged_files = Path(args.input).glob("*.faa")
        mafft_faa_command = "mafft --anysymbol {input} > {output}"
    elif args.file_extension == "fna":
        merged_files = Path(args.input).glob("*.fna")
        mafft_faa_command = "mafft {input} > {output}"
    else:
        print ("your extension is wrong!")

    counter = 0
    for file in merged_files:
        counter += 1
        mafft_command_outname = file.stem + ".mafft." + args.file_extension
        mafft_command_output = args.mafft_command_outdir / mafft_command_outname
        if counter % args.amount == 0:
            outfile = "mafft.tasks.%s.sh" % str(counter)
            outfile_path = outdir / outfile
        with open(outfile_path, 'a') as out:
            out.write(mafft_faa_command.format(input=file, output=mafft_command_output) + "\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="script for getting bash scripts for mafft multiple alignment")
    group_required = parser.add_argument_group('Required options')
    group_required.add_argument('-i', '--input', type=str,
                                help="merged_sequences directory (fna and faa files)")
    group_required.add_argument('-o', '--outdir', type=str, help="output directory name")
    group_required.add_argument('-f', '--file-extension', type=str, help="'faa' or 'fna'", default='faa')
    group_additional = parser.add_argument_group('Additional options')
    group_additional.add_argument('-a', '--amount', type=int, 
                             default=20, help="amount of mafft commands in the file")
    group_additional.add_argument('-m', '--mafft_command_outdir', type=Path, help="output directory in maft command")
    args = parser.parse_args()
    main()