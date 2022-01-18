#!/usr/bin/env python3

from pathlib import Path
import argparse

def get_id(header):
    id = header[1:].split("|")[0]
    if "_" in id:
        id = id.split("_")[0]
    return id

def split_fasta_to_directory(fasta, outdir, single_copy_ids, ext):
    opened_file_flag = False
    single_copy_file_flag = True
    with open(fasta, "r") as infile:
        for line in infile:
            if line[0] == ">":
                if opened_file_flag:
                    outfile.close()
                opened_file_flag = True
                id = get_id(line)
                if id not in single_copy_ids:
                    single_copy_file_flag = False
                    opened_file_flag = False
                    continue
                else:
                    single_copy_file_flag = True
                filepath = Path(outdir / f"{id}.{ext}")
                outfile = open(filepath, "w")
            if single_copy_file_flag:
                outfile.write(line)
        outfile.close()


def main():
    # set PATHs
    initial_codon_fasta = list(Path(args.initial_results).glob("*.codon.fas"))[0]
    initial_protein_fasta = list(Path(args.initial_results).glob("*a.fas"))[0] # *.fasta.fas or *.fa.fas
    rerun_codon_fasta = list(Path(args.rerun_results).glob("*.codon.fas"))[0]
    rerun_protein_fasta = list(Path(args.rerun_results).glob("*a.fas"))[0] # *.fasta.fas or *.fa.fas
    single_copy_ids = [id.stem for id in list(Path(args.single_copy_busco_sequences).glob("*.faa"))]
    print("Single copy IDs: ", len(single_copy_ids))
    outdir = Path(args.outdir)
    outdir.mkdir()
    # write FASTA files with single copy CDs sequences to output directory
    split_fasta_to_directory(initial_codon_fasta, outdir, single_copy_ids, "fna")
    split_fasta_to_directory(initial_protein_fasta, outdir, single_copy_ids, "faa")
    split_fasta_to_directory(rerun_codon_fasta, outdir, single_copy_ids, "fna")
    split_fasta_to_directory(rerun_protein_fasta, outdir, single_copy_ids, "faa")
    single_copy_CDs = [cd.stem for cd in list(outdir.glob("*.faa"))]
    print(len("Single copy CDs: ", single_copy_CDs))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="script to get directory with single copy CDs sequences from Metaeuk output")
    group_required = parser.add_argument_group('Required options')
    group_required.add_argument('--initial_results', type=str,
                                help="initial_results directory path")
    group_required.add_argument('--rerun_results', type=str,
                                help="rerun_results directory path")
    group_required.add_argument('--single_copy_busco_sequences', type=str,
                                help="single_copy_busco_sequences directory path")
    group_required.add_argument('--outdir', type=str, help="single copy CDs directory path")
    args = parser.parse_args()
    main()