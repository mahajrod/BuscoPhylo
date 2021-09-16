localrules: gblocks

rule gblocks:
    input:
        # fna=mafft_dir_path / "{sample}.fna"
        mafft_dir=directory(mafft_dir_path),
        names_dir=directory(output_dir_path / "tmp" / "{sample}")
    output:
        # gb=gblocks_dir_path / "{sample}.fna-gb",
        # gb_txt=gblocks_dir_path / "{sample}.fna-gb.txt"
        temp(directory(output_dir_path / "gblocks_tmp" / "{sample}"))
    params:
        gblocks_path=config["gblocks_path"],
        gblocks_dna_flags="-t=d -p=t",
        gblocks_protein_flags="-t=p -p=t"
    log:
        std=log_dir_path / "{sample}.gblocks.log",
        cluster_log=cluster_log_dir_path / "{sample}.gblocks.cluster.log",
        cluster_err=cluster_log_dir_path / "{sample}.gblocks.cluster.err"
    benchmark:
        benchmark_dir_path / "{sample}.gblocks.benchmark.txt"
    # conda:
    #     "../../%s" % config["conda_config"]
    resources:
        cpus=config["gblocks_threads"],
        time=config["gblocks_time"],
        mem=config["gblocks_mem_mb"]
    shell:
        "sleep 10; "
        "mkdir -p {output}; "
        "for FILE in `ls -d {input.names_dir}/*`; do "
        "FILE=$(basename $FILE); "
        "{params.gblocks_path}/Gblocks {input.mafft_dir}/merged_$FILE.fna {params.gblocks_dna_flags} 1> {log.std} || true; "
        "{params.gblocks_path}/Gblocks {input.mafft_dir}/merged_$FILE.faa {params.gblocks_protein_flags} 1> {log.std} || true; "
        "sleep 10; "
        "mv {input.mafft_dir}/merged_$FILE.fna-gb {output}/; mv {input.mafft_dir}/merged_$FILE.fna-gb.txt {output}/; "
        "mv {input.mafft_dir}/merged_$FILE.faa-gb {output}/; mv {input.mafft_dir}/merged_$FILE.faa-gb.txt {output}/; "
        "done"