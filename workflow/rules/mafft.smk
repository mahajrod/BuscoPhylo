localrules: merged_sequences, mafft_dna, mafft_protein, finish


checkpoint merged_sequences:
    input:
        common_ids=output_dir_path / "single_copy_busco_sequences.common.dir/single_copy_busco_sequences.common.ids{N}"
    output:
        merged_ids=directory(output_dir_path / "merged_sequences/{N}")
    params:
        single_copy_files=expand(busco_dir_path / "{species}" / "single_copy_busco_sequences", species=config["species_list"])
    log:
        std=log_dir_path / "{N}.merged_ids.log",
        cluster_log=cluster_log_dir_path / "{N}.merged_ids.cluster.log",
        cluster_err=cluster_log_dir_path / "{N}.merged_ids.cluster.err"
    benchmark:
        benchmark_dir_path / "{N}.merged_ids.benchmark.txt"
    resources:
        cpus=config["common_ids_threads"],
        time=config["common_ids_threads"],
        mem=config["common_ids_threads"]
    shell:
        "workflow/scripts/merged_sequences.py "
        "--input {input.common_ids} "
        "--single_copy_files {params.single_copy_files} "
        "--outdir {output.merged_ids} 2> {log.std}"


rule mafft_dna:
    input:
        fna=output_dir_path / "merged_sequences/{N}/"
    output:
        outdir=output_dir_path / "mafft_tmp" / "{N}"
    params:
        mafft_path=config["mafft_path"]
    log:
        std=log_dir_path / "{N}.fna.mafft.log",
        cluster_log=cluster_log_dir_path / "{N}.fna.mafft.cluster.log",
        cluster_err=cluster_log_dir_path / "{N}.fna.mafft.cluster.err"
    benchmark:
        benchmark_dir_path / "{N}.fna.mafft.benchmark.txt"
    # conda:
    #     "../../%s" % config["conda_config"]
    resources:
        cpus=config["mafft_threads"],
        time=config["mafft_time"],
        mem=config["mafft_mem_mb"]
    threads:
        config["mafft_threads"]
    shell:
        "mkdir -p {output.outdir}; "
        "for FILE in `ls {input.fna}/*`; do "
        "{params.mafft_path}/mafft --thread {threads} $FILE > {output.outdir}/$(basename $FILE); "
        "done"

# rule mafft:
#     input:
#         directory(output_dir_path / "tmp" / "{sample}")
#     output:
#         outdir=temp(directory(output_dir_path / "mafft_tmp" / "{sample}"))
#     params:
#         mafft_path=config["mafft_path"]
#     log:
#         std=log_dir_path / "{sample}.mafft.log",
#         cluster_log=cluster_log_dir_path / "{sample}.mafft.cluster.log",
#         cluster_err=cluster_log_dir_path / "{sample}.mafft.cluster.err"
#     benchmark:
#         benchmark_dir_path / "{sample}.mafft.benchmark.txt"
#     # conda:
#     #     "../../%s" % config["conda_config"]
#     resources:
#         cpus=config["mafft_threads"],
#         time=config["mafft_time"],
#         mem=config["mafft_mem_mb"]
#     threads:
#         config["mafft_threads"]
#     shell:
#         "mkdir -p {output.outdir}; "
#         "for FILE in `ls -d {input}/*`; do "
#         "FILE=$(basename $FILE); "
#         "{params.mafft_path}/mafft --thread {threads} results/busco/merged_sequences/merged_$FILE.fna > {output.outdir}/merged_$FILE.fna 2> {log.std}; "
#         "{params.mafft_path}/mafft --thread {threads} --anysymbol results/busco/merged_sequences/merged_$FILE.faa > {output.outdir}/merged_$FILE.faa 2> {log.std}; "
#         "done"

rule mafft_protein:
    input:
        faa=output_dir_path / "merged_sequences/{N}/merged_{sample}.faa"
    output:
        outfile=temp(output_dir_path / "mafft_tmp" / "{N}/{sample}.faa")
    params:
        mafft_path=config["mafft_path"]
    log:
        std=log_dir_path / "{N}.{sample}.faa.mafft.log",
        cluster_log=cluster_log_dir_path / "{N}.{sample}.faa.mafft.cluster.log",
        cluster_err=cluster_log_dir_path / "{N}.{sample}.faa.mafft.cluster.err"
    benchmark:
        benchmark_dir_path / "{N}.{sample}.faa.mafft.benchmark.txt"
    # conda:
    #     "../../%s" % config["conda_config"]
    resources:
        cpus=config["mafft_threads"],
        time=config["mafft_time"],
        mem=config["mafft_mem_mb"]
    threads:
        config["mafft_threads"]
    shell:
        "{params.mafft_path}/mafft --anysymbol --thread {threads} {input.faa} > {output.outfile} 2> {log.std}"