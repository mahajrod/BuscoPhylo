import os
localrules: merged_sequences
ruleorder:merged_sequences > mafft_fna


checkpoint merged_sequences:
    input:
        common_ids=busco_dir_path / "single_copy_busco_sequences.common.ids",
        single_copy_files=expand(busco_dir_path / "{species}/single_copy_busco_sequences", species=config["species_list"])
    output:
        merged_ids=directory(busco_dir_path / "merged_sequences")
    log:
        std=log_dir_path / "merged_ids.log",
        cluster_log=cluster_log_dir_path / "merged_ids.cluster.log",
        cluster_err=cluster_log_dir_path / "merged_ids.cluster.err"
    benchmark:
        benchmark_dir_path / "merged_ids.benchmark.txt"
    resources:
        cpus=config["common_ids_threads"],
        time=config["common_ids_threads"],
        mem=config["common_ids_threads"]
    shell:
        "workflow/scripts/merged_sequences.py "
        "--input {input.common_ids} "
        "--single_copy_files {input.single_copy_files} "
        "--outdir {output.merged_ids} 2> {log.std}"

def mafft_input(wildcards):
    checkpoint_output = checkpoints.merged_sequences.get(**wildcards).output[0]
    return expand(busco_dir_path / "merged_sequences" / "merged_{i}.fna",
           i=glob_wildcards(os.path.join(checkpoint_output, '{i}.fna')).i)


rule mafft_fna:
    input:
        mafft_input
    output:
        outfile=directory(mafft_dir_path)
    params:
        outfile=mafft_dir_path / "{sample}",
        mafft_path=config["mafft_path"]
    log:
        std=log_dir_path / "{sample}.mafft_fna.log",
        cluster_log=cluster_log_dir_path / "{sample}.mafft_fna.cluster.log",
        cluster_err=cluster_log_dir_path / "{sample}.mafft_fna.cluster.err"
    benchmark:
        benchmark_dir_path / "{sample}.mafft_fna.benchmark.txt"
    conda:
        "../../%s" % config["conda_config"]
    resources:
        cpus=config["mafft_threads"],
        time=config["mafft_time"],
        mem=config["mafft_mem_mb"]
    threads:
        config["mafft_threads"]
    shell:
        "{params.mafft_path}/mafft --thread {threads} {input} > {params.outfile} 2> {log.std}"




# rule merged_sequences:
#     input:
#         common_ids=busco_dir_path / "single_copy_busco_sequences.common.ids",
#         single_copy_files=expand(busco_dir_path / "{species}/single_copy_busco_sequences", species=config["species_list"])
#     output:
#         merged_ids=directory(busco_dir_path / "merged_sequences")
#     log:
#         std=log_dir_path / "merged_ids.log",
#         cluster_log=cluster_log_dir_path / "merged_ids.cluster.log",
#         cluster_err=cluster_log_dir_path / "merged_ids.cluster.err"
#     benchmark:
#         benchmark_dir_path / "merged_ids.benchmark.txt"
#     resources:
#         cpus=config["common_ids_threads"],
#         time=config["common_ids_threads"],
#         mem=config["common_ids_threads"]
#     shell:
#         "workflow/scripts/merged_sequences.py "
#         "--input {input.common_ids} "
#         "--single_copy_files {input.single_copy_files} "
#         "--outdir {output.merged_ids} 2> {log.std}"


# rule mafft_fna:
#     input:
#         fna=busco_dir_path / "merged_sequences" / "merged_{sample}.fna"
#     output:
#         outfile=mafft_dir_path / "{sample}"
#     params:
#         mafft_path=config["mafft_path"]
#     log:
#         std=log_dir_path / "{sample}.mafft_fna.log",
#         cluster_log=cluster_log_dir_path / "{sample}.mafft_fna.cluster.log",
#         cluster_err=cluster_log_dir_path / "{sample}.mafft_fna.cluster.err"
#     benchmark:
#         benchmark_dir_path / "{sample}.mafft_fna.benchmark.txt"
#     conda:
#         "../../%s" % config["conda_config"]
#     resources:
#         cpus=config["mafft_threads"],
#         time=config["mafft_time"],
#         mem=config["mafft_mem_mb"]
#     threads:
#         config["mafft_threads"]
#     shell:
#         "{params.mafft_path}/mafft --thread {threads} {input.fna} > {output.outfile} 2> {log.std}"