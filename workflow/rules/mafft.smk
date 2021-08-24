import os
localrules: merged_sequences
ruleorder:merged_sequences > mafft_fna

rule merged_sequences:
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

def ids_list(single_copy_busco_sequences):
    result = []
    with open(single_copy_busco_sequences, 'r') as file:
        for line in file:
            result.append(line.strip())
    return result

rule mafft_fna:
    input:
        fna=expand(busco_dir_path / "merged_sequences" / "merged_{sample}.fna", sample=ids_list(busco_dir_path / "single_copy_busco_sequences.common.ids"))
    output:
        directory(mafft_dir_path)
    params:
        outfile=mafft_dir_path / "{sample}",
        mafft_path=config["mafft_path"]
    log:
        std=log_dir_path / "{sample}.mafft_fna.log",
        cluster_log=cluster_log_dir_path / "{sample}.mafft_fna.cluster.log",
        cluster_err=cluster_log_dir_path / "{sample}.mafft_fna.cluster.err"
    benchmark:
        benchmark_dir_path / "mafft_fna.benchmark.txt"
    conda:
        "../../%s" % config["conda_config"]
    resources:
        cpus=config["mafft_threads"],
        time=config["mafft_time"],
        mem=config["mafft_mem_mb"],
    shell:
        "mkdir -p {output}; {params.mafft_path}/mafft {input.fna} > {params.outfile} 2> {log.std}"


# checkpoint mafft_tasks_list:
#     input:
#         merged_ids=directory(busco_dir_path / "merged_sequences")
#     output:
#         mafft_tasks=directory(mafft_dir_path / "slurm")
#     params:
#         amount_of_tasks = 20,
#         file_extension = "faa",
#         mafft_command_outdir = mafft_dir_path / "output"
#     log:
#         std=log_dir_path / "mafft_tasks_list.log",
#         cluster_log=cluster_log_dir_path / "mafft_tasks_list.cluster.log",
#         cluster_err=cluster_log_dir_path / "mafft_tasks_list.cluster.err"
#     benchmark:
#         benchmark_dir_path / "mafft_tasks_list.benchmark.txt"
#     resources:
#         cpus=config["common_ids_threads"],
#         time=config["common_ids_threads"],
#         mem=config["common_ids_threads"]
#     shell:
#         "workflow/scripts/mafft_tasks_list.py "
#         "--input {input.merged_ids} "
#         "--file-extension {params.file_extension} "
#         "--amount {params.amount_of_tasks} "
#         "--mafft_command_outdir {params.mafft_command_outdir} "
#         "--outdir {output.mafft_tasks} > {log.std} 2>&1"