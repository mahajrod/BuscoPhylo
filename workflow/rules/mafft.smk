import os
localrules: merged_sequences


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


rule mafft:
    input:
        fna=busco_dir_path / "merged_sequences" / "merged_{sample}.{extension}"
    output:
        outfile=mafft_dir_path / "{sample}.{extension}"
    params:
        mafft_path=config["mafft_path"]
    log:
        std=log_dir_path / "{sample}.{extension}.mafft.log",
        cluster_log=cluster_log_dir_path / "{sample}.{extension}.mafft.cluster.log",
        cluster_err=cluster_log_dir_path / "{sample}.{extension}.mafft.cluster.err"
    benchmark:
        benchmark_dir_path / "{sample}.{extension}.mafft.benchmark.txt"
    # conda:
    #     "../../%s" % config["conda_config"]
    resources:
        cpus=config["mafft_threads"],
        time=config["mafft_time"],
        mem=config["mafft_mem_mb"]
    threads:
        config["mafft_threads"]
    shell:
        "{params.mafft_path}/mafft --thread {threads} {input.fna} > {output.outfile} 2> {log.std}"
