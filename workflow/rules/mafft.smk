localrules: merged_sequences, directories_with_sample_names, mafft_results_to_one_directory, mafft

rule merged_sequences:
    input:
        common_ids=busco_dir_path / "single_copy_busco_sequences.common.ids"
    output:
        merged_ids=directory(merged_sequences_dir_path)
    params:
        single_copy_files=expand(busco_dir_path / "{species}" / "single_copy_busco_sequences", species=config["species_list"])
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
        "--single_copy_files {params.single_copy_files} "
        "--outdir {output.merged_ids} 2> {log.std}"


checkpoint directories_with_sample_names:
    input:
        rules.merged_sequences.output.merged_ids
    output:
        temp(ancient(directory(output_dir_path / "tmp")))
    params:
        number_of_files=20
    log:
        std=log_dir_path / "directories_with_sample_names.log",
    shell:
        "workflow/scripts/get_directories_with_sample_names.sh "
        "-i {input} "
        "-o {output} "
        "-f {params.number_of_files} 1> {log.std} 2> {log.std}"


rule mafft:
    input:
        ancient(directory(output_dir_path / "tmp" / "{sample}"))
    output:
        outdir=temp(directory(output_dir_path / "mafft_tmp" / "{sample}"))
    params:
        mafft_path=config["mafft_path"]
    log:
        std=log_dir_path / "{sample}.mafft.log",
        cluster_log=cluster_log_dir_path / "{sample}.mafft.cluster.log",
        cluster_err=cluster_log_dir_path / "{sample}.mafft.cluster.err"
    benchmark:
        benchmark_dir_path / "{sample}.mafft.benchmark.txt"
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
        "for FILE in `ls -d {input}/*`; do "
        "FILE=$(basename $FILE); "
        "{params.mafft_path}/mafft --thread {threads} results/busco/merged_sequences/merged_$FILE.fna > {output.outdir}/merged_$FILE.fna 2> {log.std}; "
        "{params.mafft_path}/mafft --thread {threads} --anysymbol results/busco/merged_sequences/merged_$FILE.faa > {output.outdir}/merged_$FILE.faa 2> {log.std}; "
        "done"