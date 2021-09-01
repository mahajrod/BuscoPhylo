# localrules: merged_sequences
ruleorder: merged_sequences > tmp > mafft_dna > mafft_protein

checkpoint merged_sequences:
    input:
        common_ids=busco_dir_path / "single_copy_busco_sequences.common.ids"
    output:
        merged_ids=directory(busco_dir_path / "merged_sequences")
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


rule tmp:
    input:
        fna = busco_dir_path / "merged_sequences" / "merged_{sample}.fna",
        faa=busco_dir_path / "merged_sequences" / "merged_{sample}.faa"
    output:
        directory("merged_sequences_tmp")
    shell:
        "mv {input.faa} {output}; "
        "mv {input.fna} {output} "

rule mafft_dna:
    input:
        rules.tmp.output
    output:
        outfile=mafft_dir_path / "{sample}.fna"
    params:
        fna = expand("merged_sequences_tmp/merged_{sample}.fna", sample = rules.tmp.output),
        mafft_path=config["mafft_path"]
    log:
        std=log_dir_path / "{sample}.fna.mafft.log",
        cluster_log=cluster_log_dir_path / "{sample}.fna.mafft.cluster.log",
        cluster_err=cluster_log_dir_path / "{sample}.fna.mafft.cluster.err"
    benchmark:
        benchmark_dir_path / "{sample}.fna.mafft.benchmark.txt"
    # conda:
    #     "../../%s" % config["conda_config"]
    resources:
        cpus=config["mafft_threads"],
        time=config["mafft_time"],
        mem=config["mafft_mem_mb"]
    threads:
        config["mafft_threads"]
    shell:
        "{params.mafft_path}/mafft --thread {threads} {params.fna} > {output.outfile} 2> {log.std}"


rule mafft_protein:
    input:
        rules.tmp.output
    output:
        outfile=mafft_dir_path / "{sample}.faa"
    params:
        faa = expand("merged_sequences_tmp/merged_{sample}.faa", sample = [os.path.splitext(filename)[0] for filename in os.listdir(rules.tmp.output)]),
        mafft_path=config["mafft_path"]
    log:
        std=log_dir_path / "{sample}.faa.mafft.log",
        cluster_log=cluster_log_dir_path / "{sample}.faa.mafft.cluster.log",
        cluster_err=cluster_log_dir_path / "{sample}.faa.mafft.cluster.err"
    benchmark:
        benchmark_dir_path / "{sample}.faa.mafft.benchmark.txt"
    # conda:
    #     "../../%s" % config["conda_config"]
    resources:
        cpus=config["mafft_threads"],
        time=config["mafft_time"],
        mem=config["mafft_mem_mb"]
    threads:
        config["mafft_threads"]
    shell:
        "{params.mafft_path}/mafft --anysymbol --thread {threads} {params.faa} > {output.outfile} 2> {log.std}"