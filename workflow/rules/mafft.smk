localrules: merged_sequences, combine_samples


checkpoint merged_sequences:
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

rule combine_samples:
    input:
        fna=lambda w: expand_template_from_merged_sequences(w, merged_sequences_dir_path / "merged_{sample}.fna"),
        faa=lambda w: expand_template_from_merged_sequences(w, merged_sequences_dir_path / "merged_{sample}.faa"),
    output:
        fna=temp("fna.txt"),
        faa=temp("faa.txt"),
    shell:
        "cat {input.fna} > {output.fna}"
        "cat {input.faa} > {output.faa}"

rule mafft_dna:
    input:
        f=rules.combine_samples.output.fna,
        fna=merged_sequences_dir_path / "merged_{sample}.fna"
    output:
        outfile=mafft_dir_path / "{sample}.fna"
    params:
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
        "{params.mafft_path}/mafft --thread {threads} {input.fna} > {output.outfile} 2> {log.std}"


rule mafft_protein:
    input:
        f = rules.combine_samples.output.faa,
        faa=merged_sequences_dir_path / "merged_{sample}.faa"
    output:
        outfile=mafft_dir_path / "{sample}.faa"
    params:
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
        "{params.mafft_path}/mafft --anysymbol --thread {threads} {input.faa} > {output.outfile} 2> {log.std}"


