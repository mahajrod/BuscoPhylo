localrules: merged_sequences, tmp
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

def mafft_dna_input(wildcards):
    checkpoint_output = checkpoints.merged_sequences.get(**wildcards).output[0]
    file_names = expand(mafft_dir_path / "{sample}.fna",
                        sample = glob_wildcards(os.path.join(checkpoint_output, "merged_{sample}.fna")).sample)
    return file_names

def mafft_protein_input(wildcards):
    checkpoint_output = checkpoints.merged_sequences.get(**wildcards).output[0]
    file_names = expand(mafft_dir_path / "{sample}.faa",
                        sample = glob_wildcards(os.path.join(checkpoint_output, "merged_{sample}.faa")).sample)
    return file_names

rule tmp:
    input:
        fna = busco_dir_path / "merged_sequences" / "merged_{sample}.fna",
        faa=busco_dir_path / "merged_sequences" / "merged_{sample}.faa"
    output:
        fna="merged_sequences_tmp/{sample}.fna",
        faa="merged_sequences_tmp/{sample}.faa"
    shell:
        "mv {input.fna} {output.fna}; "
        "mv {input.faa} {output.faa} "

rule mafft_dna:
    input:
        rules.tmp.output.fna
    output:
        outfile=mafft_dir_path / "{sample}.fna"
    params:
        fna = expand("merged_sequences_tmp/merged_{sample}", sample = rules.tmp.output.fna),
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
        rules.tmp.output.faa
    output:
        outfile=mafft_dir_path / "{sample}.faa"
    params:
        faa = expand("merged_sequences_tmp/merged_{sample}", sample = rules.tmp.output.faa),
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

rule finished:
    input:
        mafft_dna_input,
        mafft_protein_input
    output:
        "finished.txt"
    shell:
        "touch {output}"