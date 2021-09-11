# localrules: merged_sequences


checkpoint merged_sequences: #TODO: добавить в скрипт разделение по дирректориям
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
    group:
        "mafft_group0"
    resources:
        cpus=config["common_ids_threads"],
        time=config["common_ids_threads"],
        mem=config["common_ids_threads"]
    shell:
        "workflow/scripts/merged_sequences.py "
        "--input {input.common_ids} "
        "--single_copy_files {params.single_copy_files} "
        "--outdir {output.merged_ids} 2> {log.std}"


rule mafft_dna: #TODO: sample - это дирректория, в которой опред. кол-во файлов. Тут просто цикл по этим файлам
    input:
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
    group:
        "mafft_group1"
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
    group:
        "mafft_group1"
    resources:
        cpus=config["mafft_threads"],
        time=config["mafft_time"],
        mem=config["mafft_mem_mb"]
    threads:
        config["mafft_threads"]
    shell:
        "{params.mafft_path}/mafft --anysymbol --thread {threads} {input.faa} > {output.outfile} 2> {log.std}"

def expand_template_from_merged_sequences(wildcards, template):
    checkpoint_output = checkpoints.merged_sequences.get(**wildcards).output[0]
    sample, = glob_wildcards(os.path.join(checkpoint_output, "merged_{sample}.fna"))
    return expand(str(template), sample=sample)

rule crutch:
    input:
        lambda w: expand_template_from_merged_sequences(w, mafft_dir_path / "{sample}.fna"),
        lambda w: expand_template_from_merged_sequences(w, mafft_dir_path / "{sample}.faa"),
    output:
        "tmp.txt"
    group:
        "mafft_group1"
    shell:
        "touch {output}"