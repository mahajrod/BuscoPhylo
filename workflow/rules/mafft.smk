localrules: merged_sequences, directories_with_sample_names, mafft_one_directory

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
        temp(directory(output_dir_path / "tmp"))
    params:
        number_of_files=20
    log:
        std=log_dir_path / "directories_with_sample_names.log",
    shell:
        "workflow/scripts/get_directories_with_sample_names.sh "
        "-i {input} "
        "-o {output} "
        "-f {params.number_of_files}"


def expand_template_from_directories_with_sample_names(wildcards, template):
    checkpoint_output = checkpoints.directories_with_sample_names.get(**wildcards).output[0]
    sample, = glob_wildcards(os.path.join(checkpoint_output, "{sample}/{file}"))
    return expand(str(template), sample=sample, file="*")


rule mafft:
    input:
        # fna=merged_sequences_dir_path / "merged_{sample}.fna"
        output_dir_path / "tmp" / "{sample}"
    output:
        # outfile=mafft_dir_path / "{sample}.fna"
        outdir=mafft_dir_path / "{sample}"
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
        "for FILE in `ls {input}`; do"
        "{params.mafft_path}/mafft --thread {threads} $FILE.fna > {output.outdir}/$FILE.fna 2> {log.std} ;"
        "{params.mafft_path}/mafft --thread {threads} $FILE.faa > {output.outdir}/$FILE.faa 2> {log.std} ;"
        "done"

checkpoint mafft_one_directory:
    input:
        lambda w: expand_template_from_directories_with_sample_names(w,mafft_dir_path / "{sample}")
    output:
        directory(mafft_dir_path)
    shell:
        "for i in `ls {input}`; do"
        "for j in `ls`; do "
        "mv $i/$j/* {output}; "
        "done; "
        "done"


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
    resources:
        cpus=config["mafft_threads"],
        time=config["mafft_time"],
        mem=config["mafft_mem_mb"]
    threads:
        config["mafft_threads"]
    shell:
        "{params.mafft_path}/mafft --anysymbol --thread {threads} {input.faa} > {output.outfile} 2> {log.std}"
