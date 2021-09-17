localrules: merged_sequences, mafft_dna, mafft_protein, mafft_finish


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
        outdir=directory(mafft_dir_path / "fna" / "{N}")
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
        "{params.mafft_path}/mafft --thread {threads} ${{FILE%.*}}.fna > {output.outdir}/$(basename ${{FILE%.*}}.fna) 2> {log.std}; "
        "done; "
        "sleep 10; "


rule mafft_protein:
    input:
        faa=output_dir_path / "merged_sequences/{N}/"
    output:
        outdir=directory(mafft_dir_path / "faa" / "{N}")
    params:
        mafft_path=config["mafft_path"]
    log:
        std=log_dir_path / "{N}.faa.mafft.log",
        cluster_log=cluster_log_dir_path / "{N}.faa.mafft.cluster.log",
        cluster_err=cluster_log_dir_path / "{N}.faa.mafft.cluster.err"
    benchmark:
        benchmark_dir_path / "{N}.faa.mafft.benchmark.txt"
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
        "for FILE in `ls {input.faa}/*`; do "
        "{params.mafft_path}/mafft --thread {threads} ${{FILE%.*}}.faa > {output.outdir}/$(basename ${{FILE%.*}}.faa) 2> {log.std}; "
        "done; "
        "sleep 10; "


