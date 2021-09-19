localrules: mafft_dna, mafft_protein


rule mafft_dna:
    input:
        fna=merged_sequences_dir_path / "{N}/"
    output:
        outdir=temp(directory(mafft_dir_path / "fna_tmp" / "{N}"))
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


rule mafft_protein:
    input:
        faa=merged_sequences_dir_path / "{N}/"
    output:
        outdir=temp(directory(mafft_dir_path / "faa_tmp" / "{N}"))
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


