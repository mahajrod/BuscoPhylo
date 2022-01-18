# localrules: mafft_dna, mafft_protein


rule mafft_dna:
    input:
        fna_list=merged_sequences_dir_path / "{N}/"
    output:
        temp(directory(mafft_dir_path / "fna_tmp" / "{N}"))
    params:
        mafft_path=config["mafft_path"],
        options=config["mafft_dna_params"]
    log:
        std=log_dir_path / "mafft_dna.{N}.log",
        cluster_log=cluster_log_dir_path / "mafft_dna.{N}.cluster.log",
        cluster_err=cluster_log_dir_path / "mafft_dna.{N}.cluster.err"
    benchmark:
        benchmark_dir_path / "mafft_dna.{N}.benchmark.txt"
    conda:
        "../../%s" % config["conda_config"]
    resources:
        cpus=config["mafft_threads"],
        time=config["mafft_time"],
        mem=config["mafft_mem_mb"]
    threads:
        config["mafft_threads"]
    shell:
        "mkdir -p {output}; "
        "for FILE in `ls {input.fna_list}/*`; do "
        "mafft --thread {threads} {params.options} ${{FILE%.*}}.fna > {output}/$(basename ${{FILE%.*}}.fna) 2> {log.std}; "
        "done; "


rule mafft_protein:
    input:
        faa_list=merged_sequences_dir_path / "{N}/"
    output:
        temp(directory(mafft_dir_path / "faa_tmp" / "{N}"))
    params:
        mafft_path=config["mafft_path"],
        options=config["mafft_protein_params"]
    log:
        std=log_dir_path / "mafft_protein.{N}.log",
        cluster_log=cluster_log_dir_path / "mafft_protein.{N}.cluster.log",
        cluster_err=cluster_log_dir_path / "mafft_protein.{N}.cluster.err"
    benchmark:
        benchmark_dir_path / "mafft_protein.{N}.benchmark.txt"
    conda:
        "../../%s" % config["conda_config"]
    resources:
        cpus=config["mafft_threads"],
        time=config["mafft_time"],
        mem=config["mafft_mem_mb"]
    threads:
        config["mafft_threads"]
    shell:
        "mkdir -p {output}; "
        "for FILE in `ls {input.faa_list}/*`; do "
        "mafft --thread {threads} {params.options} ${{FILE%.*}}.faa > {output}/$(basename ${{FILE%.*}}.faa) 2> {log.std}; "
        "done; "


