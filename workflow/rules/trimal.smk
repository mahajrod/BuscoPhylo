localrules: trimal_dna, trimal_protein

rule trimal_dna:
    input:
        fna=directory(mafft_dir_path / "fna_tmp" / "{N}")
    output:
        temp(directory(trimal_dir_path / "fna_tmp" /"{N}"))
    params:
        trimal_path=config["trimal_path"],
        trimal_dna_flags="-gt 0.9 -cons 60 "
    log:
        std=log_dir_path / "{N}.fna.trimal.log",
        cluster_log=cluster_log_dir_path / "{N}.fna.trimal.cluster.log",
        cluster_err=cluster_log_dir_path / "{N}.fna.trimal.cluster.err"
    benchmark:
        benchmark_dir_path / "{N}.fna.trimal.benchmark.txt"
    # conda:
    #     "../../%s" % config["conda_config"]
    resources:
        cpus=config["trimal_threads"],
        time=config["trimal_time"],
        mem=config["trimal_mem_mb"]
    shell:
        "mkdir -p {output}; "
        "for FILE in `ls {input.fna}/*`; do "
        "{params.trimal_path}/trimal -in ${{FILE%.*}}.fna -out {output}/$(basename ${{FILE%.*}}.fna) {params.trimal_dna_flags} 1> {log.std} 2> {log.std}; "
        "done"


rule trimal_protein:
    input:
        faa=directory(mafft_dir_path / "faa_tmp" / "{N}")
    output:
        temp(directory(trimal_dir_path / "faa_tmp" /"{N}"))
    params:
        trimal_path=config["trimal_path"],
        trimal_protein_flags="-gt 0.9 -cons 60 "
    log:
        std=log_dir_path / "{N}.faa.trimal.log",
        cluster_log=cluster_log_dir_path / "{N}.faa.trimal.cluster.log",
        cluster_err=cluster_log_dir_path / "{N}.faa.trimal.cluster.err"
    benchmark:
        benchmark_dir_path / "{N}.faa.trimal.benchmark.txt"
    # conda:
    #     "../../%s" % config["conda_config"]
    resources:
        cpus=config["trimal_threads"],
        time=config["trimal_time"],
        mem=config["trimal_mem_mb"]
    shell:
        "mkdir -p {output}; "
        "for FILE in `ls {input.faa}/*`; do "
        "{params.trimal_path}/trimal -in ${{FILE%.*}}.faa -out {output}/$(basename ${{FILE%.*}}.faa) {params.trimal_protein_flags} 1> {log.std} 2> {log.std}; "
        "done"