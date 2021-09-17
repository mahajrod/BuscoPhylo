localrules: trimal_dna, trimal_protein

rule trimal_dna:
    input:
        fna=directory(mafft_dir_path / "fna" / "{N}")
    output:
        directory(trimal_dir_path / "fna" /"{N}")
    params:
        trimal_path=config["trimal_path"],
        trimal_dna_flags="-t=d -p=t"
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
        "{params.gblocks_path}/Gblocks ${{FILE%.*}}.fna {params.gblocks_dna_flags} 1> {log.std} 2> {log.std}; "
        "sleep 10; "
        "mv ${{FILE%.*}}.fna-gb {output}/; mv ${{FILE%.*}}.fna-gb.txt {output}/; "
        "done"


rule trimal_protein:
    input:
        faa=directory(mafft_dir_path / "faa" / "{N}")
    output:
        directory(trimal_dir_path / "faa" /"{N}")
    params:
        trimal_path=config["trimal_path"],
        trimal_protein_flags="-t=p -p=t"
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
        "{params.gblocks_path}/Gblocks ${{FILE%.*}}.faa {params.gblocks_protein_flags} 1> {log.std} 2> {log.std}; "
        "sleep 10; "
        "mv ${{FILE%.*}}.faa-gb {output}/; mv ${{FILE%.*}}.faa-gb.txt {output}/; "
        "done"