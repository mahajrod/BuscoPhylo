localrules: gblocks_dna, gblocks_protein

rule gblocks_dna:
    input:
        fna=directory(mafft_dir_path / "fna" / "{N}")
    output:
        directory(gblocks_dir_path / "fna" /"{N}")
    params:
        gblocks_path=config["gblocks_path"],
        gblocks_dna_flags="-t=d -p=t"
    log:
        std=log_dir_path / "{N}.fna.gblocks.log",
        cluster_log=cluster_log_dir_path / "{N}.fna.gblocks.cluster.log",
        cluster_err=cluster_log_dir_path / "{N}.fna.gblocks.cluster.err"
    benchmark:
        benchmark_dir_path / "{N}.fna.gblocks.benchmark.txt"
    # conda:
    #     "../../%s" % config["conda_config"]
    resources:
        cpus=config["gblocks_threads"],
        time=config["gblocks_time"],
        mem=config["gblocks_mem_mb"]
    shell:
        "mkdir -p {output}; "
        "sleep 10; "
        "for FILE in `ls {input.fna}/*`; do "
        "{params.gblocks_path}/Gblocks ${{FILE%.*}}.fna {params.gblocks_dna_flags}; "
        "mv ${{FILE%.*}}.fna-gb {output}/; mv ${{FILE%.*}}.fna-gb.txt {output}/; "
        "done"


rule gblocks_protein:
    input:
        faa=directory(mafft_dir_path / "faa" / "{N}")
    output:
        directory(gblocks_dir_path / "faa" /"{N}")
    params:
        gblocks_path=config["gblocks_path"],
        gblocks_protein_flags="-t=p -p=t"
    log:
        std=log_dir_path / "{N}.faa.gblocks.log",
        cluster_log=cluster_log_dir_path / "{N}.faa.gblocks.cluster.log",
        cluster_err=cluster_log_dir_path / "{N}.faa.gblocks.cluster.err"
    benchmark:
        benchmark_dir_path / "{N}.faa.gblocks.benchmark.txt"
    # conda:
    #     "../../%s" % config["conda_config"]
    resources:
        cpus=config["gblocks_threads"],
        time=config["gblocks_time"],
        mem=config["gblocks_mem_mb"]
    shell:
        "mkdir -p {output}; "
        "sleep 10; "
        "for FILE in `ls {input.faa}/*`; do "
        "{params.gblocks_path}/Gblocks ${{FILE%.*}}.faa {params.gblocks_protein_flags}; "
        "mv ${{FILE%.*}}.faa-gb {output}/; mv ${{FILE%.*}}.faa-gb.txt {output}/; "
        "done"