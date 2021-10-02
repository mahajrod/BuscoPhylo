# localrules: gblocks_dna, gblocks_protein

rule gblocks_dna:
    input:
        fna=directory(mafft_dir_path / "fna_tmp" / "{N}")
    output:
        directory(gblocks_dir_path / "fna_tmp" /"{N}")
    params:
        gblocks_path=config["gblocks_path"],
        options=config["gblocks_dna_params"]
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
        "for FILE in `ls {input.fna}/*`; do "
        "{params.gblocks_path}/Gblocks ${{FILE%.*}}.fna {params.options} 1> {log.std} 2>&1; "
        "touch ${{FILE%.*}}.fna-gb; touch ${{FILE%.*}}.fna-gb.txt; "
        "mv ${{FILE%.*}}.fna-gb {output}/; mv ${{FILE%.*}}.fna-gb.txt {output}/; "
        "done"


rule gblocks_protein:
    input:
        faa=directory(mafft_dir_path / "faa_tmp" / "{N}")
    output:
        directory(gblocks_dir_path / "faa_tmp" /"{N}")
    params:
        gblocks_path=config["gblocks_path"],
        options=config["gblocks_protein_params"]
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
        "for FILE in `ls {input.faa}/*`; do "
        "{params.gblocks_path}/Gblocks ${{FILE%.*}}.faa {params.options} 1> {log.std} 2>&1; "
        "touch ${{FILE%.*}}.faa-gb; touch ${{FILE%.*}}.faa-gb.txt; "
        "mv ${{FILE%.*}}.faa-gb {output}/; mv ${{FILE%.*}}.faa-gb.txt {output}/; "
        "done"