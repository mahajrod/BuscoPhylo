# localrules: mrbayes_dna, mrbayes_protein


rule mrbayes_dna:
    input:
        concat_aligments_dir_path / config["concat_fna_filename"]
    output:
        directory(mrbayes_dir_path / "fna")
    params:
        mrbayes_path=config["mrbayes_path"],
        prefix=config["mrbayes_fna_prefix"],
        options=config["mrbayes_dna_params"]
    log:
        std=log_dir_path / "fna.mrbayes.log",
        cluster_log=cluster_log_dir_path / "fna.mrbayes.cluster.log",
        cluster_err=cluster_log_dir_path / "fna.mrbayes.cluster.err"
    benchmark:
        benchmark_dir_path / "fna.mrbayes.benchmark.txt"
    # conda:
    #     "../../%s" % config["conda_config"]
    resources:
        cpus=config["mrbayes_threads"],
        time=config["mrbayes_time"],
        mem=config["mrbayes_mem_mb"]
    shell:
        "mkdir -p {output}; "
        "{params.mrbayes_path}/mrbayes -s {input} -pre {params.prefix} -nt {resources.cpus} {params.options} 1> {log.std} 2>&1; "
        "mv {params.prefix}.bionj {output}; "
        "mv {params.prefix}.ckp.gz {output}; "
        "mv {params.prefix}.log {output}; "
        "mv {params.prefix}.mldist {output}; "
        "mv {params.prefix}.model.gz {output}; "
        "mv {params.prefix}.treefile {output}; "
        "mv {params.prefix}.mrbayes {output}; "