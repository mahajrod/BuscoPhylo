# localrules: iqtree_dna, iqtree_protein

rule iqtree_dna:
    input:
        concat_aligments_dir_path / concat_fna_filename
    output:
        directory(iqtree_dir_path / "fna")
    params:
        iqtree_path=config["iqtree_path"],
        prefix=config["iqtree_fna_prefix"],
        options=config["iqtree_dna_params"]
    log:
        std=log_dir_path / "fna.iqtree.log",
        cluster_log=cluster_log_dir_path / "fna.iqtree.cluster.log",
        cluster_err=cluster_log_dir_path / "fna.iqtree.cluster.err"
    benchmark:
        benchmark_dir_path / "fna.iqtree.benchmark.txt"
    # conda:
    #     "../../%s" % config["conda_config"]
    resources:
        cpus=config["iqtree_threads"],
        time=config["iqtree_time"],
        mem=config["iqtree_mem_mb"]
    shell:
        "mkdir -p {output}; "
        "DIR=$(dirname {input}); "
        "{params.iqtree_path}/iqtree -s {input} -pre {params.prefix} -nt {resources.cpus} {params.options} 1> {log.std} 2>&1; "
        "mv $DIR/{params.prefix}.bionj {output}; "
        "mv $DIR/{params.prefix}.ckp.gz {output}; "
        "mv $DIR/{params.prefix}.log {output}; "
        "mv $DIR/{params.prefix}.mldist {output}; "
        "mv $DIR/{params.prefix}.model.gz {output}; "
        "mv $DIR/{params.prefix}.treefile {output}; "
        "mv $DIR/{params.prefix}.iqtree {output}; "


rule iqtree_protein:
    input:
        concat_aligments_dir_path / concat_faa_filename
    output:
        directory(iqtree_dir_path / "faa")
    params:
        iqtree_path=config["iqtree_path"],
        prefix=config["iqtree_faa_prefix"],
        options=config["iqtree_protein_params"]
    log:
        std=log_dir_path / "faa.iqtree.log",
        cluster_log=cluster_log_dir_path / "faa.iqtree.cluster.log",
        cluster_err=cluster_log_dir_path / "faa.iqtree.cluster.err"
    benchmark:
        benchmark_dir_path / "faa.iqtree.benchmark.txt"
    # conda:
    #     "../../%s" % config["conda_config"]
    resources:
        cpus=config["iqtree_threads"],
        time=config["iqtree_time"],
        mem=config["iqtree_mem_mb"]
    shell:
        "mkdir -p {output}; "
        "DIR=$(dirname {input}); "
        "{params.iqtree_path}/iqtree -s {input} -pre {params.prefix} -nt {resources.cpus} {params.options} 1> {log.std} 2>&1; "
        "mv $DIR/{params.prefix}.bionj {output}; "
        "mv $DIR/{params.prefix}.ckp.gz {output}; "
        "mv $DIR/{params.prefix}.log {output}; "
        "mv $DIR/{params.prefix}.mldist {output}; "
        "mv $DIR/{params.prefix}.model.gz {output}; "
        "mv $DIR/{params.prefix}.treefile {output}; "
        "mv $DIR/{params.prefix}.iqtree {output}; "
