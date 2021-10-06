# localrules: mrbayes_dna, mrbayes_protein


rule mrbayes_dna:
    input:
        concat_aligments_dir_path / "concat.aln.fna.nexus"
    output:
        directory(mrbayes_dir_path / "fna")
    params:
        mrbayes_path=config["mrbayes_path"],
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
        "mpirun -np {resources.cpus} {params.mrbayes_path}/mb {input} {params.options} 1> {log.std} 2>&1; "
        "mv {input}.* {output}; "