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
        std=log_dir_path / "fna.mrbayes_dna.log",
        cluster_log=cluster_log_dir_path / "fna.mrbayes_dna.cluster.log",
        cluster_err=cluster_log_dir_path / "fna.mrbayes_dna.cluster.err"
    benchmark:
        benchmark_dir_path / "fna.mrbayes_dna.benchmark.txt"
    # conda:
    #     "../../%s" % config["conda_config"]
    resources:
        cpus=config["mrbayes_threads"],
        time=config["mrbayes_time"],
        mem=config["mrbayes_mem_mb"]
    shell:
        "mkdir -p {output}; "
        "{params.mrbayes_path}/mb {input} {params.options} 1> {log.std} 2>&1; " #mpirun -np {resources.cpus} 
        "mv {input}.* {output}/; "