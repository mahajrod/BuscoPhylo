rule gblocks_dna:
    input:
        fna=mafft_dir_path / "{sample}.fna"
    output:
        gb=gblocks_dir_path / "{sample}.fna-gb",
        gb_txt=gblocks_dir_path / "{sample}.fna-gb.txt"
    params:
        gblocks_path=config["gblocks_path"]
    log:
        std=log_dir_path / "{sample}.fna.gblocks.log",
        cluster_log=cluster_log_dir_path / "{sample}.fna.gblocks.cluster.log",
        cluster_err=cluster_log_dir_path / "{sample}.fna.gblocks.cluster.err"
    benchmark:
        benchmark_dir_path / "{sample}.fna.gblocks.benchmark.txt"
    # conda:
    #     "../../%s" % config["conda_config"]
    resources:
        cpus=config["gblocks_threads"],
        time=config["gblocks_time"],
        mem=config["gblocks_mem_mb"]
    # threads:
    #     config["gblocks_threads"]
    shell:
        "{params.gblocks_path}/Gblocks {input.fna} -t=d -p=t 1> {log.std}; "
        "mv {input.fna}-gb {output.gb};"
        "mv {input.fna}-gb.txt {output.gb_txt}"