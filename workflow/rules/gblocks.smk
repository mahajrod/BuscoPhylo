rule gblocks_dna:
    input:
        fna=mafft_dir_path / "{sample}.fna"
    output:
        gb=gblocks_dir_path / "{sample}.fna-gb",
        gb_txt=gblocks_dir_path / "{sample}.fna-gb.txt"
    params:
        gblocks_path=config["gblocks_path"],
        gb_tmp=mafft_dir_path / "{sample}.fna-gb",
        gb_txt_tmp=mafft_dir_path / "{sample}.fna-gb.txt"
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
        "mv {params.gb_tmp} {output.gb};"
        "mv {params.gb_txt_tmp} {output.gb_txt}"