localrules: gblocks

rule gblocks_dna:
    input:
        fna=directory(mafft_dir_path / "fna" / "{N}")
    output:
        directory(gblocks_dir_path / "fna" /"{N}")
    params:
        gblocks_path=config["gblocks_path"],
        gblocks_dna_flags="-t=d -p=t",
        gblocks_protein_flags="-t=p -p=t"
    log:
        std=log_dir_path / "{N}.gblocks.log",
        cluster_log=cluster_log_dir_path / "{N}.gblocks.cluster.log",
        cluster_err=cluster_log_dir_path / "{N}.gblocks.cluster.err"
    benchmark:
        benchmark_dir_path / "{N}.gblocks.benchmark.txt"
    # conda:
    #     "../../%s" % config["conda_config"]
    resources:
        cpus=config["gblocks_threads"],
        time=config["gblocks_time"],
        mem=config["gblocks_mem_mb"]
    shell:
        "mkdir -p {output}; "
        "for FILE in `ls {input.fna}/*`; do "
        "{params.gblocks_path}/Gblocks ${{FILE%.*}}.fna {params.gblocks_dna_flags} || true; "
        "mv ${{FILE%.*}}.fna-gb {output}/; mv ${{FILE%.*}}.fna-gb.txt {output}/; "
        "done"
#
# rule gblocks_dna:
#     input:
#         fna=mafft_dir_path / "{sample}.fna"
#     output:
#         gb=gblocks_dir_path / "{sample}.fna-gb",
#         gb_txt=gblocks_dir_path / "{sample}.fna-gb.txt"
#     params:
#         gblocks_dir = directory(gblocks_dir_path),
#         gblocks_path=config["gblocks_path"],
#         gblocks_flags="-t=d -p=t"
#     log:
#         std=log_dir_path / "{sample}.fna.gblocks.log",
#         cluster_log=cluster_log_dir_path / "{sample}.fna.gblocks.cluster.log",
#         cluster_err=cluster_log_dir_path / "{sample}.fna.gblocks.cluster.err"
#     benchmark:
#         benchmark_dir_path / "{sample}.fna.gblocks.benchmark.txt"
#     # conda:
#     #     "../../%s" % config["conda_config"]
#     resources:
#         cpus=config["gblocks_threads"],
#         time=config["gblocks_time"],
#         mem=config["gblocks_mem_mb"]
#     shell:
#         "mkdir -p {params.gblocks_dir}; "
#         "{params.gblocks_path}/Gblocks {input.fna} {params.gblocks_flags} 1> {log.std}; "
#         "mv {input.fna}-gb {output.gb}; "
#         "mv {input.fna}-gb.txt {output.gb_txt}"


rule gblocks_protein:
    input:
        faa=mafft_dir_path / "{sample}.faa"
    output:
        gb=gblocks_dir_path / "{sample}.faa-gb",
        gb_txt=gblocks_dir_path / "{sample}.faa-gb.txt"
    params:
        gblocks_dir = directory(gblocks_dir_path),
        gblocks_path=config["gblocks_path"],
        gblocks_flags="-t=p -p=t"
    log:
        std=log_dir_path / "{sample}.faa.gblocks.log",
        cluster_log=cluster_log_dir_path / "{sample}.faa.gblocks.cluster.log",
        cluster_err=cluster_log_dir_path / "{sample}.faa.gblocks.cluster.err"
    benchmark:
        benchmark_dir_path / "{sample}.faa.gblocks.benchmark.txt"
    # conda:
    #     "../../%s" % config["conda_config"]
    resources:
        cpus=config["gblocks_threads"],
        time=config["gblocks_time"],
        mem=config["gblocks_mem_mb"]
    shell:
        "mkdir -p {params.gblocks_dir}; "
        "{params.gblocks_path}/Gblocks {input.faa} {params.gblocks_flags} 1> {log.std}; "
        "mv {input.faa}-gb {output.gb}; "
        "mv {input.faa}-gb.txt {output.gb_txt}"