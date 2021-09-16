localrules: gblocks, gblocks_results_to_one_directory

rule gblocks:
    input:
        # fna=mafft_dir_path / "{sample}.fna"
        directory(output_dir_path / "tmp" / "{sample}")
    output:
        # gb=gblocks_dir_path / "{sample}.fna-gb",
        # gb_txt=gblocks_dir_path / "{sample}.fna-gb.txt"
        temp(directory(output_dir_path / "gblocks_tmp" / "{sample}"))
    params:
        mafft_dir=directory(mafft_dir_path),
        gblocks_path=config["gblocks_path"],
        gblocks_dna_flags="-t=d -p=t",
        gblocks_protein_flags="-t=p -p=t"
    log:
        std=log_dir_path / "{sample}.gblocks.log",
        cluster_log=cluster_log_dir_path / "{sample}.gblocks.cluster.log",
        cluster_err=cluster_log_dir_path / "{sample}.gblocks.cluster.err"
    benchmark:
        benchmark_dir_path / "{sample}.gblocks.benchmark.txt"
    # conda:
    #     "../../%s" % config["conda_config"]
    resources:
        cpus=config["gblocks_threads"],
        time=config["gblocks_time"],
        mem=config["gblocks_mem_mb"]
    shell:
        "mkdir -p {output}; "
        "for FILE in `ls -d {input}/*`; do "
        "FILE=$(basename $FILE); "
        "{params.gblocks_path}/Gblocks {params.mafft_dir}/merged_$FILE.fna {params.gblocks_dna_flags} 1> {log.std} || true; "
        "{params.gblocks_path}/Gblocks {params.mafft_dir}/merged_$FILE.faa {params.gblocks_protein_flags} 1> {log.std} || true; "
        "mv {params.mafft_dir}/merged_$FILE.fna-gb {output}/; mv {params.mafft_dir}/merged_$FILE.fna-gb.txt {output}/; "
        "mv {params.mafft_dir}/merged_$FILE.faa-gb {output}/; mv {params.mafft_dir}/merged_$FILE.faa-gb.txt {output}/; "
        "done"


def expand_template_from_directories_with_sample_names(wildcards, template):
    checkpoint_output = checkpoints.directories_with_sample_names.get(**wildcards).output[0]
    sample, = glob_wildcards(os.path.join(checkpoint_output, "{sample}"))
    sample = list(set([i.split('/')[0] for i in sample]))
    return expand(str(template), sample=sample)

rule gblocks_results_to_one_directory:
    input:
        lambda w: expand_template_from_directories_with_sample_names(w, output_dir_path / "gblocks_tmp" / "{sample}")
    output:
        directory(gblocks_dir_path)
    shell:
        "mkdir -p {output}; "
        "for i in {input}/*; do "
        "mv $i {output}/; "
        "done; "