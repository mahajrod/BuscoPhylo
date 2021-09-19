localrules: common_ids, species_ids, merged_sequences


rule species_ids:
    input:
        busco_outdir=directory(busco_dir_path / "{species}")
    output:
        ids=common_ids_dir_path / "{species}.ids"
    params:
        single_copy_files="single_copy_busco_sequences/"
    log:
        std=log_dir_path / "{species}.species_ids.log",
        cluster_log=cluster_log_dir_path / "{species}.species_ids.cluster.log",
        cluster_err=cluster_log_dir_path / "{species}.species_ids.cluster.err"
    benchmark:
        benchmark_dir_path / "{species}/species_ids.benchmark.txt"
    resources:
        cpus=config["common_ids_threads"],
        time=config["common_ids_time"],
        mem=config["common_ids_mem_mb"]
    shell:
        "ls {input.busco_outdir}/{params.single_copy_files} | grep -P '.fna$' | sed 's/.fna//' > {output.ids} 2> {log.std}"


checkpoint common_ids:
    input:
        ids=expand(common_ids_dir_path / "{species}.ids", species=config["species_list"])
    output:
        directory(single_copy_busco_sequences_dir_path)
    params:
        nfiles=len(config["species_list"]),
        prefix="common.ids",
        split_size=50
    log:
        std=log_dir_path / "common_ids.log",
        cluster_log=cluster_log_dir_path / "common_ids.cluster.log",
        cluster_err=cluster_log_dir_path / "common_ids.cluster.err"
    benchmark:
        benchmark_dir_path / "common_ids.benchmark.txt"
    resources:
        cpus=config["common_ids_threads"],
        time=config["common_ids_time"],
        mem=config["common_ids_mem_mb"]
    shell:
        "mkdir -p {output}; "
        "cat {input.ids} | "
        "sort | uniq -c | awk '{{if($1=={params.nfiles}){{print $2}}}}' | "
        "split -l {params.split_size} --numeric-suffixes - {output}/{params.prefix}"


checkpoint merged_sequences:
    input:
        common_ids=single_copy_busco_sequences_dir_path / "common.ids{N}"
    output:
        merged_ids=directory(merged_sequences_dir_path / "{N}/")
    params:
        single_copy_files=expand(busco_dir_path / "{species}" / "single_copy_busco_sequences", species=config["species_list"])
    log:
        std=log_dir_path / "{N}.merged_ids.log",
        cluster_log=cluster_log_dir_path / "{N}.merged_ids.cluster.log",
        cluster_err=cluster_log_dir_path / "{N}.merged_ids.cluster.err"
    benchmark:
        benchmark_dir_path / "{N}.merged_ids.benchmark.txt"
    resources:
        cpus=config["common_ids_threads"],
        time=config["common_ids_threads"],
        mem=config["common_ids_threads"]
    shell:
        "workflow/scripts/merged_sequences.py "
        "--input {input.common_ids} "
        "--single_copy_files {params.single_copy_files} "
        "--outdir {output.merged_ids} 2> {log.std}"