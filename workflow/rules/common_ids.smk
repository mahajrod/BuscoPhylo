localrules: common_ids, species_ids, merged_sequences


rule species_ids: # get files with IDs for each species
    input:
        directory(busco_dir_path / "{species}/single_copy_busco_sequences")
    output:
        temp(common_ids_dir_path / "{species}.ids")
    log:
        std=log_dir_path / "species_ids.{species}.log",
        cluster_log=cluster_log_dir_path / "species_ids.{species}.cluster.log",
        cluster_err=cluster_log_dir_path / "species_ids.{species}.cluster.err"
    benchmark:
        benchmark_dir_path / "species_ids.{species}.benchmark.txt"
    resources:
        cpus=config["common_ids_threads"],
        time=config["common_ids_time"],
        mem=config["common_ids_mem_mb"]
    shell:
        "ls {input} | grep -P '.fna$' | sed 's/.fna//' > {output} 2> {log.std}"


checkpoint common_ids: # get common IDs for all species and split them into files
    input:
        expand(common_ids_dir_path / "{species}.ids", species=config["species_list"])
    output:
        directory(single_copy_busco_sequences_dir_path)
    params:
        nfiles=len(config["species_list"]),
        prefix="common.ids",
        split_size=config["split_size"]
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
        "mkdir -p {output}; cat {input} | sort | uniq -c | awk '{{if($1=={params.nfiles}){{print $2}}}}' | "
        "split -l {params.split_size} --numeric-suffixes - {output}/{params.prefix} 1> {log.std} 2>&1"


checkpoint merged_sequences: # get merged sequences by common IDs
    input:
        single_copy_busco_sequences_dir_path / "{N}"
    output:
        directory(merged_sequences_dir_path / "{N}")
    params:
        single_copy_files=expand(busco_dir_path / "{species}" / "single_copy_busco_sequences", species=config["species_list"])
    log:
        std=log_dir_path / "merged_sequences.{N}.log",
        cluster_log=cluster_log_dir_path / "merged_sequences.{N}.cluster.log",
        cluster_err=cluster_log_dir_path / "merged_sequences.{N}.cluster.err"
    benchmark:
        benchmark_dir_path / "merged_sequences.{N}.benchmark.txt"
    resources:
        cpus=config["common_ids_threads"],
        time=config["common_ids_threads"],
        mem=config["common_ids_threads"]
    shell:
        "workflow/scripts/merged_sequences.py "
        "--common_ids {input} "
        "--single_copy_files {params.single_copy_files} "
        "--outdir {output} 1> {log.std} 2>&1"