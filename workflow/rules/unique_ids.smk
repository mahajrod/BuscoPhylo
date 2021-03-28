localrules: unique_ids, species_ids

rule unique_ids:
    input:
        id_files=expand(busco_dir_path/ "{species}.ids", species=config["species_list"])
    output:
        unique_ids=busco_dir_path / "single_copy_busco_sequences.uniq.ids"
    params:
        nfiles=len(config["species_list"])
    log:
        std=log_dir_path / "unique_ids.log",
        cluster_log=cluster_log_dir_path / "unique_ids.cluster.log",
        cluster_err=cluster_log_dir_path / "unique_ids.cluster.err"
    benchmark:
        benchmark_dir_path / "unique_ids.benchmark.txt"
    resources:
        cpus=config["unique_ids_threads"],
        time=config["unique_ids_time"],
        mem=config["unique_ids_mem_mb"]
    shell:
        "cat {input.id_files} |"
        "sort | uniq -c | awk '{{if($1=={params.nfiles}){{print $2}}}}' > {output.unique_ids} 2> {log.std}"


rule species_ids:
    input:
        single_copy_dir=directory(busco_dir_path / "{species}/single_copy_busco_sequences")
    output:
        ids=temp(busco_dir_path / "{species}.ids")
    log:
        std=log_dir_path / "{species}.species_ids.log",
        cluster_log=cluster_log_dir_path / "{species}.species_ids.cluster.log",
        cluster_err=cluster_log_dir_path / "{species}.species_ids.cluster.err"
    benchmark:
        benchmark_dir_path / "{species}/species_ids.benchmark.txt"
    resources:
        cpus=config["unique_ids_threads"],
        time=config["unique_ids_time"],
        mem=config["unique_ids_mem_mb"]
    shell:
        "ls {input.single_copy_dir}/*.fna | sed 's/.fna//' > {output.ids} 2> {log.std}"