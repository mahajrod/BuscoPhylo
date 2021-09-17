localrules: common_ids, species_ids
ruleorder: species_ids > common_ids

rule species_ids:
    input:
        busco_outdir=directory(busco_dir_path / "{species}")
    output:
        ids=busco_dir_path / "{species}.ids"
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
        id_files=expand(busco_dir_path/ "{species}.ids", species=config["species_list"])
    output:
        directory(output_dir_path / "single_copy_busco_sequences.common.dir")
    params:
        nfiles=len(config["species_list"]),
        prefix="single_copy_busco_sequences.common.ids",
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
        "cat {input.id_files} | "
        "sort | uniq -c | awk '{{if($1=={params.nfiles}){{print $2}}}}' | "
        "split -l {params.split_size} --numeric-suffixes - {output}/{params.prefix}"