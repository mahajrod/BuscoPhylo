localrules: merged_sequences

rule merged_sequences:
    input:
        unique_ids=busco_dir_path / "single_copy_busco_sequences.uniq.ids",
        single_copy_dir=expand(busco_dir_path / "{species}/single_copy_busco_sequences", species=config["species_list"])
    output:
        merged_ids=directory(busco_dir_path / "merged_sequences")
    # params:
    #     nfiles=len(config["species_list"])
    log:
        std=log_dir_path / "merged_ids.log",
        cluster_log=cluster_log_dir_path / "merged_ids.cluster.log",
        cluster_err=cluster_log_dir_path / "merged_ids.cluster.err"
    benchmark:
        benchmark_dir_path / "merged_ids.benchmark.txt"
    resources:
        cpus=config["unique_ids_threads"],
        time=config["unique_ids_time"],
        mem=config["unique_ids_mem_mb"]
    shell:
        "workflow/scripts/merged_sequences.py "
        "--input {input.unique_ids} "
        "--species-directories {input.single_copy_dir} "
        "--outdir {output.merged_ids} "