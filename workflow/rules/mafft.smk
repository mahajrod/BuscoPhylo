localrules: mafft_tasks_list


checkpoint mafft:
    input:
        mafft_task= aggregate_input # mafft_dir_path / 'slurm/mafft.tasks.{i}.sh'
    output:
        mafft_outpath=directory(mafft_dir_path / "output")
    log:
        std=log_dir_path / "mafft.log",
        cluster_log=cluster_log_dir_path / "mafft.cluster.log",
        cluster_err=cluster_log_dir_path / "mafft.cluster.err"
    benchmark:
        benchmark_dir_path / "mafft.benchmark.txt"
    resources:
        cpus=config["mafft_threads"],
        time=config["mafft_time"],
        mem=config["mafft_mem_mb"],
    shell:
        "bash {input.mafft_task} > {log.std} 2>&1 "


checkpoint mafft_tasks_list:
    input:
        merged_ids=directory(busco_dir_path / "merged_sequences")
    output:
        mafft_tasks=directory(mafft_dir_path / "slurm")
    params:
        amount_of_tasks = 20,
        file_extension = "faa",
        mafft_command_outdir = mafft_dir_path / "output"
    log:
        std=log_dir_path / "mafft_tasks_list.log",
        cluster_log=cluster_log_dir_path / "mafft_tasks_list.cluster.log",
        cluster_err=cluster_log_dir_path / "mafft_tasks_list.cluster.err"
    benchmark:
        benchmark_dir_path / "mafft_tasks_list.benchmark.txt"
    resources:
        cpus=config["common_ids_threads"],
        time=config["common_ids_threads"],
        mem=config["common_ids_threads"]
    shell:
        "workflow/scripts/mafft_tasks_list.py "
        "--input {input.merged_ids} "
        "--file-extension {params.file_extension} "
        "--amount {params.amount_of_tasks} "
        "--mafft_command_outdir {params.mafft_command_outdir} "
        "--outdir {output.mafft_tasks} > {log.std} 2>&1"