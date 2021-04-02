localrules: mafft_tasks_list, mafft

rule mafft:
    input:
        mafft_tasks=directory(mafft_dir_path / "slurm")
    output:
        mafft_tasks=directory(mafft_dir_path / "output")
    log:
        std=log_dir_path / "mafft.log",
        cluster_log=cluster_log_dir_path / "mafft.cluster.log",
        cluster_err=cluster_log_dir_path / "mafft.cluster.err"
    benchmark:
        benchmark_dir_path / "mafft.benchmark.txt"
    resources:
        cpus=config["common_ids_threads"],
        time=config["common_ids_threads"],
        mem=config["common_ids_threads"]
    shell:
        "for task in `ls {input.mafft_tasks}`; do ./$task; done"


rule mafft_tasks_list:
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
        "--outdir {output.mafft_tasks}"