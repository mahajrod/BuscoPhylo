localrules: mafft_tasks_list

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
        "for task in `ls {input.mafft_tasks}`; do"
        "bash $task; "
        "done"


rule mafft_tasks_list:
    input:
        merged_ids=directory(busco_dir_path / "merged_sequences")
    output:
        mafft_tasks=directory(mafft_dir_path / "slurm")
    params:
        number_of_tasks = 50, 
        merged_ids_path = busco_dir_path / "merged_sequences",
        mafft_outpath = mafft_dir_path / "output"
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
        "counter=0; filename='mafft_task'; "
        "for i in `ls {input.merged_ids}.faa`; do "
        "(( counter++ )); "
        "echo -e \"mafft --anysymbol {params.merged_ids_path}/$i > {params.mafft_outpath}mafft.$i\" >> {output.mafft_tasks}/${filename}.sh; "
        "if [ $[$C % {params.number_of_tasks}] -eq '0']; "
        "then filename=mafft_task_${C}.sh; fi; "
        "done"