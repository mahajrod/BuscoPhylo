if config['busco_version'] == 3:
    rule busco3:
        input:
            fasta=genome_dir_path / "{species}.fasta"
        output:
            busco_outdir=directory(busco_dir_path / "{species}"),
            summary=busco_dir_path / "{species}/short_summary_{species}.txt"
        params:
            busco_path=config["busco_path"],
            mode=config["busco_mode"],
            species=config["augustus_species"],
            busco_dataset_path=config["busco_dataset_path"],
            output_prefix="{species}"
        log:
            std=log_dir_path / "{species}.busco.log",
            cluster_log=cluster_log_dir_path / "{species}.busco.cluster.log",
            cluster_err=cluster_log_dir_path / "{species}.busco.cluster.err"
        benchmark:
            benchmark_dir_path / "{species}/busco.benchmark.txt"
        resources:
            cpus=config["busco_threads"],
            time=config["busco_time"],
            mem=config["busco_mem_mb"],
        threads:
            config["busco_threads"]
        shell:
            "mkdir -p {output.busco_outdir}; cd {output.busco_outdir}; {params.busco_path}/run_BUSCO.py -m {params.mode} -sp {params.species}"
            " -i {input.fasta} -c {threads} -l {params.busco_dataset_path} -o {params.output_prefix} 1>../../../{log.std} 2>&1;"
            " mv run_{params.output_prefix}/* ./; rm -r run_{params.output_prefix} tmp/"

elif config['busco_version'] == 5:
    rule busco5:
        input:
            fasta=genome_dir_path / "{species}.fasta"
        output:
            busco_outdir=directory(busco_dir_path / "{species}"),
            summary=busco_dir_path / "{species}/short_summary_{species}.txt"
        params:
            mode=config["busco_mode"],
            busco_dataset_path=config["busco_dataset_path"],
            output_prefix="{species}"
        log:
            std=log_dir_path / "{species}.busco.log",
            cluster_log=cluster_log_dir_path / "{species}.busco.cluster.log",
            cluster_err=cluster_log_dir_path / "{species}.busco.cluster.err"
        benchmark:
            benchmark_dir_path / "{species}/busco.benchmark.txt"
        conda:
            config["conda_config"]
        resources:
            cpus=config["busco_threads"],
            time=config["busco_time"],
            mem=config["busco_mem_mb"],
        threads:
            config["busco_threads"]
        shell:
            "mkdir -p {output.busco_outdir}; cd {output.busco_outdir}; "
            "busco -m {params.mode} -i {input.fasta} -c {threads} -l {params.busco_dataset_path} -o {params.output_prefix} 1>../../../{log.std} 2>&1; "
            "rm -r busco_downloads/; "
            "mv {params.output_prefix}/* ./; rm -r {params.output_prefix}/; "
            "mv run_/* ./; rm -r run_/; "
            "mv short_summary.txt short_summary_{params.output_prefix}.txt; "
            "mv full_table.tsv full_table_{params.output_prefix}.tsv; "
            "mv missing_busco_list.tsv missing_busco_list_{params.output_prefix}.tsv; "
            "ln -s busco_sequences/single_copy_busco_sequences single_copy_busco_sequences; "
else:
    print("Specify the version of BUSCO in 'busco_version' parameter!")