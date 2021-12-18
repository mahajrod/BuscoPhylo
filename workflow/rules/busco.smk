if config['busco_version'] == 3:
    rule busco3:
        input:
            fasta=genome_dir_path / "{species}.fasta"
        output:
            busco_outdir=directory(busco_dir_path / "{species}"),
            single_copy_files_dir = directory(busco_dir_path / "{species}/single_copy_busco_sequences"),
            summary=busco_dir_path / "{species}/short_summary_{species}.txt"
        params:
            busco_path=config["busco_path"],
            mode=config["busco_mode"],
            species=config["augustus_species"],
            busco_dataset_path=config["busco_dataset_path"],
            output_prefix="{species}"
        log:
            std=log_dir_path / "busco/{species}.log",
            cluster_log=cluster_log_dir_path / "busco/{species}.cluster.log",
            cluster_err=cluster_log_dir_path / "busco/{species}.cluster.err"
        benchmark:
            benchmark_dir_path / "busco/{species}.benchmark.txt"
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
            single_copy_files_dir = temp(directory(busco_dir_path / "{species}/busco_sequences/single_copy_busco_sequences")),
            summary=busco_dir_path / "{species}/short_summary_{species}.txt"
        params:
            mode=config["busco_mode"],
            busco_dataset_path=config["busco_dataset_path"],
            output_prefix="{species}"
        log:
            std=log_dir_path / "busco/{species}.log",
            cluster_log=cluster_log_dir_path / "busco/{species}.cluster.log",
            cluster_err=cluster_log_dir_path / "busco/{species}.cluster.err"
        benchmark:
            benchmark_dir_path / "busco/{species}.benchmark.txt"
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
            "mv {params.output_prefix}/* ./; rm -r {params.output_prefix}/ 1>../../../{log.std} 2>&1; "
            "mv run*/* ./; rm -r run* 1>../../../{log.std} 2>&1; "
            "mv full_table.tsv full_table_{params.output_prefix}.tsv 1>../../../{log.std} 2>&1; "
            "mv missing_busco_list.tsv missing_busco_list_{params.output_prefix}.tsv 1>../../../{log.std} 2>&1; "
            "mv short_summary.txt short_summary_{params.output_prefix}.txt 1>../../../{log.std} 2>&1; "


    rule get_fna_sequences:
        input:
            single_copy_files_dir=temp(directory(busco_dir_path / "{species}/busco_sequences/single_copy_busco_sequences")),
            fasta=genome_dir_path / "{species}.fasta"
        output:
            directory(busco_dir_path / "{species}/single_copy_busco_sequences")
        log:
            std=log_dir_path / "get_fna_sequences/{species}.log",
            cluster_log=cluster_log_dir_path / "get_fna_sequences/{species}.cluster.log",
            cluster_err=cluster_log_dir_path / "get_fna_sequences/{species}.cluster.err"
        benchmark:
            benchmark_dir_path / "get_fna_sequences/{species}.benchmark.txt"
        conda:
            config["conda_config"]
        resources:
            cpus=config["common_ids_threads"],
            time=config["common_ids_time"],
            mem=config["common_ids_mem_mb"],
        threads:
            config["common_ids_threads"]
        shell:
            "mkdir -p {output}; "
            "for FILE in `ls {input.single_copy_files_dir}`; do "
            "HEADER=$(head -n 1 {input.single_copy_files_dir}/$FILE | sed 's/^.//'); " # without '>'
            "samtools faidx {input.fasta} $HEADER >> {output}/${{FILE%.*}}.fna 2> {log.std}; "
            "mv {input.single_copy_files_dir}/$FILE {output}/; "
            "done"
else:
    print("Specify the version of BUSCO in 'busco_version' parameter!")