if config["alignment_tool"] == "prank":
    rule prank_dna:
        input:
            merged_sequences_dir_path / "{N}"
        output:
            temp(directory(alignment_dir_path / "fna_tmp" / "{N}"))
        params:
            options=config["prank_dna_params"]
        log:
            std=log_dir_path / "prank_dna.{N}.log",
            cluster_log=cluster_log_dir_path / "prank_dna.{N}.cluster.log",
            cluster_err=cluster_log_dir_path / "prank_dna.{N}.cluster.err"
        benchmark:
            benchmark_dir_path / "prank_dna.{N}.benchmark.txt"
        conda:
            "../../%s" % config["conda_config"]
        resources:
            cpus=config["prank_threads"],
            time=config["prank_time"],
            mem=config["prank_mem_mb"]
        threads:
            config["prank_threads"]
        shell:
            "mkdir -p {output}; "
            "for FILE in `ls {input}/*.fna`; do "
            "prank -d=$FILE -o={output}/$(basename $FILE) -codon -F > {log.std} 2>&1; "
            "mv {output}/$(basename $FILE).best.fas {output}/$(basename $FILE); "
            "done; "


    rule prank_protein:
        input:
            merged_sequences_dir_path / "{N}"
        output:
            temp(directory(alignment_dir_path / "faa_tmp" / "{N}"))
        params:
            options=config["prank_protein_params"]
        log:
            std=log_dir_path / "prank_protein.{N}.log",
            cluster_log=cluster_log_dir_path / "prank_protein.{N}.cluster.log",
            cluster_err=cluster_log_dir_path / "prank_protein.{N}.cluster.err"
        benchmark:
            benchmark_dir_path / "prank_protein.{N}.benchmark.txt"
        conda:
            "../../%s" % config["conda_config"]
        resources:
            cpus=config["prank_threads"],
            time=config["prank_time"],
            mem=config["prank_mem_mb"]
        threads:
            config["prank_threads"]
        shell:
            "mkdir -p {output}; cd {output}; "
            "for FILE in `ls ../../../../../{input}/*.fna`; do "
            "prank -d=$FILE -o=$(basename $FILE) -translate -F > {log.std} 2>&1; "
            "mv $(basename $FILE).best.pep.fas $(basename $FILE); "
            "rm $(basename $FILE).best.nuc.fas; "
            "done; "

        # "mkdir -p {output}; cd {output}; "
        # "for FILE in `ls ../../../../../{input}/*.fna`; do "
        # "prank -d=$FILE -o={output}/$(basename $FILE) -translate -F > {log.std} 2>&1; "
        # "mv {output}/$(basename $FILE).best.pep.fas {output}/$(basename $FILE); "
        # "rm {output}/$(basename $FILE).best.nuc.fas; "
        # "done; "