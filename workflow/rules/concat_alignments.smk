localrules: concat_fasta_dna, concat_fasta_protein, concat_nexus_dna, concat_nexus_protein


rule concat_fasta_dna:
    input:
        directory(trimal_dir_path / "fna")
    output:
        concat_alignments_dir_path / fasta_dna_filename
    log:
        std=log_dir_path / "concat_fasta_dna.log",
        cluster_log=cluster_log_dir_path / "concat_fasta_dna.cluster.log",
        cluster_err=cluster_log_dir_path / "concat_fasta_dna.cluster.err"
    benchmark:
        benchmark_dir_path / "concat_fasta_dna.benchmark.txt"
    shell:
        "cat {input}/*.fna | workflow/scripts/concat_fasta.py -o {output} 1> {log.std} 2>&1"


rule concat_fasta_protein:
    input:
        directory(trimal_dir_path / "faa")
    output:
        concat_alignments_dir_path / fasta_protein_filename
    log:
        std=log_dir_path / "concat_fasta_protein.log",
        cluster_log=cluster_log_dir_path / "concat_fasta_protein.cluster.log",
        cluster_err=cluster_log_dir_path / "concat_fasta_protein.cluster.err"
    benchmark:
        benchmark_dir_path / "concat_fasta_protein.benchmark.txt"
    shell:
        "cat {input}/*.faa | workflow/scripts/concat_fasta.py -o {output} 1> {log.std} 2>&1"


rule concat_nexus_dna:
    input:
        concat_alignments_dir_path / fasta_dna_filename
    output:
        concat_alignments_dir_path / nexus_dna_filename
    params:
        type="DNA",
        block=config["mrbayes_block"]
    log:
        std=log_dir_path / "concat_nexus_dna.log",
        cluster_log=cluster_log_dir_path / "concat_nexus_dna.cluster.log",
        cluster_err=cluster_log_dir_path / "concat_nexus_dna.cluster.err"
    benchmark:
        benchmark_dir_path / "concat_nexus_dna.benchmark.txt"
    shell:
        "workflow/scripts/fasta_to_nexus.py -i {input} -t {params.type} -b {params.block} -o {output} 1> {log.std} 2>&1"


rule concat_nexus_protein:
    input:
        concat_alignments_dir_path / fasta_protein_filename
    output:
        concat_alignments_dir_path / nexus_protein_filename
    params:
        type="protein",
        block=config["mrbayes_block"]
    log:
        std=log_dir_path / "concat_nexus_protein.log",
        cluster_log=cluster_log_dir_path / "concat_nexus_protein.cluster.log",
        cluster_err=cluster_log_dir_path / "concat_nexus_protein.cluster.err"
    benchmark:
        benchmark_dir_path / "concat_nexus_protein.benchmark.txt"
    shell:
        "workflow/scripts/fasta_to_nexus.py -i {input} -t {params.type} -b {params.block} -o {output} 1> {log.std} 2>&1"