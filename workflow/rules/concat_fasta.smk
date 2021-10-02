localrules: concat_fasta_dna, concat_fasta_protein


rule concat_fasta_dna:
    input:
        directory(trimal_dir_path / "fna")
    output:
        concat_aligments_dir_path / config["concat_fna_filename"]
    log:
        std=log_dir_path / "fna.concat_fasta.log",
        cluster_log=cluster_log_dir_path / "fna.concat_fasta.log",
        cluster_err=cluster_log_dir_path / "fna.concat_fasta.err"
    benchmark:
        benchmark_dir_path / "fna.concat_fasta.benchmark.txt"
    shell:
        "cat {input}/*.fna | workflow/scripts/concat_fasta.py -o {output} 1> {log.std} 2>&1"


rule concat_fasta_protein:
    input:
        directory(trimal_dir_path / "faa")
    output:
        concat_aligments_dir_path / config["concat_faa_filename"]
    log:
        std=log_dir_path / "faa.concat_fasta.log",
        cluster_log=cluster_log_dir_path / "faa.concat_fasta.log",
        cluster_err=cluster_log_dir_path / "faa.concat_fasta.err"
    benchmark:
        benchmark_dir_path / "faa.concat_fasta.benchmark.txt"
    shell:
        "cat {input}/*.faa | workflow/scripts/concat_fasta.py -o {output} 1> {log.std} 2>&1"