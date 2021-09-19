localrules: concat_fasta_dna, concat_fasta_protein


rule concat_fasta_dna:
    input:
        directory(trimal_dir_path / "fna")
    output:
        output_dir_path / "concat.fna"
    log:
        std=log_dir_path / "concat_fasta_dna.log",
        cluster_log=cluster_log_dir_path / "concat_fasta_dna.log",
        cluster_err=cluster_log_dir_path / "concat_fasta_dna.err"
    benchmark:
        benchmark_dir_path / "concat_fasta_dna.benchmark.txt"
    shell:
        "cat {input}/*.fna | workflow/scripts/merged_sequences.py -o {output}"

rule concat_fasta_protein:
    input:
        directory(trimal_dir_path / "faa")
    output:
        output_dir_path / "concat.faa"
    log:
        std=log_dir_path / "concat_fasta_protein.log",
        cluster_log=cluster_log_dir_path / "concat_fasta_protein.log",
        cluster_err=cluster_log_dir_path / "concat_fasta_protein.err"
    benchmark:
        benchmark_dir_path / "concat_fasta_protein.benchmark.txt"
    shell:
        "cat {input}/*.faa | workflow/scripts/merged_sequences.py -o {output}"