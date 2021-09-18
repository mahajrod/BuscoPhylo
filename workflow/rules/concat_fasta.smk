localrules: concat_fasta_dna


rule concat_fasta_dna:
    input:
        directory(trimal_dir_path / "fna" / "{N}")
    output:
        output_dir_path / "concat.fna"
    log:
        std=log_dir_path / "{N}.concat_fasta_dna.log",
        cluster_log=cluster_log_dir_path / "{N}.concat_fasta_dna.log",
        cluster_err=cluster_log_dir_path / "{N}.concat_fasta_dna.err"
    benchmark:
        benchmark_dir_path / "{N}.concat_fasta_dna.benchmark.txt"
    shell:
        ""

rule concat_fasta_protein:
    input:
        directory(trimal_dir_path / "faa" / "{N}")
    output:
        output_dir_path / "concat.faa"
    log:
        std=log_dir_path / "{N}.concat_fasta_protein.log",
        cluster_log=cluster_log_dir_path / "{N}.concat_fasta_protein.log",
        cluster_err=cluster_log_dir_path / "{N}.concat_fasta_protein.err"
    benchmark:
        benchmark_dir_path / "{N}.concat_fasta_protein.benchmark.txt"
    shell:
        ""