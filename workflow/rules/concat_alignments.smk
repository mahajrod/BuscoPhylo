localrules: concat_fasta_dna, concat_fasta_protein, concat_nexus_dna, concat_nexus_protein


rule concat_fasta_dna:
    input:
        directory(trimal_dir_path / "fna")
    output:
        concat_aligments_dir_path / config["concat_fna_filename"]
    log:
        std=log_dir_path / "fna.concat_fasta_dna.log"
    benchmark:
        benchmark_dir_path / "fna.concat_fasta_dna.benchmark.txt"
    shell:
        "cat {input}/*.fna | workflow/scripts/concat_fasta.py -o {output} 1> {log.std} 2>&1"


rule concat_fasta_protein:
    input:
        directory(trimal_dir_path / "faa")
    output:
        concat_aligments_dir_path / config["concat_faa_filename"]
    log:
        std=log_dir_path / "faa.concat_fasta_protein.log"
    benchmark:
        benchmark_dir_path / "faa.concat_fasta_protein.benchmark.txt"
    shell:
        "cat {input}/*.faa | workflow/scripts/concat_fasta.py -o {output} 1> {log.std} 2>&1"


rule concat_nexus_dna:
    input:
        concat_aligments_dir_path / config["concat_fna_filename"]
    output:
        concat_aligments_dir_path / "concat.aln.fna.nexus"
    params:
        type="DNA",
        block="../../%s" % {config["mrbayes_block"]}
    log:
        std=log_dir_path / "fna.concat_nexus_dna.log"
    benchmark:
        benchmark_dir_path / "fna.concat_nexus_dna.benchmark.txt"
    shell:
        "workflow/scripts/fasta_to_nexus.py -i {input} -t {params.type} -b {params.block} -o {output} 1> {log.std} 2>&1"


rule concat_nexus_protein:
    input:
        concat_aligments_dir_path / config["concat_faa_filename"]
    output:
        concat_aligments_dir_path / "concat.aln.faa.nexus"
    params:
        type="protein",
        block="../../%s" % {config["mrbayes_block"]}
    log:
        std=log_dir_path / "faa.concat_nexus_protein.log"
    benchmark:
        benchmark_dir_path / "faa.concat_nexus_protein.benchmark.txt"
    shell:
        "workflow/scripts/fasta_to_nexus.py -i {input} -t {params.type} -b {params.block} -o {output} 1> {log.std} 2>&1"