# BuscoPhylo

Pipeline to construct species phylogenies using BUSCOs.

![BuscoPhylo pipeline](./pipeline.png)

The pipeline automatically creates the directory structure from the config file. Default directory structure:

```
|-results
    |- busco
    |- ids
    |   |- species_ids
    |   |- common_ids
    |   |- merged_sequences  
    |- alignments
    |   |- raw
    |   |- trimal
    |- concat_alignments
    |- phylogeny
    |   |- iqtree
    |   |- mrbayes
```

### Usage example

```
snakemake --snakefile Snakemake --profile profile/slurm/ --configfile config/default.yaml \
    --config genome_dir="" busco_dataset_path="" trimal_path="" iqtree_path="" mrbayes_path="" \
    --printshellcmds --latency-wait 60
```

### Requirements
* [Snakemake](https://snakemake.github.io/)
* [BUSCO v3.0.2](https://busco-archive.ezlab.org/) or [BUSCO v5.2.2](https://busco.ezlab.org/)
* [MAFFT](https://mafft.cbrc.jp/alignment/software/)
* [PRANK](http://wasabiapp.org/software/prank/)
* [TrimAl](http://trimal.cgenomics.org/)
* [IQ-TREE](http://www.iqtree.org/)
* [MrBayes](https://nbisweden.github.io/MrBayes/index.html)

You can select the BUSCO version (v3.0.2 or v5.2.2). To do this, specify the value of `busco_version` 3 or 5 in config file.
You can also choose MAFFT or PRANK to align sequences by specifying 'mafft' or 'prank' in the `alignment_tool` option in config file. 
PRANK is only used for DNA sequences! Protein sequences are processed with MAFFT in both cases.

`BUSCO v3.0.2`, `TrimAl`, `IQ-TREE` and `MrBayes` paths should be in config file or specified at startup.
`BUSCO v5.2.2`, `MAFFT` and `PRANK` are installed using a conda.
