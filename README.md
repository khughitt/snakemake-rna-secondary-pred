RNA Secondary Structure Prediction Pipeline
===========================================

Overview
--------

This repository contains a simple pipeline for performing RNA secondary structure
prediction for individual gene features (introns, exons, 5' and 3'UTRs) for an entire
genome.

It requires a reference genome and annotations, and is currently designed to work with
those provided by [Gencode](https://www.gencodegenes.org/human/).

First, [GenomeTools](http://genometools.org/tools.html) is used to clean up the Gencode
GFF file and fix some small issues that prevent it from being used in later steps in the
pipline..

Next, GenomeTools is used to infer the locations of introns, which are not included in
the original Gencode GFF files, and create a modified version of the GFF file including
these features.

[bedtools](https://bedtools.readthedocs.io/en/latest/) is then used to create a fasta
file with sequences for each intron, exon, etc. in the genome, and this file is split
into several files, each containing sequences associated with one feature type.

Finally, [RNAfold](https://www.tbi.univie.ac.at/RNA/RNAfold.1.html) is used to infer
secondary structure for each of these sequences. In the future, other structure
prediction methods may be included as well.

Note: in order to prevent RNAfold from filling up the directory with postscript
structure plots, RNAfold is currently executed with an option to disable postscript
generation.

Usage
-----

To begin, create and activate conda environment with the needed requirements:

```
conda create -n rna2pred --file requirements.txt
conda activate rna2pred
```

Next, download the latest version of human reference genome and annotations from the
[gencode FTP site](ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/). A helper
script is provided in the `extra/` directory which can be used for this purpose.

Modify `config/config-v0.1.yml` to point to the locations of the reference files.

To run the pipeline, one can then do:

```
snakemake --configfile=config/config-v0.1.yml -j4
```

The above command will run the pipeline with 4 threads.

Issues
------

- Currently, RNAfold fails for some sequences with an error along the lines of
  "ERROR: Unexpected large magnitude discriminant..."; this causes the snakemake rule to
  fail for those features.
  - RNAfold doesn't appear to have any option to ignore such problematic sequences.
  - As such, I will either need to track down and exclude the problematic sequences 
    manually.
  - It may be helpful to split the feature-specific fasta files (e.g.
    `GRCh38.primary_assembly.genome.five_prime_UTR.fa`) into multiple sub-files, each
    containing the sequence for a single feature.
- The pipeline currently generates fasta files with entry headers corresponding to the
  _coordinates_ of the features in the original fasta file, e.g.
  ">five_prime_UTR::chr1:65418-65433".
  - It may be useful modify this to include the gene identifiers as well.

