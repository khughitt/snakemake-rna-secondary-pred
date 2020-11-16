#!/bin/env bash
#
# Download human reference genome sequence and annotations from Gencode
# KH (Nov 2020)
#
REL="35"

# make output dir
OUTDIR="/data/ref/gencode/human/v$REL"

mkdir -p $OUTDIR/genome/gff3
mkdir -p $OUTDIR/genome/fasta

# ftp base url
BASE_URL="ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_$REL"

# gff
GFF="$OUTDIR/genome/gff3/gencode.v$REL.annotation.gff3.gz"

if [ ! -f $GFF ]; then
    echo "Downloading $GFF..."
    curl "$BASE_URL/gencode.v$REL.annotation.gff3.gz" -o $GFF 
fi

# genome fasta
GENOME_FASTA="$OUTDIR/genome/fasta/GRCh38.primary_assembly.genome.fa.gz"

if [ ! -f $GENOME_FASTA ]; then
    echo "Downloading $GENOME_FASTA..."
    curl "$BASE_URL/GRCh38.primary_assembly.genome.fa.gz" -o $GENOME_FASTA
fi

# checksums
curl "$BASE_URL/MD5SUMS" -o $OUTDIR/MD5SUMS
