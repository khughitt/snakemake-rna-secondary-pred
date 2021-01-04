#!/bin/env Rscript
#
# Creates a list of Ensembl gene identifiers corresponding to MM25 genes of interest
#
library(biomaRt)
library(readr)
library(annotables)

# extract top 100 genes from MM25
mm25 <- read_tsv(snakemake@config[['mm25']], col_types = cols())
ensgenes <- grch38$ensgene[match(head(mm25$symbol, snakemake@config[['num_genes']]), grch38$symbol)]

# convert to transcript ids;
# required to handle intron gff entries which only have tx ids associated with them
ensembl <- useEnsembl(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl")

res <- getBM(attributes = c('ensembl_gene_id', 'ensembl_transcript_id'), 
             filters = "ensembl_gene_id",
             values = ensgenes, 
             mart = ensembl)

enstranscripts = res[, 'ensembl_transcript_id']

writeLines(enstranscripts, snakemake@output[[1]])
