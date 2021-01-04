"""
Snakemake RNA Secondary Structure Prediction Pipeline
KH (Nov 2020)
"""

# output files
tidy_gff = config['reference']['gff'].replace('.gff3', '.tidy.gff3')
intron_gff = tidy_gff.replace('.gff3', '.incl.introns.gff3')

out_dir = os.path.join(config['out_dir'], config['version'])

mm25_genes = os.path.join(out_dir, 'mm25_genes.txt')
mm25_gff   = os.path.join(out_dir, 'gff', 'mm25_top_genes.gff')

combined_fasta = os.path.join(out_dir, "fasta", 
        os.path.basename(config['reference']['fasta']).replace('.fa', '.features.fa'))
rnafold_outputs = os.path.join(out_dir, "rnafold", 
        os.path.basename(combined_fasta).replace(".features.fa", ".{feature}.rnafold.out"))

# feature types and corresponding outputs
feature_types = ['five_prime_UTR', 'three_prime_UTR', 'exon', 'intron']
feature_fastas = combined_fasta.replace(".features.fa", ".{feature}.fa")
feature_fastas_filtered = combined_fasta.replace(".features.fa", ".{feature}.filtered.fa")

#
# rules
#
rule all:
    input:
        expand(rnafold_outputs, feature=feature_types)

# run rnafold on feature sequences
rule run_rnafold:
    input:
        feature_fastas_filtered
    output:
        rnafold_outputs
    shell:
        """
        ps_dir="`dirname {output}`/ps/{wildcards.feature}"

        mkdir -p $ps_dir
        cd $ps_dir

        RNAfold --verbose -t 4 -i {input} -g > {output}
        """

# filter multifasta files to remove sequences longer than the specified amount
rule filter_fasta:
    input:
        feature_fastas
    output:
        feature_fastas_filtered
    shell:
        """
        cat {{input}} | seqkit seq -m{} -M{} > {{output}} 
        """.format(config['filtering']['min_len'], config['filtering']['max_len'])

# split combined fasta file into feature-specific sub-files
# https://www.biostars.org/p/392018/#392021
rule split_fasta:
    input:
        combined_fasta
    output:
        expand(combined_fasta.replace('.features.fa', '.{feature}.fa'), feature=feature_types)
    shell:
        """
        awk '{{if(substr($0,1,1) == \">\"){{split(substr($0,2,length($0)),a,/::/);outfile=a[1]}};print $0 > outfile \".tmp.fa\" }}' {input}

        prefix={input}
        prefix=${{prefix/.features.fa/}}

        for x in *.tmp.fa; do
            mv $x $prefix.${{x/.tmp.fa/.fa}}
        done

        # delete unneeded files
        for x in "CDS" "gene" "start_codon" "stop_codon" "stop_codon_redefined_as_selenocysteine" "transcript"; do
            rm -f "$prefix.$x.fa"
        done
        """

# parse chromosome-level fasta into individual exon, intron, 5'utr, etc. sequences
rule extract_feature_seqs:
    input:
        fasta=config['reference']['fasta'],
        gff=mm25_gff
    output:
        combined_fasta
    shell:
        "bedtools getfasta -fi {input.fasta} -bed {input.gff} -split -name | fold -w 100 > {output}"

# create a version of the gff with only our genes of interest
rule extract_mm25_genes:
    input:
        gff=intron_gff,
        genes=mm25_genes
    output:
        mm25_gff
    shell:
        """
        grep -f {input.genes} --color='never' {input.gff} > {output}
        """

# builds a list of transcript ids associated with the top mm25 genes
rule get_mm25_top_genes:
    output:
        mm25_genes
    script:
        "script/create_mm25_gene_list.R"

# use genometools to infer introns from gencode gff
rule infer_introns:
    input:
        tidy_gff
    output:
        intron_gff
    shell:
        "gt gff3 -retainids -addintrons -o {output} {input}"

# fix phase in gencode gff3 so that introns can be inferred
rule create_tidy_gff3:
    input:
        config['reference']['gff']
    output:
        tidy_gff
    shell:
        "gt gff3 -sort -tidy -retainids -o {output} {input}"
