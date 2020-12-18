"""
Snakemake RNA Secondary Structure Prediction Pipeline
KH (Nov 2020)
"""

# output files
tidy_gff = config['reference']['gff'].replace('.gff3', '.tidy.gff3')
intron_gff = tidy_gff.replace('.gff3', '.incl.introns.gff3')
combined_fasta = config['reference']['fasta'].replace('.fa', '.features.fa')

out_dir = os.path.join(config['out_dir'], config['version'])
rnafold_outputs = os.path.join(out_dir, os.path.basename(combined_fasta).replace(".features.fa", ".{feature}.rnafold.out"))

# feature types and corresponding outputs
feature_types = ['five_prime_UTR', 'three_prime_UTR', 'exon', 'intron']
feature_fastas = combined_fasta.replace(".features.fa", ".{feature}.fa")

#
# rules
#
rule all:
    input:
        expand(rnafold_outputs, feature=feature_types)

# run rnafold on feature sequences
rule run_rnafold:
    input:
        feature_fastas
    output:
        rnafold_outputs
    shell:
        """
        ps_dir="`dirname {output}`/ps/{wildcards.feature}"

        mkdir -p $ps_dir
        cd $ps_dir

        RNAfold --verbose -i {input} -g > {output}
        """

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
            rm "$prefix.$x.fa"
        done
        """

# parse chromosome-level fasta into individual exon, intron, 5'utr, etc. sequences
rule extract_feature_seqs:
    input:
        fasta=config['reference']['fasta'],
        gff=intron_gff
    output:
        combined_fasta
    shell:
        "bedtools getfasta -fi {input.fasta} -bed {input.gff} -split -name | fold -w 100 > {output}"

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
