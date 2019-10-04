# Snakemake rules for quantifying mRNA transcript counts

def get_feature_counts_tsv(wildcards):
    path = os.path.join("genomes/expression/gene_counts", "{genome}.tsv")
    feature_counts_tsv = expand(path, genome=get_genomes_(wildcards))

    return(feature_counts_tsv)

def feature_counts_parameters(wc,input):
    params={}

    if (str(config['feature_counts']['non_overlap']).lower() != 'no_limit'):
        params['non_overlap'] = "--nonOverlap {p}".format(p=config['feature_counts']['non_overlap'])
    else:
        params['non_overlap'] = ""

    if (str(config['feature_counts']['non_overlap_feature']).lower() != 'no_limit'):
        params['non_overlap_feature'] = "--nonOverlapFeature {p}".format(p=config['feature_counts']['non_overlap_feature'])
    else:
        params['non_overlap_feature'] = ""

    if (int(config['feature_counts']['read_extension_5']) != 0):
        params['read_extension_5'] = "--readExtension5 {p}".format(p=config['feature_counts']['read_extension_5'])
    else:
        params['read_extension_5'] = ""

    if (int(config['feature_counts']['read_extension_3']) != 0):
        params['read_extension_3'] = "--readExtension3 {p}".format(p=config['feature_counts']['read_extension_3'])
    else:
        params['read_extension_3'] = ""

    if (str(config['feature_counts']['read_to_pos']).lower() != 'na'):
        params['read_to_pos'] = "--read2pos {p}".format(p=config['feature_counts']['read_to_pos'])
    else:
        params['read_to_pos'] = ""

    params['split_only'] = "--splitOnly" if config['feature_counts']['split_only'] else ""
    params['non_split_only'] = "--nonSplitOnly" if config['feature_counts']['non_split_only'] else ""
    params['ignore_dup'] = "--ignoreDup" if config['feature_counts']['ignore_dup'] else ""
    params['need_both_paired_reads'] = "-B" if config['feature_counts']['need_both_paired_reads'] else ""

    return params

## FeatureCounts per MAG based on prodigal output
rule run_feature_counts:
    input:
        gff = "genomes/annotations/genes/{genome}.gff",
        samples_bam = expand("genomes/alignments/{sample}.bam", sample=SAMPLES)
    output:
        counts = temp(expand("genomes/expression/gene_counts/{{genome}}_raw.{ext}",
                              ext=['tsv','tsv.summary']))
        # TODO: Consider saving some of the content of the summary file
    threads:
        config['threads']
    log:
        "logs/genomes/feature_counts/{genome}.log"
    conda:
        "%s/required_packages.yaml" % CONDAENV
    shadow:
        "shallow"
    params:
        p = lambda wc,input: feature_counts_parameters(wc,input),
        min_overlap = config['feature_counts']['min_overlap'],
        frac_overlap = config['feature_counts']['frac_overlap'],
        frac_overlap_feature = config['feature_counts']['frac_overlap_feature'],
        min_mapping_qual_score = config['feature_counts']['min_mapping_qual_score'],
        max_m_cigar_operations = config['feature_counts']['max_m_cigar_operations']
    shell:
        """
            featureCounts \
            -a {input.gff} \
            -o {output.counts[0]} \
            -F GTF \
            -t CDS \
            -g ID \
            --minOverlap {params.min_overlap} \
            --fracOverlap {params.frac_overlap} \
            --fracOverlapFeature {params.frac_overlap_feature} \
            {params.p[non_overlap]} \
            {params.p[non_overlap_feature]} \
            {params.p[read_extension_5]} \
            {params.p[read_extension_3]} \
            {params.p[read_to_pos]} \
            -Q {params.min_mapping_qual_score} \
            {params.p[split_only]} \
            {params.p[non_split_only]} \
            {params.p[ignore_dup]} \
            {params.p[need_both_paired_reads]} \
            --maxMOp {params.max_m_cigar_operations} \
            --donotsort \
            -T {threads} \
            {input.samples_bam} \
            2> {log}
        """

# Now parse the output
rule parse_feature_counts_output:
    input:
        "genomes/expression/gene_counts/{genome}_raw.tsv"
    output:
        "genomes/expression/gene_counts/{genome}.tsv"
    threads:
        config['threads']
    run:
        """ INPUT FORMAT
        # Program:featureCounts v1.6.4; Command:"featureCounts" "-a" "MAG1.gff" "-o" "test.mag1.tsv" "-t" "CDS" "-g" "ID" "../../alignments/Cfx_borealis.bam"
        Geneid  Chr     Start   End     Strand  Length  genomes/alignments/sample_1.bam genomes/alignments/sample_2.bam
        1_1     MAG1_1  142     894     +       753     277                             437
        1_2     MAG1_1  979     2151    +       1173    488                             985
        1_3     MAG1_1  2182    3450    +       1269    500                             236
        1_4     MAG1_1  3492    4235    -       744     364                             2478
        1_5     MAG1_1  4386    5066    -       681     336                             347
        """

        """ DESIRED OUTPUT FORMAT
        ORF          contig  gene_length  sample_1 sample_2
        MAG1_1_1     MAG1_1  753     277      437
        MAG1_1_2     MAG1_1  1173    488      985
        MAG1_1_3     MAG1_1  1269    500      236
        MAG1_1_4     MAG1_1  744     364      2478
        MAG1_1_5     MAG1_1  681     336      347
        """

        import pandas as pd
        import os
        from multiprocessing.dummy import Pool

        def generate_ORF_ID(row_tuple):
            # Converts a row tuple from iterrows (of a featureCounts table) into a ORF ID (string)
            table_row = row_tuple[1]
            genome = str(table_row[1]).split(sep='_')[0]
            ORF_ID = "{genome}_{contig_gene}".format(genome=genome,contig_gene=table_row[0])
            return(ORF_ID)

        feature_counts_table = pd.read_csv(input[0], sep='\t', header=1)

        for bam_path in feature_counts_table.columns.values.tolist()[6:]:
            sample_ID = os.path.splitext(os.path.basename(bam_path))[0]
            feature_counts_table.rename(columns = {bam_path: sample_ID}, inplace = True)

        # Combine the MAG/contig and gene IDs together in parallel
        p = Pool(threads)
        ORF_IDs = p.map(generate_ORF_ID, feature_counts_table.iterrows())

        feature_counts_table.insert(loc=0, column='ORF', value=ORF_IDs)

        feature_counts_table.rename(columns = {'Length': 'length_bp'}, inplace = True)
        feature_counts_table.drop(columns = ['Geneid', 'Chr', 'Start', 'End', 'Strand'], inplace = True)

        pd.DataFrame.to_csv(feature_counts_table, output[0], sep = '\t', index = False)

# TODO: this is a dummy summary file for now. Make a proper summary file with processed expression data
rule summarize_gene_counts:
    input:
        get_feature_counts_tsv
    output:
        touch("genomes/expression/gene_counts.tsv")
