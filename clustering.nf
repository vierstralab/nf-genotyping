#!/usr/bin/env nextflow
nextflow.enable.dsl = 2


process clusterIndivs {

    publishDir "${params.outdir}/clustering"

    input:
        path vcf_file
    output:
        tuple path('metadata.clustered.tsv'), path('clustering.png')
    script:
    """
    plink2 --allow-extra-chr \
    --make-king square \
    --out snps.clustering \
    --vcf ${vcf_file}

    python3 $baseDir/bin/cluster_king.py --matrix snps.clustering.king \
    --matrix-ids snps.clustering.king.id \
    --meta-file ${params.samples_file} \
    --outpath ./
    """
}

workflow clustering {
    take:
        vcfs_and_index
    main:
        clusterIndivs(vcf_and_index)
    emit:
        clusterIndivs.out
}

workflow {
    vcfs_and_index = Channel.value(params.genotype_file)
    clustering(vcfs_and_index)
}