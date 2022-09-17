#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

params.conda = "$moduleDir/environment.yml"

process cluster_indivs {

    publishDir "${params.outdir}/clustering"
    conda params.conda

    output:
        tuple path('metadata.clustered.tsv'), path('clustering.png')
        tuple path('snps.clustering.king.id'), path('snps.clustering.king')
    script:
    """
    plink2 --allow-extra-chr \
    --make-king square \
    --out snps.clustering \
    --vcf ${params.genotype_file}

    python3 $moduleDir/bin/cluster_king.py --matrix snps.clustering.king \
    --matrix-ids snps.clustering.king.id \
    --meta-file ${params.samples_file} \
    --outpath ./
    """
}

workflow clusterIndivs {
    main:
        updated_meta = cluster_indivs()[0]
    emit:
        updated_meta
}

workflow {
    clusterIndivs()
}