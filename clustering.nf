#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

conda_env = "$moduleDir/environment.yml"

process cluster_indivs {

    publishDir "${params.outdir}/clustering"
    conda conda_env

    output:
        tuple path('metadata.clustered.tsv'), path('clustering.png')
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
        cluster_indivs()
    emit:
        cluster_indivs.out
}

workflow {
    clusterIndivs()
}