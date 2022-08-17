#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

conda_env = "$moduleDir/environment.yml"

process cluster_indivs {

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

    python3 $moduleDir/bin/cluster_king.py --matrix snps.clustering.king \
    --matrix-ids snps.clustering.king.id \
    --meta-file ${params.samples_file} \
    --outpath ./
    """
}

workflow clustering {
    take:
        vcfs_and_index
    main:
        cluster_indivs(vcfs_and_index)
    emit:
        cluster_indivs.out
}

workflow {
    vcfs_and_index = Channel.value(params.genotype_file)
    clustering(vcfs_and_index)
}