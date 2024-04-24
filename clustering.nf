#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

params.conda = "$moduleDir/environment.yml"


process cluster_indivs {
    publishDir "${params.outdir}/clustering"
    conda params.conda
    scratch true
    label "highmem"

    input:
        path genotype_file
        path bcftools_stats

    output:
        tuple path('metadata.clustered.tsv'), path('snps.clustering.king.id'), path('snps.clustering.king'), path('clustering.png')

    script:
    """
    plink2 --allow-extra-chr \
        --make-king square \
        --out snps.clustering \
        --vcf ${genotype_file}

    python3 $moduleDir/bin/cluster_king.py \
        snps.clustering \
        ${params.samples_file} \
        ${bcftools_stats} \
        --min-variants ${params.min_variants} \
        --outpath ./
    """
}

workflow clusterIndivs {
    take:
        genotype_file
        bcftools_stats
    main:
        updated_meta = cluster_indivs(genotype_file, bcftools_stats)[0]
    emit:
        updated_meta
}

// Path to resulting genotype and bcftools stats files a.k.a the output of genotyping.nf script (all.filtered.snps.annotated.vcf.gz)
workflow {
    clusterIndivs(params.genotype_file, params.bcftools_stats)
}