#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

params.conda = "$moduleDir/environment.yml"


process cluster_indivs {

    publishDir "${params.outdir}/clustering"
    conda params.conda
    scratch true

    input:
        path genotype_file
        path bcftools_stats

    output:
        path('metadata.clustered.tsv')
        tuple path('snps.clustering.king.id'), path('snps.clustering.king')
        path('clustering.png')
    script:
    """
    grep "PSC" ${bcftools_stats} | tail -n+2 > stats.txt

    plink2 --allow-extra-chr \
        --make-king square \
        --out snps.clustering \
        --vcf ${genotype_file}

    python3 $moduleDir/bin/cluster_king.py \
        snps.clustering \
        ${params.samples_file} \
        stats.txt \
        --min-hets ${params.min_hets} \
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