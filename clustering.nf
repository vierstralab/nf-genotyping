#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

params.conda = "$moduleDir/environment.yml"
params.genotype_file = "${launchDir}/${params.outdir}/genotypes/all.filtered.snps.vcf.gz"


process cluster_indivs {

    publishDir "${params.outdir}/clustering"
    conda params.conda

    input:
        path genotype_file
    output:
        tuple path('metadata.clustered.tsv'), path('clustering.png')
        tuple path('snps.clustering.king.id'), path('snps.clustering.king')
    script:
    """
    plink2 --allow-extra-chr \
    --make-king square \
    --out snps.clustering \
    --vcf ${genotype_file}

    python3 $moduleDir/bin/cluster_king.py --matrix snps.clustering.king \
    --matrix-ids snps.clustering.king.id \
    --meta-file ${params.samples_file} \
    --outpath ./
    """
}

workflow clusterIndivs {
    take:
        genotype_file
    main:
        updated_meta = cluster_indivs(genotype_file)[0]
    emit:
        updated_meta
}

// Path to resulting genotype file a.k.a the output of genotyping.nf script (all.filtered.snps.annotated.vcf.gz)
workflow {
    clusterIndivs(params.genotype_file)
}