#!/usr/bin/env nextflow

params.vcf_file = "/home/sabramov/projects/ENCODE4/genotyping/outdir/genotypes/all.filtered.snps.annotated.vcf.gz"

process clusterIndivs {
    publishDir "${outdir}/clustering"

    input:
        tuple path(vcf_file), path("${vcf_file}.csi")
    output:
        path "./*"
    script:
    """
    plink2 --allow-extra-chr \
    --make-king square \
    --out snps.clustering \
    --vcf ${vcf_file}

    python3 $basePath/bin/cluster_king.py --matrix snps.clustering.king \
    --matrix-ids snps.clustering.king.id \
    --meta-file ${params.samples_file}
    --outpath ./
    """
}

workflow {
    vcf_file = Channel.of([params.vcf_file, "${params.vcf_file}.csi"])
    clusterIndivs(vcf_file)
}