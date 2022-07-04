process clusterIndivs {

    publishDir "${params.outdir}/clustering"

    input:
        tuple path(vcf_file), path("${vcf_file}.csi") from FILTERED_SNPS_VCF
    output:
        path "./*"
    script:
    """
    plink2 --allow-extra-chr \
    --make-king square \
    --out snps.clustering \
    --vcf ${vcf_file}

    python3 $baseDir/bin/cluster_king.py --matrix snps.clustering.king \
    --matrix-ids snps.clustering.king.id \
    --meta-file ${params.samples_file}
    --outpath ./
    """
}