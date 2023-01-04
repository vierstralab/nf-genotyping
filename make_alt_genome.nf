#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

params.conda = "${moduleDir}/environment.yml"
params.genotype_file = "${launchDir}/${params.outdir}/genotypes/all.filtered.snps.annotated.vcf.gz"


process make_iupac_genome {
	conda "${params.conda}"
    tag "${sample_id}"
	publishDir "${params.outdir}/alt_genome"

	input:
		val sample_id

	output:
		tuple path("${name}"), path("${name}.fai")

	script:
    prefix = sample_id ?: "all_samples"
	name = "${prefix}.iupac.genome.fa"
	additional_params = sample_id ? "--sample ${sample_id}" : ""
    """
    python3 $moduleDir/bin/nonref_genome.py ${params.genome_fasta_file} ${name} ${params.genotype_file} ${additional_params}
    """
}

// Make iupac coded genome from genotype_file
workflow {
    params.sample_id = ""
	make_iupac_genome(params.sample_id)
}
