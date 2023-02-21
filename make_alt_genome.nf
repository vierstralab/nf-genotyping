#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

params.conda = "${moduleDir}/environment.yml"
params.genotype_file = "${launchDir}/${params.outdir}/genotypes/all.filtered.snps.annotated.vcf.gz"


process make_iupac_genome {
	conda "${params.conda}"
    tag "${prefix}"
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
    python3 $moduleDir/bin/nonref_genome.py ${params.genome_fasta_file} ${params.genotype_file} ${name} ${additional_params}
    """
}

process make_dhs_annotation {
	conda "${params.conda}"
    tag "${ag_id}"

	input:
		tuple val(indiv_id), val(ag_id), path(hotspots_file)

	output:
		tuple path("${name}"), path("${name}.fai")

	script:
    name = "${ag_id}.dhs_annotation.bed"
    """
	cat ${params.genotype_annotation} | awk '(\$6 == "${indiv_id}") { print; }' \
		| bedtools intersect -a ${hotspots_file} -b stdin -wa -wb > ${name}
    """
}

workflow annotateDHS {
	// TODO add a step to the pipeline
	params.genotype_annotation = "${launchDir}/${params.outdir}/genotypes/genotypes_by_indiv.bed"
	Channel.fromPath(params.samples_file)
		| splitCsv(header:true, sep:'\t')
		| map(row -> tuple(row.indiv_id, row.ag_id, file(row.hotspots_file)))
		| make_dhs_annotation
		| collectFile(storeDir: "${params.outdir}",
			name: "genotypes_annotation.bed")

}
// Make iupac coded genome from genotype_file
workflow {
    params.sample_id = ""
	make_iupac_genome(params.sample_id)
}
