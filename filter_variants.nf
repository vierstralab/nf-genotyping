#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

params.conda = "$moduleDir/environment.yml"

// Might add some filtering for different ag_ids
process filter_variants {
	tag "${outname}"
	conda params.conda
	publishDir "${params.outdir}/bed_files"

	input:
		tuple val(indiv_id), val(ag_id)

	output:
		tuple val(indiv_id), val(ag_id), path(outname), path("${outname}.tbi")

	script:
	outname = "${indiv_id}_${ag_id}.bed.gz"
	"""
	bcftools query \
		-s ${indiv_id} \
		-i'GT="alt"' \
		-f'%CHROM\\t%POS0\\t%POS\\t%ID\\t%REF\\t%ALT\\t%INFO/MAF\\t[%GT\\t%GQ\\t%DP\\t%AD{0}\\t%AD{1}]\\n' \
		${params.genotype_file} \
	| awk -v OFS="\\t" \
		-v min_GQ=${params.min_GQ} -v min_AD=${params.min_AD} -v min_DP=${params.min_DP}\
		'\$9<min_GQ { next; } \$10<min_DP { next; }\
			(\$8=="0/1" || \$8=="1/0" || \$8=="0|1" || \$8=="1|0") && (\$11<min_AD || \$12<min_AD) { \
				next; \
			} \
			{ print; }' \
	| sort-bed - \
	| grep -v chrX | grep -v chrY | grep -v chrM | grep -v _random | grep -v _alt | grep -v chrUn \
	| bgzip -c > ${outname}
	# Check if file is empty
	if [ -s  ${outname} ]; then
		tabix -f -p bed ${outname}
	else
		touch "${outname}.tbi"
	fi
	"""
}

// process extend_metadata {
// 	publishDir params.outdir

// 	input:
// 		tuple val(indiv_id), val(ag_id), path(bed_files)

// 	output:
// 		path name

// 	script:
// 	name = 'metadata.with_intervals.txt'
// 	column = bed_files.join('\n')
// 	"""
// 	echo "filtered_sites_file\n${column}" > columns.txt
// 	paste ${params.samples_file} columns.txt > ${name}
// 	"""
// }

workflow filterVariants {
	take:
		indiv_cell_types_meta
	main:
		variants_paths = filter_variants(indiv_cell_types_meta)
		//extend_metadata(variants_paths.map(it -> it[2]).collect())
	emit:
		filter_variants.out
		//extend_metadata.out
}

workflow {
	// Read samples file
	INDIV_CELL_TYPE = Channel
		.fromPath(params.samples_file)
		.splitCsv(header:true, sep:'\t')
		.map(row -> tuple(row.indiv_id, row.ag_id))
	filterVariants(INDIV_CELL_TYPE)
}