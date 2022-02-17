#!/usr/bin/env nextflow

params.samples_file=''
params.genotype_file=''
params.outdir='output'

//heterzygous filtering parameters
params.min_GQ=50 // Minimum genotype quality
params.min_DP=12 // Minimum read depth over SNP
params.min_AD=4 // Minimum reads per alleles

// Read samples file
Channel
	.fromPath(params.samples_file)
	.splitCsv(header:true, sep:'\t')
	.map{ row -> tuple( row.indiv_id, row.cell_type, row.hotspots_file ) }
	.set{ INDIV_CELL_TYPE }

process filter_variants {
	tag "${indiv_id}:${cell_type}"
	module ' bedops/2.4.40-typical'

	publishDir params.outdir, mode: 'symlink'

	input:
	set val(indiv_id), val(cell_type), val(hotspots_file) from INDIV_CELL_TYPE
	
	file genotype_file from file(params.genotype_file)
	file '*' from file("${params.genotype_file}.csi")
	
	val min_DP from params.min_DP
	val min_AD from params.min_AD
	val min_GQ from params.min_GQ

	output:
	file("${indiv_id}_${cell_type}.bed.gz")
	file("${indiv_id}_${cell_type}.bed.gz.tbi")

	script:
	"""
	bcftools query \
		-s ${indiv_id} \
		-i'GT="alt"' \
		-f'%CHROM\\t%POS0\\t%POS\\t%CHROM:%POS:%REF:%ALT\\t%REF\\t%ALT\\t[%GT\\t%GQ\\t%DP\\t%AD{0}\\t%AD{1}]\\n' \
		${genotype_file} \
	| awk -v OFS="\\t" \
		-v min_GQ=${min_GQ} -v min_AD=${min_AD} -v min_DP=${min_DP}\
		'\$8<min_GQ { next; } \$9<min_DP { next; }\
			(\$7=="0/1" || \$7=="1/0" || \$7=="0|1" || \$7=="1|0") && (\$10<min_AD || \$11<min_AD) { \
				next; \
			} \
			{ print; }' \
	| sort-bed - \
	| grep -v chrX | grep -v chrY | grep -v chrM | grep -v _random | grep -v _alt | grep -v chrUn \
	| bedops -e 1 - ${hotspots_file} \
	| bgzip -c > ${indiv_id}_${cell_type}.bed.gz

	tabix -p bed ${indiv_id}_${cell_type}.bed.gz
	"""
}

