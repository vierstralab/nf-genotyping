#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

params.conda = "${moduleDir}/environment.yml"

process merge_bamfiles {
	tag "${indiv_id}"

	conda "${params.conda}"
	cpus 2

	input:
		tuple val(indiv_id), val(bam_files)

	output:
		tuple val(indiv_id), path(name), path("${name}.*ai")

	script:
	bam_files_names = bam_files.split(' ')
	if ( bam_files_names.size() > 1 ) {
		name = "${indiv_id}.cram"
		"""
		samtools merge -f -@${task.cpus} --reference ${params.genome_fasta_file} ${name} ${bam_files}
		samtools index ${name}
		"""
	} else {
		bam_ext = file(bam_files).extension
		name = "${indiv_id}.${bam_ext}"
		"""
		ln -s ${bam_files} ${name}
		samtools index ${name}
		"""
	}
}

// Chunk genome up
process create_genome_chunks {
	executor 'local'

	output:
		stdout

	script:
	"""
	cat ${params.genome_chrom_sizes} \
	| grep -v chrX | grep -v chrY | grep -v chrM | grep -v _random | grep -v _alt | grep -v chrUn \
  	| awk -v step=${params.chunksize} -v OFS="\t" \
		'{ \
			for(i=step; i<=\$2; i+=step) { \
				print \$1":"i-step+1"-"i; \
			} \
			print \$1":"i-step+1"-"\$2; \
		}' \
	"""
} 

// Call genotypes per region
process call_genotypes {
	tag "${region}"
	scratch true
	cache true
	conda "${params.conda}"
	publishDir "${params.outdir}/region_genotypes"
	cpus 2

	input:
	    each region 
		path indiv_bams
		val n_indivs

	output:
		tuple path("${region}.filtered.annotated.vcf.gz"), path("${region}.filtered.annotated.vcf.gz.csi")

	script:
	indiv_bams_names = indiv_bams.tap { it.name }.findAll({ !it.endsWith('ai')})
	indiv_ids = indiv_bams_names.tap{ it.simpleName }
	"""
	# Workaround
	export OMP_NUM_THREADS=1
	export USE_SIMPLE_THREADED_LEVEL3=1

	echo "${indiv_ids}" | tr " " "\n" > samples.txt
	echo "${indiv_bams_names}" | tr " " "\n" > filelist.txt

	bcftools mpileup \
		--regions ${region} \
		--fasta-ref ${params.genome_fasta_file} \
		--redo-BAQ \
		--adjust-MQ 50 \
		--gap-frac 0.05 \
		--max-depth ${n_indivs * 100} \
		--annotate FORMAT/DP,FORMAT/AD \
		--bam-list filelist.txt \
		--output-type u \
	| bcftools call \
		--threads ${task.cpus} \
		--keep-alts \
		--multiallelic-caller \
		-a GQ \
		--output-type v \
	| bcftools filter \
		-i"INFO/DP>=${params.min_DP}" \
		--output-type z - \
	> ${region}.vcf.gz

	bcftools index ${region}.vcf.gz

	bcftools reheader \
		-s samples.txt \
	 	${region}.vcf.gz \
	| bcftools norm \
		--threads ${task.cpus} \
		--check-ref w \
		-m - \
		--fasta-ref ${params.genome_fasta_file} \
	| bcftools filter \
		-i"QUAL>=${params.min_SNPQ} & FORMAT/GQ>=${params.min_GQ} & FORMAT/DP>=${params.min_DP}" \
		--SnpGap 3 --IndelGap 10 --set-GTs . \
	| bcftools view \
		-i'GT="alt"' --trim-alt-alleles \
	| bcftools annotate -x ^INFO/DP \
	| bcftools +fill-tags \
	| bcftools filter --output-type z -e"INFO/HWE<${params.hwe_cutoff}" \
	> ${region}.filtered.vcf.gz

	bcftools index ${region}.filtered.vcf.gz

	bcftools annotate \
		-r ${region} \
		-a ${params.dbsnp_file} \
		--columns ID,CAF,TOPMED \
		--output-type z \
  		${region}.filtered.vcf.gz \
	> ${region}.filtered.annotated.vcf.gz

	bcftools index ${region}.filtered.annotated.vcf.gz		
	"""
}

// Merge VCF chunks and add ancestral allele information
process merge_vcfs {
	scratch true
	conda "${params.conda}"
	publishDir params.outdir + '/genotypes'

	input:
		val region_vcfs

	output:
		tuple path('all.filtered.snps.annotated.vcf.gz'), path('all.filtered.snps.annotated.vcf.gz.csi')

	script:
	region_vcf_names = region_vcfs.join('\n')
	"""
	# Concatenate files
	echo "${region_vcf_names}" > files.txt
	cat files.txt | sed 's!.*/!!' | tr ":-" "\\t" | tr "." "\\t" | cut -f1-3 | paste - files.txt | sort-bed - | awk '{ print \$NF; }' > mergelist.txt

	bcftools concat \
		--output-type z \
		-f mergelist.txt \
	> all.filtered.vcf.gz
	
	bcftools index all.filtered.vcf.gz
	
	# Output only SNPs
	bcftools view \
		-m2 -M2 -v snps \
		--output-type z \
		all.filtered.vcf.gz \
	> all.filtered.snps.vcf.gz
	
	bcftools index all.filtered.snps.vcf.gz	

	# Annotate ancestral allele
	
	echo '##INFO=<ID=AA,Number=1,Type=String,Description="Inferred ancestral allele -- EPO/PECAN alignments">' > header.txt

	# Contigs with ancestral allele information
	faidx -i chromsizes ${params.genome_ancestral_fasta_file} | cut -f1 > chroms.txt
	
	# Get SNPs in BED-format; remove contigs with no FASTA sequence
	bcftools query  -f "%CHROM\t%POS0\t%POS\t%REF\t%ALT\n" all.filtered.snps.vcf.gz \
	| grep -w -f chroms.txt \
	> all.filtered.snps.bed

	# Get ancestral allele from FASTA file and make a TABIX file
	faidx -i transposed \
		-b all.filtered.snps.bed \
		${params.genome_ancestral_fasta_file} \
	| paste - all.filtered.snps.bed	\
	| awk -v OFS="\t" '{ print \$5, \$7, \$4; }' \
	| bgzip -c > all.filtered.snps.ancestral.tab.gz

	tabix -b 2 -e -2 all.filtered.snps.ancestral.tab.gz

	# Annotate VCF
	bcftools annotate \
		--output-type z \
		-h header.txt \
		-a all.filtered.snps.ancestral.tab.gz \
		-c CHROM,POS,INFO/AA \
		all.filtered.snps.vcf.gz \
	> all.filtered.snps.annotated.vcf.gz
	
	bcftools index all.filtered.snps.annotated.vcf.gz
	"""
}

workflow genotyping {
	take: 
		samples_aggregations
	main:
		bam_files = samples_aggregations
			.map(it -> tuple(it[0], it[1].join(' ')))
		merged_bamfiles = merge_bamfiles(bam_files)
			.flatMap(i -> tuple(i[1], i[2]))
			.collect()
		n_indivs = merged_bamfiles.size()
		// Workaround. Collect uses flatMap, which won't work here
		genome_chunks = create_genome_chunks()
			.flatMap(n -> n.split())
		region_genotypes = call_genotypes(genome_chunks, merged_bamfiles, n_indivs)
		merge_vcfs(region_genotypes.map(p -> p[0]).collect())
	emit:
		merge_vcfs.out
}


workflow {
	SAMPLES_AGGREGATIONS_MERGE = Channel
		.fromPath(params.samples_file)
		.splitCsv(header:true, sep:'\t')
		.map(row -> tuple( row.indiv_id, row.bam_file ))
		.groupTuple()
	genotyping(SAMPLES_AGGREGATIONS_MERGE)

}