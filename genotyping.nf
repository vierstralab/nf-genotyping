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
		ln -s ${bam_files}.crai ${name}.crai
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
	conda "${params.conda}"
	publishDir "${params.outdir}/region_genotypes"
	cpus 2

	input:
	    each region 
		val indiv_bams
		val n_indivs

	output:
		tuple path("${region}.filtered.annotated.vcf.gz"), path("${region}.filtered.annotated.vcf.gz.csi")

	script:
	indiv_bam_paths = indiv_bams.join('\n')
	"""
	# Workaround
	export OMP_NUM_THREADS=1
	export USE_SIMPLE_THREADED_LEVEL3=1
	echo "${indiv_bam_paths}" | grep -v "crai" > filelist.txt
	cat filelist.txt | xargs -I file basename file | cut -d"." -f1 > samples.txt

	bcftools mpileup \
		--regions ${region} \
		--fasta-ref ${params.genome_fasta_file} \
		--redo-BAQ \
		--adjust-MQ 50 \
		--gap-frac 0.05 \
		--max-depth ${n_indivs * 10} \
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
	vcfs = region_vcfs.join('\n')
	"""
	# Concatenate files
	echo "${vcfs}" > files.txt
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
			.flatMap(it -> tuple(it[1], it[2]))
			.collectFile(sort: true, newLine: true)
			.toList()
		merged_bamfiles.view()
		genome_chunks = create_genome_chunks()
			.flatMap(n -> n.split()).take(5)
		region_genotypes = call_genotypes(genome_chunks, merged_bamfiles, merged_bamfiles.size() / 2)
		genotypes_paths = region_genotypes.map(p -> p[0])
			.collectFile(newLine: true).toList()
		merge_vcfs(genotypes_paths)
	emit:
		merge_vcfs.out
}


workflow {
	SAMPLES_AGGREGATIONS_MERGE = Channel
		.fromPath(params.samples_file)
		.splitCsv(header:true, sep:'\t')
		.map(row -> tuple( row.indiv_id, row.bam_file ))
		.groupTuple().take(5)
	genotyping(SAMPLES_AGGREGATIONS_MERGE)

}