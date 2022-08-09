#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

genome_fasta_file="${params.genome}"  + ".fa"
genome_chrom_sizes_file="${params.genome}"  + ".chrom_sizes"

// Merge BAM files by inidividual

process merge_bamfiles {
	tag "${indiv_id}"

	cpus 2

	publishDir params.outdir + '/merged', mode: 'symlink'

	input:
	tuple val(indiv_id), val(bam_files)

	output:
	tuple val(indiv_id), path("${indiv_id}.bam"), path("${indiv_id}.bam.bai")
	script:
	"""
	samtools merge -f -@${task.cpus} --reference ${genome_fasta_file} ${indiv_id}.bam ${bam_files}
	samtools index ${indiv_id}.bam
	"""
}

// Chunk genome up
process create_genome_chunks {
	executor 'local'

	output:
		stdout

	script:
	"""
	cat ${genome_chrom_sizes_file} \
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
	tag "Region: ${region}"
	
	scratch true
	cpus 2

	input:
		tuple val(region), val(indiv_ids), val(indiv_bams), val(indiv_vcf_indices), val(n_indivs)

	output:
		tuple path("${region}.filtered.annotated.vcf.gz"), path("${region}.filtered.annotated.vcf.gz.csi")

	script:
	"""
	# Workaround
	export OMP_NUM_THREADS=1
	export USE_SIMPLE_THREADED_LEVEL3=1

	echo "${indiv_ids}" > samples.txt
	echo "${indiv_bams}" > filelist.txt

	bcftools mpileup \
		--regions ${region} \
		--fasta-ref ${genome_fasta_file} \
		--redo-BAQ \
		--adjust-MQ 50 \
		--gap-frac 0.05 \
		--max-depth ${n_indivs * 500} \
		--annotate FORMAT/DP,FORMAT/AD \
		--bam-list filelist.txt \
		--output-type u \
	| bcftools call \
		--threads ${task.cpus} \
		--keep-alts \
		--multiallelic-caller \
		--format-fields GQ \
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
		--fasta-ref ${genome_fasta_file} \
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
	publishDir params.outdir + '/genotypes', mode: 'symlink'

	input:
		val(region_vcfs)

	output:
		tuple path('all.filtered.snps.annotated.vcf.gz'), path('all.filtered.snps.annotated.vcf.gz.csi')
	script:
	// TODO, use region for sorting
	region_vcf_files = region_vcfs.map(t -> t.getName()).join('\n')
	"""
	# Concatenate files
	echo ${region_vcfs} > files.txt
	cat files.txt | tr ":-" "\\t" | tr "." "\\t" | cut -f1-3 | paste - files.txt | sort-bed - | awk '{ print \$NF; }' > mergelist.txt

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
		merged_bamfiles = merge_bamfiles(
			samples_aggregations
				.map(it -> tuple(it[0], it[1].join(' ')))
		).toList()
		n_indivs = merged_bamfiles.size()
		all_merged_files = merged_bamfiles.transpose()
			.map(it -> it.join('\n'))
			.toList()
		all_merged_files.view()
		genome_chunks = create_genome_chunks()
			.flatMap( it ->  it.split() )
			.combine(all_merged_files).combine(n_indivs)
		region_genotypes = call_genotypes(genome_chunks)
		merge_vcfs(region_genotypes.map(it -> it[0]).collect())
	emit:
		merge_vcfs.out
}


workflow {
	SAMPLES_AGGREGATIONS_MERGE = Channel
		.fromPath(params.samples_file)
		.splitCsv(header:true, sep:'\t')
		.map(row -> tuple( row.indiv_id, row.bam_file ))
		.groupTuple()
		.map( it -> tuple(it[0], it[1]))
	genotyping(SAMPLES_AGGREGATIONS_MERGE)

}