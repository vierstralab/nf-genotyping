#!/usr/bin/env nextflow
nextflow.enable.dsl = 1

genome_fasta_file="${params.genome}"  + ".fa"
genome_chrom_sizes_file="${params.genome}"  + ".chrom_sizes"

// Read samples file
Channel
	.fromPath(params.samples_file)
	.splitCsv(header:true, sep:'\t')
	.map(row -> tuple( row.indiv_id, row.bam_file ))
	.groupTuple()
	.map( it -> tuple(it[0], it[1].join(' ')))
	.set{ SAMPLES_AGGREGATIONS_MERGE }
// Merge BAM files by inidividual

process merge_bamfiles {
	tag "${indiv_id}"

	cpus 2

	publishDir params.outdir + '/merged', mode: 'symlink'

	input:
	set val(indiv_id), val(bam_files) from SAMPLES_AGGREGATIONS_MERGE

	output:
	set val(indiv_id), file("${indiv_id}.bam"), file("${indiv_id}.bam.bai") into INDIV_MERGED_FILES, INDIV_MERGED_LIST

	script:
	"""
	samtools merge -f -@${task.cpus} --reference ${genome_fasta_file} ${indiv_id}.bam ${bam_files}
	samtools index ${indiv_id}.bam
	"""
}


// Create sample map file which is used by bcftools to specific input BAM files
INDIV_MERGED_LIST
	.map{ [it[0], file(it[1]).name].join("\t") }
	.collectFile(
		name: 'sample_indiv_map.tsv',
		newLine: true
	)
	.first()
	.set{ INDIV_MERGED_SAMPLE_MAP_FILE }

// Channel containing all files needed by bcftools to call genotypes (both VCF and its corresponding index)
INDIV_MERGED_FILES
	.flatMap{ [file(it[1]), file(it[2])] }
	.set{ INDIV_MERGED_ALL_FILES }

// Chunk genome up
process create_genome_chunks {
	executor 'local'

	input:
	file chrom_sizes from file(genome_chrom_sizes_file)
	val chunksize from params.chunksize

	output:
	stdout into GENOME_CHUNKS

	script:
	"""
	cat ${chrom_sizes} | grep -i '^chr\([1-9]\|1[0-9]\|2[0-2]\|[XY]\)\s' \
  	| awk -v step=${chunksize} -v OFS="\t" \
		'{ \
			for(i=step; i<=\$2; i+=step) { \
				print \$1":"i-step+1"-"i; \
			} \
			print \$1":"i-step+1"-"\$2; \
		}' \
	"""
} 

GENOME_CHUNKS
	.flatMap{ it.split() }
	.set{ GENOME_CHUNKS_MAP }

// Call genotypes per region
process call_genotypes {
	tag "Region: ${region}"
	
	scratch true
	cpus 2

	input:
	file genome_fasta from file(genome_fasta_file)
	file dbsnp_file from file(params.dbsnp_file)
	file dbsnp_index_file from file("${params.dbsnp_file}.tbi")
	val region from GENOME_CHUNKS_MAP
	file '*' from INDIV_MERGED_ALL_FILES.collect()
	file 'sample_indiv_map.tsv' from INDIV_MERGED_SAMPLE_MAP_FILE 

	output:
	file '*filtered.annotated.vcf.gz*' into GENOME_CHUNKS_VCF

	script:
	"""
	# Workaround
	export OMP_NUM_THREADS=1
	export USE_SIMPLE_THREADED_LEVEL3=1


	cut -f1 sample_indiv_map.tsv > samples.txt
	cut -f2 sample_indiv_map.tsv > filelist.txt

	bcftools mpileup \
		--regions ${region} \
		--fasta-ref ${genome_fasta} \
		--redo-BAQ \
		--adjust-MQ 50 \
		--gap-frac 0.05 \
		--max-depth 10000 --max-idepth 200000 \
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
		--fasta-ref ${genome_fasta} \
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
		-a ${dbsnp_file} \
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
	file '*' from GENOME_CHUNKS_VCF.collect()
	file genome_fasta_ancestral from file(params.genome_ancestral_fasta_file)

	output:
	set file('all.filtered.snps.annotated.vcf.gz'), file('all.filtered.snps.annotated.vcf.gz.csi') into FILTERED_SNPS_VCF
	script:
	"""
	# Concatenate files
	ls *.filtered.annotated.vcf.gz > files.txt
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
	faidx -i chromsizes ${genome_fasta_ancestral} | cut -f1 > chroms.txt
	
	# Get SNPs in BED-format; remove contigs with no FASTA sequence
	bcftools query  -f "%CHROM\t%POS0\t%POS\t%REF\t%ALT\n" all.filtered.snps.vcf.gz \
	| grep -w -f chroms.txt \
	> all.filtered.snps.bed

	# Get ancestral allele from FASTA file and make a TABIX file
	faidx -i transposed \
		-b all.filtered.snps.bed \
		${genome_fasta_ancestral} \
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
