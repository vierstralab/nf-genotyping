#!/usr/bin/env nextflow

params.samples_file='/net/seq/data/projects/regulotyping-h.CD3+/metadata.txt'
params.genome='/net/seq/data/genomes/human/GRCh38/noalts/GRCh38_no_alts'
params.dbsnp_file='/home/jvierstra/data/dbSNP/v151.hg38/All_20180418.fixed-chrom.vcf.gz'
params.genome_ancestral_fasta_file='/home/jvierstra/data/genomes/hg38/ancestral/homo_sapiens_ancestor_GRCh38/homo_sapiens_ancestor.fixed.fa'

//genotyping caller parameters
params.chunksize=5000000 //
params.min_SNPQ=10 // Minimum SNP quality
params.min_GQ=50 // Minimum genotype quality
params.min_DP=12 // Minimum read depth over SNP
params.hwe_cutoff=0.01 // Remove variants that are out of Hardy-Weinberg equilibrium

params.outdir='output'

//DO NOT EDIT BELOW

genome_fasta_file="$params.genome"  + ".fa"
genome_chrom_sizes_file="$params.genome"  + ".chrom_sizes"

// Read samples file
Channel
	.fromPath(params.samples_file)
	.splitCsv(header:true, sep:'\t')
	.map{ row -> tuple( row.donor_id, row.ln_number, row.bamfile ) }
	.groupTuple(by:0)
	.map{ it -> tuple(it[0], it[2].join(" ")) }
	.set{ SAMPLES_AGGREGATIONS_MERGE }

// Merge BAM files by inidividual
process merge_bamfiles {
	tag "${indiv_id}"

	cpus 2

	publishDir params.outdir + '/merged', mode: 'symlink'

	input:
	set val(indiv_id), val(bam_files) from SAMPLES_AGGREGATIONS_MERGE

	output:
	set val(indiv_id), file("${indiv_id}.bam"), file("${indiv_id}.bam.bai") into INDIVS_MERGED_FILES, INDIVS_MERGED_LIST

	script:
	"""
	samtools merge -f -@${task.cpus} ${indiv_id}.bam ${bam_files}
	samtools index ${indiv_id}.bam
	"""
}


// Create sample map file which is used by bcftools to specific input BAM files
INDIVS_MERGED_LIST
	.map{ [it[0], file(it[1]).name].join("\t") }
	.collectFile(
		name: 'sample_indiv_map.tsv',
		newLine: true
	)
	.first()
	.set{ INDIVS_MERGED_SAMPLE_MAP_FILE }

// Channel containing all files needed by bcftools to call genotypes (both VCF and its corresponding index)
INDIVS_MERGED_FILES
	.flatMap{ [file(it[1]), file(it[2])] }
	.set{ INDIVS_MERGED_ALL_FILES }

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
	cat ${chrom_sizes} \
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
	tag "region: ${region}"
	
	scratch true
	cpus 2

	input:
	file genone_fasta from file(genome_fasta_file)
	
	file dbsnp_file from file(params.dbsnp_file)
	file dbsnp_index_file from file("${params.dbsnp_file}.tbi")
 	
	val min_DP from params.min_DP
	val min_SNPQ from params.min_SNPQ
	val min_GQ from params.min_GQ
	val hwe_cutoff from params.hwe_cutoff

	val region from GENOME_CHUNKS_MAP
	file '*' from INDIVS_MERGED_ALL_FILES.collect()
	file 'sample_indiv_map.tsv' from INDIVS_MERGED_SAMPLE_MAP_FILE 

	output:
	file '*filtered.annotated.vcf.gz*' into GENOME_CHUNKS_VCF

	script:
	"""
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
		-i"INFO/DP>=${min_DP}" \
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
		-i"QUAL>=${min_SNPQ} & FORMAT/GQ>=${min_GQ} & FORMAT/DP>=${min_DP}" \
		--SnpGap 3 --IndelGap 10 --set-GTs . \
	| bcftools view \
		-i'GT="alt"' --trim-alt-alleles \
	| bcftools annotate -x ^INFO/DP \
	| bcftools +fill-tags -- -t all \
	| bcftools filter --output-type z -e"INFO/HWE<${hwe_cutoff}" \
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
	file 'all.filtered.vcf.gz*'
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
	
	echo "##INFO=<ID=AA,Number=1,Type=String,Description=\"Inferred ancestral allele -- EPO/PECAN alignments\">" > header.txt

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

