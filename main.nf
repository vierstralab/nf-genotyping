#!/usr/bin/env nextflow

params.outdir='output'
params.chromsizes_filepath='/net/seq/data/genomes/human/GRCh38/noalts/GRCh38_no_alts.chrom_sizes'
params.fasta_reference_filepath='/net/seq/data/genomes/human/GRCh38/noalts/GRCh38_no_alts.fa'
params.dbsnp_filepath='/home/jvierstra/data/dbSNP/v151.hg38/All_20180418.fixed-chrom.vcf.gz'
params.fasta_ancestral_filepath='/home/jvierstra/data/genomes/hg38/ancestral/homo_sapiens_ancestor_GRCh38/homo_sapiens_ancestor.fixed.fa'

//genotyping caller parameters
params.min_SNPQ=10 // Minimum SNP quality
params.min_GQ=50 // Minimum genotype quality
params.min_DP=12 // Minimum read depth over SNP
params.hwe_cutoff=0.01 // Remove variants that are out of Hardy-Weinbery equilibrium

Channel
	.fromPath("/net/seq/data/projects/regulotyping-h.CD3+/metadata.txt")
	.splitCsv(header:true, sep:'\t')
	.map{ row -> tuple( row.donor_id, row.ds_number, row.ln_number, row.ag_number, file(row.bamfile) ) }
	.set{ SAMPLES_AGGREGATIONS }


// This is a hack to deal with input files that have the same basename
process symlink_input_files {
	tag "${donor_id}:AG${ag_number}"
	
	executor 'local'	

	input:
	set val(donor_id), val(ds_number), val(ln_number), val(ag_number), file(bam_file) from SAMPLES_AGGREGATIONS

	output:
	set val(donor_id), file("*.bam") into SAMPLES_AGGREGATIONS_SYMLINKED

	script:
	"""
	ln -s ${bam_file} "${ag_number}.bam"
	"""
}

SAMPLES_AGGREGATIONS_SYMLINKED
	.groupTuple(by:0)
	.set{ SAMPLES_AGGREGATIONS_MERGE }

process merge_bamfiles {
	tag "${donor_id}"

	cpus 2

	input:
	set val(donor_id), file(bam_files) from SAMPLES_AGGREGATIONS_MERGE

	output:
	file '*.bam*' into DONORS_MERGED

	script:
	"""
	samtools merge -f -@${task.cpus} ${donor_id}.bam ${bam_files}
	samtools index ${donor_id}.bam
	"""
}


process create_genome_chunks {
	executor 'local'

	input:
	file chromsizes from file(params.chromsizes_filepath)

	output:
	stdout into GENOME_CHUNKS

	script:
	"""
	cat ${chromsizes} \
  	| awk -v step=5000000 -v OFS="\t" \
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

process call_genotypes {
	tag "region: ${region}"
	
	scratch true
	cpus 2

	input:
	file fasta_reference from file(params.fasta_reference_filepath)
	
	file dbsnp_file from file(params.dbsnp_filepath)
	file dbsnp_index_file from file("${params.dbsnp_filepath}.tbi")
 	
	val min_DP from params.min_DP
	val min_SNPQ from params.min_SNPQ
	val min_GQ from params.min_GQ
	val hwe_cutoff from params.hwe_cutoff

	val region from GENOME_CHUNKS_MAP
	file bam_files from DONORS_MERGED.collect()

	output:
	file '*filtered.annotated.vcf.gz*' into GENOME_CHUNKS_VCF

	script:
	"""
	ls *.bam > filelist.txt

	bcftools mpileup \
		--regions ${region} \
		--fasta-ref ${fasta_reference} \
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

	bcftools query -l ${region}.vcf.gz | grep -o "[^/]*\$" | sed "s/\\.bam//" > samples.txt

	bcftools reheader \
		-s samples.txt \
	 	${region}.vcf.gz \
	| bcftools norm \
		--threads ${task.cpus} \
		--check-ref w \
		-m - \
		--fasta-ref ${fasta_reference} \
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


process merge_vcfs {

	scratch true
	publishDir params.outdir + '/calls', mode: 'copy'

	input:
	file '*' from GENOME_CHUNKS_VCF.collect()
	file fasta_ancestral from file(params.fasta_ancestral_filepath)

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
	#
	#
	
	echo "##INFO=<ID=AA,Number=1,Type=String,Description=\"Inferred ancestral allele-- EPO/PECAN alignments\">" > header.txt

	# Contigs with ancestral allele information
	faidx -i chromsizes ${fasta_ancestral} | cut -f1 > chroms.txt
	
	# Get SNPs in BED-format; remove contigs with no FASTA sequence
	bcftools query  -f "%CHROM\t%POS0\t%POS\t%REF\t%ALT\n" all.filtered.snps.vcf.gz \
	| grep -w -f chroms.txt \
	> all.filtered.snps.bed

	# Get ancestral allele from FASTA file and make a TABIX file
	faidx -i transposed \
		-b all.filtered.snps.bed \
		${fasta_ancestral} \
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

process generate_h5_tables {

	scratch true
	publishDir params.outdir + '/h5', mode: 'copy'

	input:
		set file(vcf_file), file(vcf_index_file) from FILTERED_SNPS_VCF
		file chromsizes from file(params.chromsizes_filepath)

	output:
		file '*.h5'

	script:
	"""
	
	mkdir by_chrom

	chroms=("\$(tabix -l ${vcf_file})")
	for chrom in \${chroms[@]}; do
		bcftools view -r \${chrom} -Oz ${vcf_file} > by_chrom/\${chrom}.vcf.gz
		bcftools index by_chrom/\${chrom}.vcf.gz
	done

	gzip -c ${chromsizes} > chromsizes.txt.gz

	snp2h5 --chrom chromsizes.txt.gz \
		--format vcf \
		--haplotype haplotypes.h5 \
		--snp_index snp_index.h5 \
		--snp_tab snp_tab.h5 \
		by_chrom/*.vcf.gz

	"""

}
