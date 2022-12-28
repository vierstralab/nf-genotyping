#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

params.conda = "${moduleDir}/environment.yml"

process merge_bamfiles {
	tag "${indiv_id}:${bam_files_names.size()}"

	conda "${params.conda}"
	cpus 2
	memory 32.GB

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
	memory 500.MB
	conda params.conda
	scratch true

	output:
		stdout

	script:
	"""
	cat ${params.genome_chrom_sizes} \
	| grep -v chrX | grep -v chrY | grep -v chrM | grep -v _random | grep -v _alt | grep -v chrUn \
  	| awk -v step=${params.chunksize} -v OFS="\t" \
		'{ \
			for(i=step; i<=\$2; i+=step) { \
				print \$1"\t"i-step+1"\t"i; \
			} \
			print \$1"\t"i-step+1"\t"\$2; \
		}' > genome_chunks.bed
	bedtools subtract -a genome_chunks.bed -b ${params.encode_blacklist_regions} | awk '{ print \$1":"\$2"-"\$3 }'
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
		--max-depth 10000 \
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
	conda "${params.conda}"

	input:
		val region_vcfs

	output:
		tuple path(name), path("${name}.csi")

	script:
	vcfs = region_vcfs.join('\n')
	name = "all.filtered.snps.vcf.gz"
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
	> ${name}
	
	bcftools index ${name}
	"""
}


process annotate_vcf {
	conda "${params.conda}"
	publishDir params.outdir + '/genotypes'

	input:
		tuple path(snps_vcf), path(snps_vcf_index)

	output:
		tuple path(name), path("${name}.csi")

	script:
	name = 'all.filtered.snps.annotated.vcf.gz'
	"""
	echo '##INFO=<ID=AA,Number=1,Type=String,Description="Inferred ancestral allele -- EPO/PECAN alignments">' > header.txt
	echo '##INFO=<ID=AAF,Number=1,Type=Float,Description="TOPMED alternative allele frequency">' >> header.txt
	echo '##INFO=<ID=RAF,Number=1,Type=Float,Description="TOPMED reference allele frequency">' >> header.txt
	
	# Get SNPs in BED-format
	bcftools query -f "%CHROM\t%POS0\t%POS\t%REF\t%ALT\n" \
		${snps_vcf} > all.filtered.snps.bed

	# Get ancestral allele from FASTA file and make a TABIX file
	faidx -i chromsizes ${genome_fasta_ancestral} | cut -f1 > chroms.txt

	cat all.filtered.snps.bed | grep -w -f chroms.txt > ancestral_chrs.snps.bed
	faidx -i transposed \
		-b ancestral_chrs.snps.bed \
		${params.genome_ancestral_fasta_file} \
		| paste - ancestral_chrs.snps.bed	\
		| awk -v OFS="\t" '{ print \$5, \$6, \$7, \$8, \$9, \$4; }' \
		| bgzip -c > ancestral.allele.annotation.bed.gz
	
	bcftools query -f "%CHROM\t%POS0\t%POS\t%REF\t%ALT\t%INFO/TOPMED\n" ${params.dbsnp_file} | \
		bedtools intersect -a stdin -b all.filtered.snps.bed -sorted -wa \
		| bgzip -c > dbsnp_annotations.bed.gz

	# Add TOPMED freqs annotation
	python3 $moduleDir/bin/explode_topmed_annotations.py \
		ancestral.allele.annotation.bed.gz dbsnp_annotations.bed.gz all.filtered.snps.bed  | \
		bgzip -c > all.filtered.snps.annotations.bed.gz

	tabix all.filtered.snps.annotations.bed.gz

	# Annotate VCF
	bcftools annotate \
		--output-type z \
		-h header.txt \
		-a all.filtered.snps.annotations.bed.gz \
		-c CHROM,BEG,END,-,ALT,INFO/AAF,INFO/RAF,INFO/AA \
		${snps_vcf} \
	> ${name}
	
	bcftools index ${name}
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
			.toSortedList( { a, b -> b <=> a } )
		genome_chunks = create_genome_chunks()
			.flatMap(n -> n.split())
		region_genotypes = call_genotypes(genome_chunks, merged_bamfiles)
		genotypes_paths = region_genotypes.map(p -> p[0])
			.toSortedList( { a, b -> b <=> a } )
		out = merge_vcfs(genotypes_paths) | annotate_vcf
	emit:
		out
}


workflow {
	SAMPLES_AGGREGATIONS_MERGE = Channel
		.fromPath(params.samples_file)
		.splitCsv(header:true, sep:'\t')
		.map(row -> tuple( row.indiv_id, row.bam_file ))
		.groupTuple()
	genotyping(SAMPLES_AGGREGATIONS_MERGE)

}