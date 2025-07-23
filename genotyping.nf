#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

params.conda = "${moduleDir}/environment.yml"


def set_key_for_group_tuple(ch) {
    ch | groupTuple()
        | map(it -> tuple(groupKey(it[0], it[1].size()), *it[1..(it.size()-1)]))
        | transpose()
}


process merge_bamfiles {
	tag "${indiv_id}:${s}"

	conda "${params.conda}"
	label "medmem"
    cpus 2

	input:
		tuple val(indiv_id), path(bam_files, stageAs: "?/*"), path(bam_files_index, stageAs: "?/*")

	output:
		tuple path(name), path("${name}.crai")

	script:
	s = indiv_id.size
	name = "${indiv_id}.merged.cram"
	if ( s > 1 ) {
		"""
		samtools merge -f -@${task.cpus} --reference ${params.genome_fasta_file} ${name} ${bam_files}
		samtools index ${name}
		"""
	} else {
        if (bam_files.extension == "bam") {
            """
            samtools view -@${task.cpus} \
                -C \
                --reference ${params.genome_fasta_file} \
                -o ${name} \
                ${bam_files}
            samtools index ${name}
            """
        } else {
            """
            ln -s ${bam_files} ${name}
            ln -s ${bam_files}.crai ${name}.crai
            """
        }
	}
}

// Chunk genome up
process create_genome_chunks {
	memory 500.MB
	conda "${params.conda}"
	scratch true

	output:
		stdout

	script:
	"""
	cat ${params.genome_chrom_sizes} \
	    | grep -v chrX \
        | grep -v chrY \
        | grep -v chrM \
        | grep -v _random \
        | grep -v _alt \
        | grep -v chrUn \
        | awk -v step=${params.chunksize} -v OFS="\t" \
            '{ \
                for(i=step; i<=\$2; i+=step) { \
                    print \$1"\t"i-step+1"\t"i; \
                } \
                print \$1"\t"i-step+1"\t"\$2; \
            }' > genome_chunks.bed
	bedtools subtract \
        -a genome_chunks.bed \
        -b ${params.encode_blacklist_regions} \
        | awk '{ print \$1":"\$2"-"\$3 }'
	"""
} 

// Call genotypes per region
process call_genotypes {
	tag "${region}"
	conda "${params.conda}"
    label "highmem"
	cpus 2

	input:
	    tuple val(region), path(bamlist)

	output:
		tuple path("${region}.filtered.annotated.vcf.gz"), path("${region}.filtered.annotated.vcf.gz.csi")

	script:
	"""
	# Workaround
	export OMP_NUM_THREADS=1
	export USE_SIMPLE_THREADED_LEVEL3=1

    cat ${bamlist} | grep -v '.crai' > filelist.txt
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
	scratch true
    label "highmem"

	input:
		val region_vcfs

	output:
		tuple path(name), path("${name}.csi"), path("all.filtered.vcf.gz"), path("all.filtered.vcf.gz.csi")

	script:
	name = "all.filtered.snps.vcf.gz"
	region_vcfs_list = region_vcfs.join('\n')
	"""
	# Concatenate files
	echo "${region_vcfs_list}" | grep -v ".csi" > files.txt
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
	publishDir "${params.outdir}/genotypes"
    label "highmem"

	input:
		tuple path("snps.vcf.gz"), path("snps.vcf.gz.csi")

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
		snps.vcf.gz | sort-bed - > all.filtered.snps.bed

	# Get ancestral allele from FASTA file and make a TABIX file
	faidx -i chromsizes ${params.genome_ancestral_fasta_file} | cut -f1 > chroms.txt

	cat all.filtered.snps.bed | grep -w -f chroms.txt | sort-bed - > ancestral_chrs.snps.bed
	faidx -i transposed \
		-b ancestral_chrs.snps.bed \
		${params.genome_ancestral_fasta_file} \
		| paste - ancestral_chrs.snps.bed	\
		| awk -v OFS="\t" '{ print \$5, \$6, \$7, \$8, \$9, \$4; }' \
		| bgzip -c > ancestral.allele.annotation.bed.gz
	
	bcftools query -f "%CHROM\t%POS0\t%POS\t%REF\t%ALT\t%INFO/TOPMED\n" ${params.dbsnp_file} \
		| sort-bed - \
		| bedtools intersect -a stdin -b all.filtered.snps.bed -sorted -wa \
		| uniq \
		| bgzip -c > dbsnp_annotations.bed.gz

	# Add TOPMED freqs annotation
	python3 $moduleDir/bin/explode_topmed_annotations.py \
		ancestral.allele.annotation.bed.gz dbsnp_annotations.bed.gz all.filtered.snps.bed annot.bed
	sort-bed annot.bed \
		| bgzip -c > all.filtered.snps.annotations.bed.gz

	tabix all.filtered.snps.annotations.bed.gz

	# Annotate VCF
	bcftools annotate \
		--output-type z \
		-h header.txt \
		--pair-logic exact \
		-a all.filtered.snps.annotations.bed.gz \
		-c CHROM,BEG,END,-,ALT,INFO/AAF,INFO/RAF,INFO/AA \
		snps.vcf.gz \
	> ${name}
	
	bcftools index ${name}
	"""
}

process vcf_stats {
	conda "${params.conda}"
	publishDir "${params.outdir}", pattern: "stats/*"
	publishDir "${params.outdir}/stats", pattern: "${name}"

	input:
		tuple path(vcf), path(vcf_index)

	output:
		tuple path(name), path("stats/*")

	script:
	prefix = "stats"
	name = "bcftools.stats.txt"
	"""
	bcftools stats -s - ${vcf} > ${name}
	plot-vcfstats -p stats ${name}
	"""
}


workflow genotyping {
	take: 
		bam_files
	main:
		merged_bamfiles = bam_files
            | groupTuple()
            | merge_bamfiles
            | flatten()
			| map(it -> it.toString()) // Get the full path to the file
            | collectFile(name: 'bam_list.txt', newLine: true)

		genome_chunks = create_genome_chunks()
			| flatMap(n -> n.split())
            | combine(merged_bamfiles)

		merged_vcf = call_genotypes(genome_chunks)
			| collect(flat: true, sort: true)
			| merge_vcfs
            | map(it -> tuple(it[0], it[1]))
            | (annotate_vcf & vcf_stats)

	emit:
		annotate_vcf.out
}

workflow {
	Channel.fromPath(params.samples_file)
		| splitCsv(header:true, sep:'\t')
		| map(row -> tuple(
                row.indiv_id,
                file(row.cram_file),
                file(row?.cram_index ?: "${row.cram_file}.crai")
            )
        )
		| filter { !it[0].isEmpty() }
        | set_key_for_group_tuple
		| genotyping
}


workflow annotateVCF {
	Channel.of([file(params.genotype_file), file("${params.genotype_file}.csi")])
		| annotate_vcf
}

// Defunc
workflow annotateVCFStats {
	Channel.of([file(params.genotype_file), file("${params.genotype_file}.csi")])
	 	| vcf_stats
}
