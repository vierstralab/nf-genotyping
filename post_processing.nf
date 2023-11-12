#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

params.conda = "$moduleDir/environment.yml"

process extract_variants_from_vcf {
    params.conda

    output:
        path name

    script:
    name = "unique_variants.bed"
    """
    echo -e "#chr\tstart\tend\tID\tref\talt\tAA" > ${name}
    bcftools query -f'%CHROM\t%POS0\t%POS\t%ID\t%REF\t%ALT\t%INFO/AA\n' \
        ${params.genotype_file} >> ${name}
    """
}


process motif_counts {
    scratch true
    tag "${motif_id}"
    conda params.conda
    publishDir "${params.outdir}/motif_counts"

    input:
        tuple val(motif_id), path(pwm_path), path(moods_file), path(variants)

    output:
        tuple val(motif_id), path(counts_file)

    script:
    counts_file = "${motif_id}.counts.bed"
    """
    
    echo -e "#chr\tstart\tend\tID\tref\talt\tmotif\toffset\twithin\tstrand\tref_score\talt_score\tseq" > ${counts_file}
    zcat ${moods_file} \
        | bedmap \
            --skip-unmapped \
            --sweep-all \
            --range ${params.flank_width} \
            --delim "|" \
            --multidelim ";" \
            --echo \
            --echo-map <(cut -f1-6 ${variants} | grep -v '#') \
            - \
        | python3 $projectDir/bin/parse_variants_motifs.py \
            ${params.genome_fasta_file} \
            ${pwm_path} \
        >> ${counts_file}
    """
}

process tabix_index {
    conda params.conda
    publishDir params.outdir
    label "high_mem"
    scratch true

    input:
        path counts

    output:
        tuple path(name), path("${name}.tbi")
    script:
    name = "all_counts.merged.bed.gz"
    """
    head -1 ${counts} > tmp.bed
    sort-bed ${counts} >> tmp.bed
    bgzip -c tmp.bed > ${name}
    tabix ${name}
    """
}

process make_iupac_genome {
	conda "${params.conda}"
    tag "${prefix}"
	publishDir "${params.outdir}/alt_genome"

	input:
		val sample_id

	output:
		tuple path("${name}"), path("${name}.fai")

	script:
    prefix = sample_id ?: "all_samples"
	name = "${prefix}.iupac.genome.fa"
	additional_params = sample_id ? "--sample ${sample_id}" : ""
    """
    python3 $moduleDir/bin/nonref_genome.py \
        ${params.genome_fasta_file} \
        ${params.genotype_file} \
        ${name} \
        ${additional_params}
    """
}

process scan_with_moods {
    conda params.conda
    tag "${motif_id}"
    scratch true
    publishDir "${params.outdir}/moods_scans", pattern: "${name}"

    input:
        tuple val(motif_id), path(pwm_path), path(alt_fasta_file), path(fasta_index)

    output:
        tuple val(motif_id), path(pwm_path), path(name)
    
    script:
    name = "${motif_id}.moods.log.bed.gz"
    moods_params = file(params.bg_file).exists() ? "--lo-bg `cat ${params.bg_file}`" : ""
    """
    moods-dna.py --sep ";" -s ${alt_fasta_file} \
        --p-value ${params.motif_pval_tr} \
        ${moods_params} \
        -m "${pwm_path}" \
        -o moods.log

    
    cat moods.log | awk '{print \$1}' > chroms.txt

    cat moods.log \
        | cut -d";" -f2- \
        | sed 's/;\$//g' \
        | awk -v FS=";" -v OFS="\t" \
            '{ print \$2, \$2+length(\$5), \$1, \$4, \$3, \$5; }' \
        | sed 's/".pfm"/""/g' \
        | paste chroms.txt - \
        | sort-bed - \
        | bgzip -c \
        > ${name}
    """
}

process process_mutation_rates {
    tag "${vcf.simpleName}"
    scratch true
    conda params.conda
    label "med_mem"

    input:
        tuple path(vcf), path(variants_file)

    output:
        path name

    script:
    name = "${vcf.simpleName}.bed"
    """
    cut -f1-6 ${variants_file} > variants_pos.bed
    echo -e "#chr\tstart\tend\tID\tref\talt\tmut_rates_roulette\tmut_rates_gnomad" > ${name}
    
    bcftools query -f"%CHROM\t%POS0\t%POS\t%REF\t%ALT\t%INFO/MR\t%INFO/MG\n" ${vcf} \
        | awk '{print "chr"\$0}' \
        | bedtools intersect -a variants_pos.bed -b stdin -sorted -wa -wb \
        | awk -F'\t' '\$6 == \$11' \
        | cut -f1-6,12-13 >> ${name}
    """
}

process extract_context {
    conda params.conda
    scratch true

    input:
        path variants
    output:
        path name

    script:
    name = "variants_context.bed"
    """
    echo -e "#chr\tstart\tend\tsequence" > ${name}
    cat ${variants} \
        | awk -v OFS='\t' '{ print \$1,\$2-${params.window},\$3+${params.window} }' \
        | uniq > variants.bed 
    bedtools getfasta -fi ${params.genome_fasta_file} -bed variants.bed -bedOut \
        | awk -v OFS='\t' '{ print \$1,\$2+${params.window},\$3-${params.window},\$4 }' >> ${name}
    """
}

// Annotates with pheWAS, clinvar, finemapping, grasp, ebi-gwas phenotypes
process annotate_with_phenotypes {
    conda params.conda
    label "high_mem"

    input:
        path pval_file

    output:
        path name
        
    when:
        file(params.phenotypes_data, type: 'dir').exists()

    script:
    name = "phenotypes_ann.bed"
    """
    python3 $moduleDir/bin/annotate_with_phenotypes.py \
        ${params.phenotypes_data} \
        ${pval_file} \
        ${name}
    """
}

process merge_annotations {
    conda params.conda
    publishDir params.outdir
    scratch true

    input:
        path unique_snps
        path context
        path mutation_rates
    
    output:
        path name
    
    script:
    name = "snvs.annotations.bed.gz"
    """
    python3 $moduleDir/bin/merge_annotations.py \
        ${unique_snps} \
        ${context} \
        ${mutation_rates} \
        tmp.bed
    
    head -1 tmp.bed > res.bed
    sort-bed tmp.bed >> res.bed
    bgzip -c res.bed > ${name}
    """
}


workflow readMotifsList {
    main:
        scans = Channel.fromPath(params.motifs_list)
            | splitCsv(header:true, sep:'\t')
            | map(row -> tuple(row.motif, file(row.motif_file)))
    emit:
        scans
}

workflow mutationRates {
    take:
        data
    main:
        out = Channel.fromPath("${params.vcfs_dir}/*.vcf.gz")
            | combine(data)
            | process_mutation_rates
            | collectFile(
                name: "mut_rates.annotation.bed",
                skip: 1,
                keepHeader: true
            )
    emit:
        out
}

workflow motifCounts {
    take:
        data
    main:
        out = readMotifsList()
            | map(it -> tuple(it[0], it[1], file("${params.moods_scans_dir}/${it[0]}.moods.log.bed.gz")))
            | combine(data)
            | motif_counts
        out 
            | map(it -> it[1])
            | collectFile(
                name: "all.counts.bed",
                keepHeader: true,
                skip: 1
            ) 
            | tabix_index
    emit:
        out
}


// ------------ Entry workflows -------------------
workflow scanWithMoods {
    readMotifsList()
        | combine(make_iupac_genome())
        | scan_with_moods
}


workflow scanWithMoodsReference {
    readMotifsList()
        | combine(
            Channel.of(
                tuple(file(params.genotype_file), file("${params.genotype_file}.fai"))
            )
        )
        | scan_with_moods
}


workflow {
    params.moods_scans_dir = "${params.outdir}/moods_scans"
    extract_variants_from_vcf()
        | (motifCounts & annotate_with_phenotypes & extract_context & mutationRates)
    
    merge_annotations(
        annotate_with_phenotypes.out, 
        extract_context.out,
        mutationRates.out
    )
}



// DEFUNC
process make_dhs_annotation {
	conda "${params.conda}"
    tag "${ag_id}"
	publishDir "${params.outdir}/dhs_annotations"

	input:
		tuple val(indiv_id), val(ag_id), path(hotspots_file)

	output:
		path name

	script:
    name = "${ag_id}.dhs_annotation.bed"
    """
	cat ${params.genotype_annotation} | awk -v OFS='\t' '(\$7 == "${indiv_id}") { print; }' > gen_ann.bed
	unstarch ${hotspots_file} \
		| bedtools intersect -a stdin -b gen_ann.bed -wa -wb -sorted \
		| sed "s/\$/\t${ag_id}/" > sample.intersect.bed

	cat ${params.index_file} | awk -F'\t' -v OFS='\t' '{print \$1,\$2,\$3,\$4 }' \
		| bedtools intersect -a stdin -b sample.intersect.bed -wa -wb -sorted > ${name}
    """
}

workflow annotateDHS {
	// TODO add a step to the pipeline
	params.index_file = "/net/seq/data2/projects/ENCODE4Plus/indexes/index_altius_22-11-28/raw_masterlist/masterlist_DHSs_2902Altius-Index_nonovl_any_chunkIDs.bed"
	params.genotype_annotation = "${launchDir}/${params.outdir}/genotypes/genotypes_by_indiv_parsed.bed"
	Channel.fromPath(params.samples_file)
		| splitCsv(header:true, sep:'\t')
		| map(row -> tuple(row.indiv_id, row.ag_id, file(row.hotspots_file)))
		| make_dhs_annotation
		| collectFile(storeDir: "${params.outdir}",
			name: "genotypes_annotation.bed")
}

workflow tmp {
    params.moods_scans_dir = "${params.outdir}/moods_scans"
    Channel.fromPath("/net/seq/data2/projects/sabramov/ENCODE4/dnase-genotypes-round2/output/unique_variants.bed")
        | motifCounts
}