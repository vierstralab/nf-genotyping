#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

params.conda = "$moduleDir/environment.yml"

process motif_counts {
    // scratch true
    tag "${motif_id}"
    conda params.conda
    publishDir "${params.outdir}/motif_counts"

    input:
        tuple val(motif_id), path(pwm_path), path(moods_file), path(pval_file)

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
            --echo-map <(tail -n+2 ${pval_file}) \
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


process extract_variants_from_vcf {
    params.conda

    output:
        path name

    script:
    name = "unique_variants.bed"
    """
    bcftools query -f'%CHROM\t%POS0\t%POS\t%ID\t%REF\t%ALT\n' \
        ${params.genotype_file} > ${name}
    
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
    python3 $moduleDir/bin/nonref_genome.py ${params.genome_fasta_file} ${params.genotype_file} ${name} ${additional_params}
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


workflow readMotifsList {
    main:
        scans = Channel.fromPath(params.motifs_list)
            | splitCsv(header:true, sep:'\t')
            | map(row -> tuple(row.motif, file(row.motif_file)))
    emit:
        scans
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
        | motifCounts
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
