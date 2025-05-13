#!/usr/bin/env nextflow
include { extract_variants_from_vcf } from "./post_processing"

nextflow.enable.dsl = 2

params.conda = "$moduleDir/environment.yml"


process motif_counts {
    scratch true
    tag "${motif_id}:${genome_type}"
    conda params.conda
    publishDir "${params.outdir}/motif_counts_${genome_type}"

    input:
        tuple val(genome_type), val(motif_id), path(pwm_path), path(moods_file), path(variants)

    output:
        tuple val(genome_type), path(counts_file)

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
    tag "${genome_type}"
    label "highmem"
    scratch true

    input:
        tuple val(genome_type), path(counts)

    output:
        tuple path(name), path("${name}.tbi")

    script:
    name = "${genome_type}.motif_hits.merged.bed.gz"
    """
    head -1 ${counts[0]} > tmp.bed
    cat ${counts} \
        | awk '\$9 == "True"' > sign_hits.bed

    sort-bed sign_hits.bed >> tmp.bed

    bgzip -c tmp.bed > ${name}
    tabix ${name}
    """
}

process make_iupac_genome {
	conda "${params.conda}"
	publishDir "${params.outdir}/alt_genome"

	output:
		tuple val("alt"), path("${name}"), path("${name}.fai")

	script:
    // prefix = sample_id ?: ""
	name = "all_samples.iupac.genome.fa"
	// additional_params = sample_id ? "--sample ${sample_id}" : ""
    """
    python3 $moduleDir/bin/nonref_genome.py \
        ${params.genome_fasta_file} \
        ${params.genotype_file} \
        ${name}
    """
}

process scan_with_moods {
    conda params.conda
    tag "${motif_id}:${genome_type}"
    // scratch true
    publishDir "${params.outdir}/moods_scans_${genome_type}", pattern: "${name}"

    input:
        tuple val(genome_type), path(fasta_file), path(fasta_index), val(motif_id), path(pwm_path)

    output:
        tuple val(genome_type), val(motif_id), path(pwm_path), path(name)
    
    script:
    name = "${motif_id}.moods.log.bed.gz"
    moods_params = file(params.bg_file).exists() ? "--lo-bg `cat ${params.bg_file}`" : ""
    """
    moods-dna.py --sep ";" -s ${fasta_file} \
        --p-value ${params.motif_pval_tr} \
        ${moods_params} \
        -m "${pwm_path}" \
        -o moods.log


    cat moods.log \
        | awk -F";" -v OFS="\t" \
            '{ \
                split(\$1, a, " "); \
                chrom = a[1]; \
                print chrom, \$3, \$3+length(\$6), \$2, \$5, \$4, \$6; \
            }' \
        | sort-bed - \
        | bgzip -c \
        > ${name}
    """
}


workflow readMotifsList {
    main:
        scans = Channel.fromPath(params.motifs_list)
            | splitCsv(header:true, sep:'\t')
            | map(row -> tuple(row.motif_id, file(row.pwm)))
    emit:
        scans
}

workflow {
    Channel.of(tuple("ref", file(params.genome_fasta_file), file("${params.genome_fasta_file}.fai")))
        | mix(make_iupac_genome())
        | combine(readMotifsList())
        | scan_with_moods // genome_type, motif_id, pwm_path, moods_scans
        | combine(
            extract_variants_from_vcf()
        ) // genome_type, motif_id, pwm_path, moods_scans, variants
        | motif_counts // genome_type, counts_file
        | groupTuple() // genome_type, counts_files
        | tabix_index
}

workflow referenceScans {
    Channel.of(tuple("ref", file(params.genome_fasta_file), file("${params.genome_fasta_file}.fai")))
        | combine(readMotifsList())
        | scan_with_moods // genome_type, motif_id, pwm_path, moods_scans
}