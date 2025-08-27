include { set_key_for_group_tuple } from "./genotyping"


process merge_bam {
    conda params.conda
    publishDir "${params.outdir}/merged_bam"
    label "medmem"
    scratch true

    input:
        tuple val(indiv_id), path(bam_files), path(bam_indices)

    output:
        tuple val(indiv_id), path(name), path("${name}.bai")

    script:
    name = "${indiv_id}.merged.bam"
    """
    samtools merge -@ ${task.cpus} -f ${name} ${bam_files}
    samtools index ${name}
    """
}

process phasing {
    conda params.conda
    publishDir "${params.outdir}/phasing"
    // label "medmem"
    //scratch true
    tag "${indiv_id}"
    memory { 100.GB + 150.GB * task.attempt }

    input:
        tuple val(indiv_id), path(cram_files), path(cram_indices)

    output:
        tuple val(indiv_id), path(name), path(bed_name), path("${bed_name}.tbi")

    script:
    name = "${indiv_id}.phased.vcf.gz"
    bed_name = "${indiv_id}.phased.bed.gz"
    """
    bcftools view -s ${indiv_id} ${params.genotype_file} -Ou \
        | bcftools view -i 'GT[*]="alt"' -Ob \
         > genotypes.bcf

    bcftools index genotypes.bcf

    whatshap phase \
        --ignore-read-groups \
        --sample ${indiv_id} \
        --reference ${params.genome_fasta_file} \
        -o ${name} \
        genotypes.bcf \
        ${cram_files}

    echo "#chr\tstart\tend\tref\talt\tsample\tgenotype\tphase_set" > tmp.bed
    bcftools query \
        -f "%CHROM\t%POS0\t%POS\t%REF\t%ALT\t[%SAMPLE\t%GT\t%PS]\n" \
        ${name} >> tmp.bed

    bgzip -c tmp.bed > ${bed_name}
    tabix ${bed_name}
    """
}

process merge_bed {

    conda params.conda
    publishDir "${params.outdir}"

    label "medmem"

    input:
        path bed_files

    output:
        tuple path(name), path("${name}.tbi") 

    script:
    name = "all_phased.bed.gz"
    """
    # Keep header from first file
    ((zcat ${bed_files[0]} | head -n 1) || true) > merged.bed

    # Concat all, skip headers, and sort with sort-bed
    zcat ${bed_files} | grep -v '^#' | sort-bed - >> merged.bed

    # Compress and index if needed
    bgzip -c merged.bed > ${name}
    tabix ${name}
    """
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
        | groupTuple()
        | merge_bam
		| phasing
        | map(it -> it[2])
        | collect(sort: true, flat: true)
        | merge_bed
}

workflow mergeResults {
    Channel.fromPath(params.samples_file)
        | splitCsv(header:true, sep:'\t')
        | map(row -> tuple(
                row.indiv_id,
                file("${params.outdir}/phasing/${row.indiv_id}.phased.bed.gz"),
            )
        )
		| filter { !it[0].isEmpty() }
        | unique()
        | map(it -> it[1])
        | collect(sort: true, flat: true)
        | merge_bed
}