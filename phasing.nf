include { set_key_for_group_tuple } from "./genotyping"

process phasing {
    conda params.conda
    publishDir "${params.outdir}/genotypes", pattern: "${prefix}.vcf.gz"
    label "medmem"
    // scratch true

    input:
        tuple val(indiv_id), path(cram_files), path(cram_indices)

    output:
        tuple val(indiv_id), path(name), path(bed_name), path("${bed_name}.tbi")

    script:
    name = "${indiv_id}.phased.vcf"
    bed_name = "${indiv_id}.phased.bed.gz"
    """
    bcftools view -s ${indiv_id} -e 'GT[*]="alt"' \
        ${params.genotype_file}  > genotypes.vcf

    whatshap phase \
        --ignore-read-groups \
        --sample ${indiv_id} \
        --reference ${params.genome_fasta_file} \
        -o ${name} \
        genotypes.vcf \
        ${cram_files}

    bcftools query -f "%CHROM\t%POS0\t%REF\t%ALT\t[%SAMPLE\t%GT\t%PS]\n" \
        ${name} | bgzip -c > ${bed_name}
    tabix ${bed_name}

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
		| phasing
}