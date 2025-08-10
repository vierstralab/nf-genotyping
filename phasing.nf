include { set_key_for_group_tuple } from "./post_processing"

process phasing {
    conda '/home/jvierstra/.local/miniconda3/envs/tf-2.13.0'
    publishDir "${params.outdir}/genotypes", pattern: "${prefix}.vcf.gz"
    label "medmem"

    input:
        tuple val(indiv_id), path(cram_files), path(cram_indices)

    output:
        tuple val(indiv_id), path(name)

    script:
    name = "${indiv_id}.phased.vcf"
    """
    bcftools view -s ${indiv_id} -e 'GT[*]="alt"' \
        ${params.genotype_file} -Oz > genotypes.vcf.gz
    
    bcftools index genotypes.vcf.gz

    whatshap phase \
        --ignore-read-groups \
        --sample ${indiv_id} \
        --reference ${params.genome_fasta_file} \
        -o ${name} \
        genotypes.vcf.gz \
        ${cram_files}
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