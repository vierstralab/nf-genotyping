#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

params.conda = "$moduleDir/environment.yml"
// Add sort-bed - 
process extract_variants_from_vcf {
    conda params.conda
    
    publishDir params.outdir

    output:
        path name

    script:
    name = "unique_variants.bed"
    """
    echo -e "#chr\tstart\tend\tID\tref\talt\tRAF\tAAF\tAA" > ${name}
    bcftools query -f'%CHROM\t%POS0\t%POS\t%ID\t%REF\t%ALT\t%INFO/RAF\t%INFO/AAF\t%INFO/AA\n' \
        ${params.genotype_file} >> ${name}
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
    label "highmem"

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

process distance_to_dhs {
    conda params.conda
    publishDir params.outdir
    label "highmem"

    input:
        path variants
    
    output:
        tuple path(genotyped_dist), path(name)

    script:
    genotyped_dist = "genotyped.distance_to_dhs.bed"
    name = "dbsnp_common.distance_to_dhs.bed"
    """
    tail -n+2 ${params.dhs_masterlist} \
        | cut -f1-4  > dhs_masterlist.bed
    
    echo -e "#chr\tstart\tend\tID\tref\talt\tRAF\tAAF\tAA\tindex_chr\tindex_start\tindex_end\tdhs_id\tdistance" > ${genotyped_dist}
    tail -n+2 ${variants} \
        | closest-features \
            --closest \
            --dist \
            --delim '\t' \
            - \
            dhs_masterlist.bed \
        >> ${genotyped_dist}
    
    echo -e "#chr\tstart\tend\tID\tref\talt\ttopmed\tmaf\tindex_chr\tindex_start\tindex_end\tdhs_id\tdistance\tis_genotyped" > ${name}
    
    tail -n+2 ${params.dbsnp_common_bed} \
        | closest-features \
        --closest \
        --dist \
        --delim '\t' \
        - \
        dhs_masterlist.bed \
        | bedmap \
            --delim '\t' \
            --echo \
            --indicator \
            - \
            <(tail -n+2 ${variants}) \
        >> ${name}
    """
}

process genomic_annotations {

    conda params.conda
    publishDir params.outdir
    scratch true

    input:
        path variants
    
    output:
        path name

    script:
    name = "snvs.genomic_annotations.bed"
    """
    tail -n+2 ${variants} \
        | sort-bed - > variants.no_header.bed
    
    bash $moduleDir/bin/gencodeAnnotations.sh \
        variants.no_header.bed \
        ${params.gencode} \
        ${params.genome_chrom_sizes}

    # Pasting all together
    paste <(head -1 ${variants}) <(cat gencodeAnnotations_header.txt) > ${name}
    paste variants.no_header.bed gencode_annotations.txt >> ${name}
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
        path genomic_annotations
    
    output:
        path name
    
    script:
    name = "snvs.annotations.bed.gz"
    """
    python3 $moduleDir/bin/merge_annotations.py \
        ${unique_snps} \
        ${context} \
        ${mutation_rates} \
        ${genomic_annotations} \
        1 \
        tmp.bed
    
    head -1 tmp.bed > res.bed
    sort-bed tmp.bed >> res.bed
    bgzip -c res.bed > ${name}
    """
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


// ------------ Entry workflows -------------------
workflow {
    extract_variants_from_vcf()
        | (
            annotate_with_phenotypes 
            & extract_context 
            & mutationRates 
            & distance_to_dhs
            & genomic_annotations
        )
    
    merge_annotations(
        annotate_with_phenotypes.out, 
        extract_context.out,
        mutationRates.out,
        genomic_annotations.out
    )
}


// Extracts the initial reads for all SNPs. This is used for WASP before/after figure
process extract_initial_reads {
    conda params.conda
    publishDir params.outdir

    output:
        path name

    script:
    name = "all.snps.initial_reads.bed.gz"
    """
    bcftools query -i'GT=="alt"' \
        -f '%CHROM\t%POS0\t%POS\t%ID\t%REF\t%ALT[\t%GT,%AD{0},%AD{1}]\n' \
        ${params.genotype_file} | bgzip > ${name}
    """
}


workflow extractInitialReadsRound1 {
    extract_initial_reads()
}


// CONVERT TO PLINK BED
process convert_to_plink_bed {
    conda params.conda
    publishDir "${params.outdir}/plink"
    cpus 8

    output:
        path("${prefix}.*")
    
    script:
    prefix = "plink.no_rsid"
    """
    bcftools view ${params.genotype_file} \
        | bcftools norm --threads ${task.cpus} -m-any \
            --check-ref w -f /home/jvierstra/data/1k_genomes/GRCh38_full_analysis_set_plus_decoy_hla.fa \
        | bcftools annotate --threads ${task.cpus} -x ID -I  +'%CHROM:%POS:%REF:%ALT' - \
        | bcftools norm --threads ${task.cpus} -Ob --rm-dup both - \
        > no_rsid.bcf
    
    plink \
        --bcf no_rsid.bcf \
        --keep-allele-order \
        --set-all-var-ids @: \
        --const-fid \
        --autosome \
        --allow-extra-chr \
        --split-x b38 no-fail \
        --make-bed \
        --out ${prefix}
    """
}

workflow convertToPlinkBed {
    convert_to_plink_bed()
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
