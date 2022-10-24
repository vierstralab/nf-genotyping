#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

params.conda = "$moduleDir/environment.yml"

process download_encode {
    publishDir "${params.outdir}/${encode_id}"
    conda params.conda
    maxForks 4
    cpus 2
    scratch true

    input:
        tuple val(encode_id), val(download_path), val(md5)

    output:
        tuple val(encode_id), path(name), path("${name}.crai")

    script:
    name = "${encode_id}.cram"
    """
    # download
    wget ${download_path} -o bamfile.bam
    # check md5
    if [[ `md5sum bamfile.bam | awk '{print \$1}'` != "${md5}" ]]; then
        exit 1
    fi
    # convert to cram
    samtools view bamfile.bam \
        -C -O cram,version=3.0,level=7,lossy_names=0 \
        -T "${params.genome_fasta_file}" \
        --threads "${task.cpus}" \
        --write-index \
        -o "${name}"
    """
}

workflow downloadEncode {
    take:
        metadata
    main:
        download_encode(metadata)
    emit:
        download_encode.out
}

params.encode_meta = "/home/sabramov/encode_chipseq.tsv"
workflow {
    metadata = Channel
		.fromPath(params.encode_meta)
		.splitCsv(header:true, sep:'\t')
		.map(row -> tuple( row.id, row.link, row.md5 ))
    
    downloadEncode(metadata)
}