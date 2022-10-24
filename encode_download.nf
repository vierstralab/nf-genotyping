#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

params.conda = "$moduleDir/environment.yml"

process download_encode {
    conda params.conda
    maxForks 1
    tag "${encode_id}"
    scratch true

    input:
        tuple val(encode_id), val(download_path), val(md5)

    output:
        tuple val(encode_id), path(name)

    script:
    name = "${encode_id}.bam"
    """
    wget ${download_path} -O ${name}
    if [[ `md5sum ${name} | awk '{print \$1}'` != "${md5}" ]]; then
        exit 1
    fi
    """
}

process convert_to_cram {
  tag "${sample_id}"
  publishDir "${params.outdir}/${sample_id}"
  cpus params.threads
  conda params.conda

  input:
    tuple val(sample_id), path(bam)

  output:
    tuple val(sample_id), path(cramfile), path("${cramfile}.crai")

  script:
  cramfile = bam.baseName + ".cram"
  """
  samtools index ${bam}

  samtools view "${bam}" \
    -C -O cram,version=3.0,level=7,lossy_names=0 \
    -T "${params.genome_fasta_file}" \
    --threads "${task.cpus}" \
    --write-index \
    -o "${cramfile}"
  """
}


workflow downloadEncode {
    take:
        metadata
    main:
        download_encode(metadata) | convert_to_cram
    emit:
        convert_to_cram.out
}

params.encode_meta = "/home/sabramov/encode_chipseq.tsv"
workflow {
    metadata = Channel
		.fromPath(params.encode_meta)
		.splitCsv(header:true, sep:'\t')
		.map(row -> tuple( row.id, row.link, row.md5 ))
    
    downloadEncode(metadata)
}