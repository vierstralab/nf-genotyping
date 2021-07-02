# Epigenomics genotyping pipeline

Nextflow pipeline for genotyping from epigenomics data

## Reqiurements
- Nextflow (https://www.nextflow.io/)
- samtools (http://www.htslib.org/)
- bcftools (http://www.htslib.org/)
- pyfaidx (https://github.com/mdshw5/pyfaidx)

## Pipeline overview

Samples BAM files are merged by corresponding individual and then used for a ``bcftools``-based genotyping pipeline.

## Input

<details><summary>Sample file [--samples_file]</summary>
<p>
	A tab-delimited file containing information about each sample.
</p>
</details>

<details><summary>Genome reference [--genome]</summary>
<p>
</p>
</details>

<details><summary>dbSNP reference [--dbsnp_file]</summary>
<p>
</p>
</details>

<details><summary>Ancestral genome [--genome_ancestral_fasta_file]</summary>
<p>
</p>
</details>



### Additonal Parameters:
<details><summary>Chunk size [--chunksize 5000000]</summary>
<p>
</p>
</details>

<details><summary>SNP quality [--min_SNPQ 10]</summary>
<p>
</p>
</details>

<details><summary>Genotype quality [--min_GQ 50]</summary>
<p>
</p>
</details>

<details><summary>Sequencing depth [--min_DP 12]</summary>
<p>
</p>
</details>

<details><summary>Hardy-Weinberg equilbrium [--hwe_cutoff 0.01]</summary>
<p>
</p>
</details>

<details><summary>Output directory [--outdir .]</summary>
<p>
</p>
</details>


## Output

The pipeline outputs a single VCF-formated file containing the called and filtered genotypes for each distinct invididual in the samples file. Each variant is annotated with the following extra infornation:

- ID field - dbSNP rs number
- INFO/CAF -- 1000 genomes project allele frequency (from dbSNP annotation file)
- INFO/TOPMED -- TOPMED project allele frequency (from dbSNP annotation file)
- INFO/AA -- Inferred ancenstral allele from EPO/PECAN alignments (see below for information about how this is obtained)


