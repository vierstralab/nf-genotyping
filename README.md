# Epigenomics genotyping pipeline

Nextflow pipeline for genotyping from epigenomics data

## Reqiurements
- Nextflow (https://www.nextflow.io/)
- samtools (version X.X) (http://www.htslib.org/)
- bcftools (version 1.9) (http://www.htslib.org/)
- htslib (version 1.9) (http://www.htslib.org/)
- BEDOPS (version 2.4.35) ((http://www.htslib.org/))
- WASP (https://github.com/bmvdgeijn/WASP)

## Pipeline overview


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

## Output


