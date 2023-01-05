#!/bin/env python

import shutil
import pyfaidx
import pysam
import argparse

iupac = "XACMGRSVTWYHKDBN"
def get_iupac(ref, alts):
	i = iupac.find(ref)
	for alt in alts:
		j = iupac.find(alt)
		i = i | j
	return iupac[i]


def main(orig_fasta_file, vcf_file, out_fasta_file, sample):
    print("Copying original fasta file...")

    shutil.copy2(orig_fasta_file, out_fasta_file)
    shutil.copy2(orig_fasta_file + ".fai", out_fasta_file + ".fai")

    print('Changing nucleotides to iupac')

    with pyfaidx.Fasta(out_fasta_file, mutable=True) as fasta:
        with pysam.VariantFile(vcf_file) as vcf:
            if sample is not None:
                if sample not in  vcf.header.samples:
                    raise ValueError('Error: Sample {} was not found in header'.format(sample))
        
            for rec in vcf.fetch():
                if any([len(i) > 1 for i in rec.alleles]):
                    continue
                
                ref = rec.ref
                alts = rec.alts

                if sample is not None:
                    if '/'.join(rec.samples[sample]['GT']) == '.':
                        continue
                #if len(ref)>1 or len(alt)>1:
                #	continue

                ambig = get_iupac(ref, alts)
                try:
                    fasta[rec.contig][rec.start] = ambig
                except KeyError as e:
                    pass

                    #print fasta[rec.contig][rec.start], ref, alt, ambig


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Make alternative genome from VCF file')
    parser.add_argument('fasta', help='Path to reference fasta file')
    parser.add_argument('vcf', help='Path to VCF file with variants')
    parser.add_argument('outpath', help='Path to fasta file with SNPs coded as IUPAC symbols')
    parser.add_argument('--sample', help='Sample name to extract from VCF file', default=None)
    args = parser.parse_args()
    main(args.fasta, args.vcf, args.outpath, args.sample)
