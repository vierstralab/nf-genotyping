#!/bin/env python

import shutil
import pyfaidx
import pysam
import argparse


iupac = "XACMGRSVTWYHKDBN"
def get_iupac(alleles):
    i = None
    for allele in alleles:
        j = iupac.find(allele)
        if i is None:
            i = j
        else:
            i = i | j
    return iupac[i]


def get_sample_gt(rec, sample):
    result = set()
    gt = rec.samples[sample]['GT']
    if gt is None or gt[0] is None or gt[0] == '.':
        return result
    for allele_count in gt:
        if allele_count == 0:
            result.add(rec.ref)
        elif allele_count == 1:
            assert len(rec.alts) == 1
            assert len(rec.alts[0]) == 1
            result.add(rec.alts[0])
        else:
            raise AssertionError
    return result
    

def main(orig_fasta_file, vcf_file, out_fasta_file, sample):
    print("Copying original fasta file...")

    shutil.copy2(orig_fasta_file, out_fasta_file)
    shutil.copy2(orig_fasta_file + ".fai", out_fasta_file + ".fai")

    print('Changing nucleotides to iupac')
    old_key = None
    with pyfaidx.Fasta(out_fasta_file, mutable=True) as fasta:
        with pysam.VariantFile(vcf_file) as vcf:
            if sample is not None:
                if sample not in vcf.header.samples:
                    raise ValueError('Error: Sample {} was not found in header'.format(sample))
        
            for rec in vcf.fetch():
                if any([len(i) > 1 for i in rec.alleles]):
                    continue

                ref = rec.ref
                alts = rec.alts
                ambig = get_iupac(set([ref, *alts]))
                try:
                    fasta[rec.contig][rec.start] = ambig
                except KeyError as e:
                    print('KeyError: {}'.format(e))
                    pass
                # key = f'{rec.contig}@{rec.start}'
                # if old_key is None or key != old_key:
                #     alleles = set()
                
                # if sample is not None:
                #     alleles |= get_sample_gt(rec, sample)
                # else:
                #     for sample in vcf.header.samples:
                #         alleles |= get_sample_gt(rec, sample)

                # old_key = key
                # if len(alleles) > 0:
                #     ambig = get_iupac(alleles)
                #     try:
                #         fasta[rec.contig][rec.start] = ambig
                #     except KeyError as e:
                #         pass


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Make alternative genome from VCF file')
    parser.add_argument('fasta', help='Path to reference fasta file')
    parser.add_argument('vcf', help='Path to VCF file with variants')
    parser.add_argument('outpath', help='Path to fasta file with SNPs coded as IUPAC symbols')
    parser.add_argument('--sample', help='Sample name to extract from VCF file', default=None)
    args = parser.parse_args()
    main(args.fasta, args.vcf, args.outpath, args.sample)
 