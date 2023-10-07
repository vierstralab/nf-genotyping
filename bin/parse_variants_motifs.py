#!/bin/env python3

import sys
from glob import glob
import os
import numpy as np
import pyfaidx
from collections import namedtuple


SNV = namedtuple('SNV', ['chrom', 'start', 'end',
         'name', 'score', 'strand', 'extra'])


def read_pfm(file):
    motif_id = os.path.splitext(os.path.basename(file))[0]
    pfm = np.loadtxt(file)
    pfm += 0.001
    pfm /= pfm.sum(axis=0)[np.newaxis,:]
    return motif_id, pfm


def complement(base):
    _comp = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    return _comp[base]

def seq_logp(mat, seq, bg=None):
    if bg is None:
        bg = [0.25, 0.25, 0.25, 0.25]
    res = 0
    for i, c in enumerate(seq):
        j = 'ACGT'.find(c)
        res += np.log(mat[j,i]/bg[j])
    return res

def ddg(seq, ref, alt, offset, pfm):
    assert len(seq) == pfm.shape[1] 
    assert offset >= 0
    
    ref_seq = list(seq)
    ref_seq[offset] = ref
    ref_score = seq_logp(pfm, ref_seq)
    
    alt_seq = list(seq)
    alt_seq[offset] = alt
    alt_score = seq_logp(pfm, alt_seq)

    return ref_score, alt_score
    
def main(fasta_path, pwm_path, input_iterator):
    fasta = pyfaidx.Fasta(fasta_path, sequence_always_upper=True)
    # Load PFM
    motif_id, pfm = read_pfm(pwm_path)

    "chrom,start,end,id,...|"
    for line in input_iterator:
        fields = line.strip('\n').split("|")

        snv_info = fields[0].split("\t")
        snv_chrom = snv_info[0]
        snv_start = int(snv_info[1])
        snv_end = int(snv_info[2])
        snv_dbsnp = snv_info[3] 
        snv_ref = snv_info[4]
        snv_alt = snv_info[5]

        elems=fields[1].split(";")

        variants_dict = {}
        for elem in elems:
            elem_info = elem.split("\t")
            chrom=str(elem_info[0])
            start=int(elem_info[1])
            end=int(elem_info[2])
            name=os.path.splitext(elem_info[3])[0]
            score=float(elem_info[4])
            strand=str(elem_info[5])
            seq=str(elem_info[6])

            overlaps_variant = start <= snv_start and end >= snv_end
            variants_dict.setdefault(name, []).append(
                SNV(chrom=chrom, start=start, end=end, name=name,
                    score=score, strand=strand, 
                    extra=overlaps_variant)
                )



        for motif, scores in variants_dict.items():
            nearest = sorted(
                scores, 
                key = lambda x: (x.extra, x.score),
                reverse=True
            )[0]
            pos = snv_start - nearest.start \
                if nearest.strand == '+' else nearest.end - snv_end

            # 'extra' is 1 if variant directly overlaps motif
            if nearest.extra:
                assert motif == motif_id

                seq = fasta[nearest.chrom][nearest.start:nearest.end]
                if nearest.strand == '-':
                    seq = seq.reverse.complement.seq
                    ref = complement(snv_ref)
                    alt = complement(snv_alt)
                else:
                    seq = seq.seq
                    ref = snv_ref
                    alt = snv_alt

                ref_score, alt_score = ddg(seq, ref, alt, pos, pfm)
            else:
                ref_score = alt_score = 0
                seq = '.'
            out = '\t'.join(map(str, [snv_chrom,
                                    snv_start,
                                    snv_end, 
                                    snv_dbsnp, 
                                    snv_ref, 
                                    snv_alt, 
                                    motif, 
                                    pos, 
                                    nearest.extra, 
                                    nearest.strand,
                                    ref_score,
                                    alt_score,
                                    seq
                                    ]))
            print(out)


if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2], sys.stdin)
