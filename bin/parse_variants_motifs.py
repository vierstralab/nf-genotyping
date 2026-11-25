#!/bin/env python3

import os
import argparse
import numpy as np
import pandas as pd
from genome_tools.genomic_interval import GenomicInterval, VariantInterval
from genome_tools.data.extractors import FastaExtractor
from genome_tools.data.pwm import read_pfm, get_allelic_scores, reverse_complement


HEADER = ['#chr', 'start', 'end', 'ID', 'ref', 'alt', 'motif', 'offset',
          'within', 'strand', 'ref_score', 'alt_score', 'seq']


def parse_variant(variant_str):
    """chrom, start, end, id, ref, alt"""
    snv_info = variant_str.split("\t")
    return VariantInterval(
        chrom=snv_info[0],
        start=int(snv_info[1]),
        end=int(snv_info[2]),
        ref=snv_info[4],
        alt=snv_info[5],
        name=snv_info[3],
    )


def parse_motif_hit(hit_str, variant):
    """chrom, start, end, motif, score, strand, seq (MOODS hit)"""
    elem_info = hit_str.split("\t")
    hit = GenomicInterval(
        chrom=elem_info[0],
        start=int(elem_info[1]),
        end=int(elem_info[2]),
        name=os.path.splitext(elem_info[3])[0],
        score=float(elem_info[4]),
        strand=elem_info[5],
    )
    hit.overlaps_variant = variant.within(hit)
    return hit


def group_hits_by_motif(hits_str, variant):
    hits_by_motif = {}
    for elem in hits_str.split(";"):
        hit = parse_motif_hit(elem, variant)
        hits_by_motif.setdefault(hit.name, []).append(hit)
    return hits_by_motif


def motif_offset(hit, variant):
    """Position of the variant within the motif, in motif orientation"""
    if hit.strand == '+':
        return variant.start - hit.start
    return hit.end - variant.end


def score_hit(fasta, hit, variant, pfm, bg=None):
    """Returns ref and alt scores and the motif sequence in motif orientation"""
    seq = fasta[hit].upper()
    ref_score, alt_score = get_allelic_scores(
        pfm_matrix=pfm,
        sequence=seq,
        ref=variant.ref,
        alt=variant.alt,
        offset=variant.start - hit.start,
        orient=hit.strand,
        bg=bg,
    )
    if hit.strand == '-':
        seq = reverse_complement(seq)
    return ref_score, alt_score, seq


def main(fasta_path, pfm, motif_id, input_path, bg=None):
    rows = []
    with FastaExtractor(fasta_path) as fasta, open(input_path) as input_file:
        # input: chrom,start,end,id,ref,alt|hit1;hit2;... (bedmap --echo --echo-map)
        for line in input_file:
            variant_str, hits_str = line.strip('\n').split("|")
            variant = parse_variant(variant_str)
            hits_by_motif = group_hits_by_motif(hits_str, variant)

            for motif, hits in hits_by_motif.items():
                # prefer hits overlapping the variant, then the highest MOODS score
                hit = max(hits, key=lambda x: (x.overlaps_variant, x.score))

                if hit.overlaps_variant:
                    if motif != motif_id:
                        raise ValueError(
                            f"Motif {motif} in input does not match PWM {motif_id} "
                            f"(variant {variant.to_str()})"
                        )
                    ref_score, alt_score, seq = score_hit(fasta, hit, variant, pfm, bg=bg)
                else:
                    ref_score = alt_score = 0
                    seq = '.'

                rows.append({
                    '#chr': variant.chrom,
                    'start': variant.start,
                    'end': variant.end,
                    'ID': variant.name,
                    'ref': variant.ref,
                    'alt': variant.alt,
                    'motif': motif,
                    'offset': motif_offset(hit, variant),
                    'within': hit.overlaps_variant,
                    'strand': hit.strand,
                    'ref_score': ref_score,
                    'alt_score': alt_score,
                    'seq': seq,
                })
    return pd.DataFrame(rows, columns=HEADER)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Score ref and alt alleles of variants overlapping motif hits'
    )
    parser.add_argument('fasta', help='Path to reference fasta file')
    parser.add_argument('pwm', help='Path to PFM file; file name is used as motif id')
    parser.add_argument('input', help='bedmap output: variant|motif hits within flank')
    parser.add_argument('outpath', help='Path to output counts file')
    parser.add_argument('--bg', help='MOODS background file: A C G T frequencies',
                        default=None)
    args = parser.parse_args()

    bg = np.loadtxt(args.bg) if args.bg is not None else None
    motif_id = os.path.splitext(os.path.basename(args.pwm))[0]
    pfm = read_pfm(args.pwm)
    df = main(args.fasta, pfm, motif_id, args.input, bg=bg)
    df.to_csv(args.outpath, sep='\t', index=False)
