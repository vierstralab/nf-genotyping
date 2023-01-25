import pandas as pd
import sys

nucleotides = ('A', 'T', 'G', 'C')


def alt_str_has_single(alts_str):
    return sum([(len(alt) == 1) and (alt in nucleotides) for alt in alts_str.split(',')]) > 0


def main(snps, annotations, aa_anotation):   
    annotations = annotations[(annotations['topmed'] != '.') &
                           (annotations['topmed'].notna()) &
                           (annotations['alts'].apply(alt_str_has_single)) &
                           (annotations['ref'].isin(nucleotides))
    ]

    assert annotations.value_counts(['#chr', 'start', 'end', 'ref']).max() == 1

    merged = snps.merge(annotations, 
        on=['#chr', 'start', 'end', 'ref'],
        how='left').merge(aa_anotation,
        on=['#chr', 'start', 'end', 'ref', 'alt'], how='left')
    merged['RAF'] = merged['topmed'].apply(lambda x: '.' if pd.isna(x)
                                           else float(x.split(',')[0]))
    merged['AAF'] = merged.apply(
        lambda row:
        '.' if row['RAF'] == '.' else
        dict(zip(row['alts'].split(','), row['topmed'].split(',')[1:])).get(row['alt'], '.'),
        axis=1
    )
    return merged[['#chr', 'start', 'end', 'ref', 'alt', 'AAF', 'RAF', 'aa']]


if __name__ == '__main__':
    aa_anotation = pd.read_table(sys.argv[1],
        header=None, names=['#chr', 'start', 'end', 'ref', 'alt', 'aa'])

    dbsnp_annotation = pd.read_table(sys.argv[2],
        header=None, names=['#chr', 'start', 'end', 'ref', 'alts', 'topmed'])

    snps_to_annotate = pd.read_table(sys.argv[3],
        header=None, names=['#chr', 'start', 'end', 'ref', 'alt'])
    df = main(snps_to_annotate, dbsnp_annotation, aa_anotation)
    df.to_csv(sys.argv[4], sep='\t', index=False, header=None)
