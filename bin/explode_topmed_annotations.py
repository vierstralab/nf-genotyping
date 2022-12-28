import sys
import pandas as pd

def main(snps, annotations, aa_anotation):
    merged = snps.merge(annotations, 
        on=['chr', 'start', 'end', 'ref'],
        how='left').merge(aa_anotation,
        on=['chr', 'start', 'end', 'ref'], how='left')
    merged['aa'].fillna('.', inplace=True)
    merged['raf'] = merged['topmed'].apply(lambda x: '.' if pd.isna(x) or x == '.'
                                           else float(x.split(',')[0]))
    merged['aaf'] = merged.apply(
        lambda row:
        '.' if row['raf'] == '.' else
        dict(zip(row['alts'].split(','), row['topmed'].split(',')[1:]))[row['alt']],
        axis=1
    )
    merged[['chr', 'start', 'end', 'ref', 'alt', 'aaf', 'raf', 'aa']].to_csv(sys.stdout, sep='\t', index=False, header=None)


if __name__ == '__main__':
    dbsnp_annotation = pd.read_table(sys.argv[1],
        header=None, names=['chr', 'start', 'end', 'ref', 'alts', 'topmed'])
    aa_anotation = pd.read_table(sys.argv[2],
        header=None, names=['chr', 'start', 'end', 'ref', 'alt', 'aa'])
    snps_to_annotate = pd.read_table(sys.argv[2],
        header=None, names=['chr', 'start', 'end', 'ref', 'alt'])
    main(snps_to_annotate, dbsnp_annotation, aa_anotation)
