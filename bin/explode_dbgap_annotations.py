import pandas as pd
import sys

nucleotides = ('A', 'T', 'G', 'C')


def alt_str_has_single(alts_str):
    return sum([(len(alt) == 1) and (alt in nucleotides) for alt in alts_str.split(',')]) > 0


def main(snps, annotations, aa_anotation):   
    annotations = annotations[(annotations['dbgap'] != '.') &
                           (annotations['dbgap'].notna()) &
                           (annotations['alts'].apply(alt_str_has_single)) &
                           (annotations['ref'].isin(nucleotides))
    ].copy()
    
    annotations['n_alts'] = annotations['alts'].str.split(',').str.len()

    # Find duplicates
    dup_mask = annotations.duplicated(subset=['#chr', 'start', 'end', 'ref'], keep=False)
    dups = annotations[dup_mask].copy()

    # Classify each duplicate group
    def classify_group(group):
        alts = group['alts'].unique()
        if len(alts) == 1:
            return 'case1_same_alts_diff_freq'
        elif any(',' in a for a in alts) and any(',' not in a for a in alts):
            return 'case2_biallelic_vs_multiallelic'
        else:
            return 'case3_different_alts'

    case_counts = dups.groupby(['#chr', 'start', 'end', 'ref']).apply(classify_group).value_counts()
    print("Duplicate site classifications:", file=sys.stderr)
    print(f"  Case 1 - same alts, different frequencies: {case_counts.get('case1_same_alts_diff_freq', 0)}", file=sys.stderr)
    print(f"  Case 2 - biallelic vs multiallelic:        {case_counts.get('case2_biallelic_vs_multiallelic', 0)}", file=sys.stderr)
    print(f"  Case 3 - different alts entirely:          {case_counts.get('case3_different_alts', 0)}", file=sys.stderr)

    # Deduplicate keeping multiallelic first
    before = len(annotations)
    annotations = annotations.sort_values(['#chr', 'start', 'end', 'ref', 'n_alts'],
                                       ascending=[True, True, True, True, False])
    annotations = annotations.drop_duplicates(subset=['#chr', 'start', 'end', 'ref'], keep='first')
    annotations = annotations.drop(columns='n_alts')
    print(f"Dropped {before - len(annotations)} duplicate positions out of {before}", file=sys.stderr)



    assert annotations.value_counts(['#chr', 'start', 'end', 'ref']).max() == 1

    merged = snps.merge(annotations, 
        on=['#chr', 'start', 'end', 'ref'],
        how='left').merge(aa_anotation,
        on=['#chr', 'start', 'end', 'ref', 'alt'], how='left')
    
    merged['RAF'] = merged['dbgap'].apply(lambda x: '.' if pd.isna(x)
                                       else ('.' if x.split(',')[0] == '.' else float(x.split(',')[0])))
    merged['AAF'] = merged.apply(
        lambda row:
        '.' if row['RAF'] == '.' else
        dict(zip(row['alts'].split(','), row['dbgap'].split(',')[1:])).get(row['alt'], '.'),
        axis=1
    )

    return merged[['#chr', 'start', 'end', 'ref', 'alt', 'AAF', 'RAF', 'aa']]


if __name__ == '__main__':
    aa_anotation = pd.read_table(sys.argv[1],
        header=None, names=['#chr', 'start', 'end', 'ref', 'alt', 'aa'])

    dbsnp_annotation = pd.read_table(sys.argv[2],
        header=None, names=['#chr', 'start', 'end', 'rsid', 'ref', 'alts', 'dbgap'])
    dbsnp_annotation = dbsnp_annotation.drop(columns='rsid')
 
    snps_to_annotate = pd.read_table(sys.argv[3],
        header=None, names=['#chr', 'start', 'end', 'ref', 'alt'])
    df = main(snps_to_annotate, dbsnp_annotation, aa_anotation)
    df.to_csv(sys.argv[4], sep='\t', index=False, header=None)
