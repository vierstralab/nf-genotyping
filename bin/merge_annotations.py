import pandas as pd
import sys
from tqdm import tqdm
import numpy as np
tqdm.pandas()


mutations = ['C>A', 'C>T', 'C>G', 'T>A']
_comp = {
        'A': 'T',
        'C': 'G',
        'G': 'C',
        'T': 'A'
    }

def revcomp(s):
    return ''.join([_comp.get(x, 'N') for x in s[::-1]])

def palindromic(ref, alt):
    return {ref, alt} == {'A', 'T'} or {ref, alt} == {'G', 'C'}

def is_palindromic(ref, alt):
    
    return ref.mappalindromic_pairs.get() == alt

def get_mutation_stats(df):
    sequence = df['sequence'].str
    is_palindromic = df['ref'].map(_comp) == df['alt']
    sequence_minus1 = sequence[19]
    sequence_plus1 = sequence[21]
    assert np.all(df['ref'] == sequence[20])
    # is_palindromic == True
    df['fwd'] = (df["ref"] + '>' + df["alt"]).isin(mutations)

    is_palindromic3 = sequence_minus1.map(_comp) == sequence_plus1
    fwd3 = (sequence_minus1 != 'T') & (sequence_plus1 != 'T') & ~(sequence_minus1 == sequence_plus1 == 'G')
    df['ref_orient'] = np.where(
        is_palindromic & ~is_palindromic3,
        fwd3 == df['fwd'],
        True
    )
    df['sub'] = np.where(
        df['fwd'], 
        df["ref"] + '>' + df["alt"],
        df["alt"] + '>' + df["ref"]
    )

    # is_palindromic == False
    for mut, ref_orient, fwd  in [
        (df["ref"]+ '>' + df["alt"], True, True),
        (df["ref"].map(_comp)  + '>' + df["alt"].map(_comp), True, False),
        (df["alt"] + '>' + df["ref"], False, True),
        (df["alt"].map(_comp)  + '>' + df["ref"].map(_comp), False, False)
    ]:
        df.loc[~df['is_palindormic'] & mut.isin(mutations), 'sub'] = mut
        df.loc[~df['is_palindormic'] & mut.isin(mutations), 'ref_orient'] = ref_orient
        df.loc[~df['is_palindormic'] & mut.isin(mutations), 'fwd'] = fwd


    # final processing
    sequence_minus1_flipped = np.where(df['fwd'], sequence_minus1, sequence_plus1.map(_comp))
    sequence_plus1_flipped = np.where(df['fwd'], sequence_plus1, sequence_minus1.map(_comp))
    
    df['signature1'] = sequence_minus1_flipped + '[' + df['sub'] + ']' + sequence_plus1_flipped
        
    return df


def main(unique_snps, context, mutation_rates):
    result = unique_snps.merge(context)
    assert len(result.index) == len(unique_snps.index)
    
    result = result.merge(
        mutation_rates, how='left'
    )
    result = get_mutation_stats(result)
    result['cpg'] = ((result['sub'] != 'A>T') & (result['signature1'].str.endswith('G'))) | ((result['sub'] == 'C>G') & (result['signature1'].str.startswith('C')))

    return result


if __name__ == '__main__':
    unique_snps = pd.read_table(sys.argv[1])
    context = pd.read_table(sys.argv[2])
    mutation_rates = pd.read_table(sys.argv[3])
    
    main(unique_snps, context, mutation_rates).to_csv(
        sys.argv[4], sep='\t', index=False
    )
