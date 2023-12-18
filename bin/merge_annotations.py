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

def get_mutation_stats(df, window_size):
    sequence = df['sequence'].str
    is_palindromic = df['ref'].map(_comp) == df['alt']
    preciding1 = sequence[window_size - 1]
    following1 = sequence[window_size + 1]
    assert np.all(df['ref'] == sequence[window_size])
    # is_palindromic == True
    df['fwd'] = (df["ref"] + '>' + df["alt"]).isin(mutations)

    is_palindromic3 = preciding1.map(_comp) == following1
    fwd3 = (following1 != 'T') & (following1 != 'T') & ~(preciding1 == following1 == 'G')
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
    preciding1_flipped = np.where(df['fwd'], preciding1, following1.map(_comp))
    following1_flipped = np.where(df['fwd'], preciding1, following1.map(_comp))
    df['signature1'] = preciding1_flipped + '[' + df['sub'] + ']' + following1_flipped
    return df


def main(unique_snps, context, mutation_rates, window_size):
    result = unique_snps.merge(context)
    assert len(result.index) == len(unique_snps.index)
    
    result = result.merge(
        mutation_rates, how='left'
    )
    result = get_mutation_stats(result, window_size)
    result['cpg'] = ((result['sub'] != 'A>T') & (result['signature1'].str.endswith('G'))) | ((result['sub'] == 'C>G') & (result['signature1'].str.startswith('C')))

    return result


if __name__ == '__main__':
    unique_snps = pd.read_table(sys.argv[1])
    context = pd.read_table(sys.argv[2])
    mutation_rates = pd.read_table(sys.argv[3])
    window_size = int(sys.argv[4])
    
    main(unique_snps, context, mutation_rates, window_size).to_csv(sys.argv[5], sep='\t', index=False)
