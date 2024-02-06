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
fwd_mutations = mutations + [f"{mut[2]}>{mut[0]}" for mut in mutations]
nucleotide_pairs_priority = {
    'AA', 'CC', 'CA', 'AC', 'GA', 'AG'
}

import pandas as pd
import numpy as np

def reverse_complement(seq):
    seq_array = np.array(list(seq))
    comp_array = np.vectorize(_comp.get)(seq_array)
    rev_comp_array = comp_array[::-1]
    return ''.join(rev_comp_array)

def find_palindrome_length(seq):
    rev_comp_seq = reverse_complement(seq)
    length = len(seq)
    mid = length // 2

    for i in range(mid):
        if seq[mid + i + 1] != rev_comp_seq[mid + i + 1]:
            return i

    return mid


def get_mutation_stats(df, window_size):
    sequence = df['sequence'].str.upper()
    if not np.all(df['ref'] == sequence[window_size]):
        print(df['ref'], sequence[window_size], sequence)
        raise ValueError('Reference allele does not match the sequence')
    df['palindrome_length'] = sequence.apply(find_palindrome_length)
    is_palindromic = (df['ref'].map(_comp) == df['alt']) & (df['palindrome_length'] < window_size)
    palindrome_orient = np.where(
        is_palindromic.astype(int),
        (sequence[window_size - df['palindrome_length'] - 1]
        + sequence[window_size + df['palindrome_length'] + 1]
        ).isin(nucleotide_pairs_priority),
        True
    )
    df['fwd'] = (df["ref"] + '>' + df["alt"]).isin(fwd_mutations)
    fwd_sub = np.where(
        df['fwd'],
        df["ref"] + '>' + df["alt"],
        df["alt"].map(_comp) + '>' + df["ref"].map(_comp)
    )
    df['ref_orient'] = ~(np.isin(fwd_sub, mutations) ^ palindrome_orient)
    df['sub'] = np.where(
        df['ref_orient'], 
        fwd_sub,
        fwd_sub.map({
            sub: sub[2] + '>' + sub[0] for sub in mutations
        })
    )

    # final processing
    preceding1 = sequence[window_size - 1]
    following1 = sequence[window_size + 1]
    preciding1_flipped = np.where(df['fwd'], preceding1, following1.map(_comp))
    following1_flipped = np.where(df['fwd'], following1, preceding1.map(_comp))
    df['signature1'] = preciding1_flipped + '[' + df['sub'] + ']' + following1_flipped
    return df


def main(unique_snps, context, mutation_rates, window_size):
    result = unique_snps.merge(context)
    assert len(result.index) == len(unique_snps.index)

    result = get_mutation_stats(result, window_size)
    result['cpg'] = (
        (
            (result['sub'] != 'A>T') 
            & (result['signature1'].str.endswith('G'))
        ) | (
            (result['sub'] == 'C>G') 
            & (result['signature1'].str.startswith('C'))
            )
        )
    result = result.merge(
        mutation_rates, how='left'
    )
    return result


if __name__ == '__main__':
    unique_snps = pd.read_table(sys.argv[1])
    context = pd.read_table(sys.argv[2])
    mutation_rates = pd.read_table(sys.argv[3])
    window_size = int(sys.argv[4])
    
    main(unique_snps, context, mutation_rates, window_size).to_csv(sys.argv[5], sep='\t', index=False)
