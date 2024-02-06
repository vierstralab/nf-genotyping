import pandas as pd
import sys
from tqdm import tqdm
import numpy as np
from matplotlib import pyplot as plt
tqdm.pandas()


mutations = ['C>A', 'C>T', 'C>G', 'T>A']
comp_vectorized = np.vectorize({
        'A': 'T',
        'C': 'G',
        'G': 'C',
        'T': 'A'
    }.get)
fwd_mutations = mutations + [f"{mut[2]}>{mut[0]}" for mut in mutations]
reverse_mapping = np.vectorize({sub: sub[2] + '>' + sub[0] for sub in fwd_mutations}.get)
nucleotide_pairs_priority = {
    'AA', 'CC', 'CA', 'AC', 'GA', 'AG'
}

def find_palindrome_lengths(sequences, mid):
    left_half = sequences[:, :mid]
    right_half = sequences[:, mid + 1:]

    rev_comp_right_half = comp_vectorized(right_half.astype(str)).astype('S1')[:, ::-1]

    mismatches = left_half != rev_comp_right_half
    first_mismatch_indices = np.argmax(mismatches, axis=1)
    
    no_mismatch = ~mismatches.any(axis=1)
    first_mismatch_indices[no_mismatch] = mid

    return first_mismatch_indices


def get_mutation_stats(df, window_size):
    sequence = df['sequence'].str
    if not np.all(df['ref'] == sequence[window_size]):
        print(df['ref'], sequence[window_size], sequence)
        raise ValueError('Reference allele does not match the sequence')
        
    char_array = np.frombuffer(''.join(df['sequence'].values).encode(), dtype='S1').reshape(len(df), -1)
    
    # df['palindrome_length'] = df['sequence'].apply(find_palindrome_length)
    df['palindrome_length'] = find_palindrome_lengths(char_array, window_size)
    is_palindromic = (comp_vectorized(df['ref']) == df['alt']) & (df['palindrome_length'] < window_size)
    
    preceding_index = window_size - df['palindrome_length'] - 1
    preceding_index = np.maximum(preceding_index, 0)
    following_index = window_size + df['palindrome_length'] + 1
    following_index = np.minimum(following_index, 2 * window_size)
    preceding_nucleotide = char_array[np.arange(len(char_array)), preceding_index].astype(str)
    following_nucleotide = char_array[np.arange(len(char_array)), following_index].astype(str)
    resolved_pair = np.core.defchararray.add(preceding_nucleotide, following_nucleotide)

    palindrome_orient = np.where(
        is_palindromic,
        np.isin(resolved_pair, nucleotide_pairs_priority),
        True
    )

    df['fwd'] = (df["ref"] + '>' + df["alt"]).isin(fwd_mutations)
    fwd_sub = np.where(
        df['fwd'],
        df["ref"] + '>' + df["alt"],
        np.char.add(comp_vectorized(df["alt"].values),
                    np.char.add('>', comp_vectorized(df["ref"].values))
                   )
    )
    df['ref_orient'] = ~(np.isin(fwd_sub, mutations).squeeze() ^ palindrome_orient)
    df['sub'] = np.where(
        df['ref_orient'], 
        fwd_sub,
        reverse_mapping(fwd_sub)
    )

    # final processing
    preceding1 = sequence[window_size - 1]
    following1 = sequence[window_size + 1]
    preciding1_flipped = np.where(df['fwd'], preceding1, comp_vectorized(following1))
    following1_flipped = np.where(df['fwd'], following1, comp_vectorized(preceding1))
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
