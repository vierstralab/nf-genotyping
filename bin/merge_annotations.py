import pandas as pd
import sys
from tqdm import tqdm
import numpy as np
from matplotlib import pyplot as plt
tqdm.pandas()


mutations = ['C/A', 'C/T', 'C/G', 'T/A']
comp = {
        'A': 'T',
        'C': 'G',
        'G': 'C',
        'T': 'A'
    }

comp_vectorized = np.vectorize(comp.get)
fwd_mutations = mutations + [f"{mut[2]}/{mut[0]}" for mut in mutations]
reverse_mapping = np.vectorize({sub: sub[2] + '/' + sub[0] for sub in set(fwd_mutations) - set(mutations)}.get)
nucleotide_pairs_priority = [
    'AA', 'CC', 'CA', 'AC', 'GA', 'AG'
]

def rc(char_array):
    return comp_vectorized(char_array.astype(str)).astype('S1')[:, ::-1]

def find_palindrome_lengths(sequences):
    mid = sequences.shape[1] // 2

    left_half = sequences[:, :mid]
    right_half = sequences[:, mid + 1:]

    rev_comp_left_half = rc(left_half)

    mismatches = rev_comp_left_half != right_half
    first_mismatch_indices = np.argmax(mismatches, axis=1)
    
    no_mismatch = ~mismatches.any(axis=1)
    first_mismatch_indices[no_mismatch] = mid

    return first_mismatch_indices


def transform_into_str(char_array):
    return np.copy(char_array, order='C').view('S' + str(char_array.shape[1])).ravel().astype(str)


def get_mutation_stats(df, window_size):
    l = len(df['sequence'].iloc[0]) # assume that all length are the same
    window_size = min(window_size, l//2)
    df['short_sequence'] = df['sequence'].str[l//2 - window_size: l//2 + window_size + 1]
    sequence = df['short_sequence'].str
    if not np.all(df['ref'] == sequence[window_size]):
        print(df['ref'], sequence[window_size], sequence)
        raise ValueError('Reference allele does not match the sequence')
        
    char_array = np.frombuffer(''.join(df['short_sequence'].values).encode(), dtype='S1').reshape(len(df), -1)
    
    df['palindrome_length'] = find_palindrome_lengths(char_array)
    is_palindromic = (df['ref'].map(comp) == df['alt']) & (df['palindrome_length'] < window_size)
    
    preceding_index = window_size - df['palindrome_length'] - 1
    preceding_index = np.maximum(preceding_index, 0)
    following_index = window_size + df['palindrome_length'] + 1
    following_index = np.minimum(following_index, 2 * window_size)
    preceding_nucleotide = char_array[np.arange(len(char_array)), preceding_index].astype(str)
    following_nucleotide = char_array[np.arange(len(char_array)), following_index].astype(str)
    resolved_pair = np.core.defchararray.add(preceding_nucleotide, following_nucleotide).astype(str)

    palindrome_orient = np.where(
        is_palindromic,
        np.isin(resolved_pair, nucleotide_pairs_priority),
        True
    )

    initial_fwd = (df["ref"] + '/' + df["alt"]).isin(fwd_mutations)
    df['fwd'] = initial_fwd & palindrome_orient

    
    fwd_sub = np.where(
        (df["ref"] + '/' + df["alt"]).isin(fwd_mutations),
        df["ref"] + '/' + df["alt"],
        df["ref"].map(comp) + '/' + df["alt"].map(comp)
    )
    df['ref_orient'] = np.isin(fwd_sub, mutations)
    df['sub'] = np.where(
        df['ref_orient'], 
        fwd_sub,
        reverse_mapping(fwd_sub)
    )

    # final processing
    preceding_rc = transform_into_str(rc(char_array[:, :window_size]))
    following_rc = transform_into_str(rc(char_array[:, window_size + 1:]))
    preciding1_flipped = np.where(df['fwd'], sequence[:window_size], following_rc)
    following1_flipped = np.where(df['fwd'], sequence[window_size+1:], preceding_rc)
    df[f'signature{window_size}'] = preciding1_flipped + '[' + df['sub'] + ']' + following1_flipped
    return df


def main(unique_snps, context, mutation_rates, window_size):
    result = unique_snps.merge(context)
    assert len(result.index) == len(unique_snps.index)

    result = get_mutation_stats(result, window_size)
    result['cpg'] = (
        (
            (result['sub'] != 'T/A') 
            & (result[f'signature{window_size}'].str.contains('\]G'))
        ) | (
            (result['sub'] == 'C/G') 
            & (result[f'signature{window_size}'].str.contains('C\['))
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
