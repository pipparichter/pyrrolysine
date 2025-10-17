import pandas as pd 
import numpy as np 
from itertools import permutations, product


aas = 'ACDEFGHIKLMNPQRSTVWY' 
# aas = ''.join(list(dayhoff.values()))

def get_alphabet(order:int=1):
    # alphabet = list(permutations(aas, order))
    alphabet = [''.join(pair) for pair in product(aas, repeat=order)]
    # alphabet = [''.join(pair) for pair in alphabet]
    return sorted(alphabet)

def get_probabilities(col, alphabet:list=None):
    symbols, counts = np.unique(col, return_counts=True)
    counts = dict(zip(symbols, counts))
    counts = np.array([counts.get(a, 0) for a in alphabet])
    return counts / counts.sum() if (counts.sum() > 0) else np.array([np.nan] * len(alphabet))


# For each phylum, model each residue as a distribution of possible values, e.g. D = {A:0.4, T:0.3, ...} which should sum to one. 
# The average should still sum to one, as sum(D1, D2, ..., Dn) = n / n = 1

# def get_divergence_scores(alignment_df, alignment_length:int=alignment_length):

#     # alignment_df = alignment_df.replace(dayhoff).copy()
#     assert 'order' in alignment_df.columns, 'get_scores: Taxonomy needs to be specified.'
#     assert 'has_pyl' in alignment_df.columns, 'get_scores: Whether or not the organism is Pyl+ needs to be specified.'
#     dfs = dict()
#     for has_pyl, df in alignment_df.groupby('has_pyl'):
#         # How to merge across phyla? Maybe compute probabilities independently and then take the average so that
#         # each phylum contributes equally to the final probabilities. 
#         df = sum([get_probabilities(df_, alignment_length=alignment_length) for _, df_ in df.groupby('order')]) 
#         df = df / np.sum(df.values, keepdims=True, axis=1)
#         # df = df[df.sum(axis=1) != 0].copy()
#         dfs[has_pyl] = df

#     # scores = [chisquare(dfs[True].iloc[i].values.ravel(), dfs[False].iloc[i].values.ravel()) for i in range(alignment_length)]
#     scores = [np.linalg.norm(dfs[True].iloc[i].values.ravel() - dfs[False].iloc[i].values.ravel()) for i in range(alignment_length)]
#     # I want the pairwise distances between columns, so shape should be alignment_length, alignment_length
#     return scores


def _get_scores_order_1(alignment_df):

    assert 'has_pyl' in alignment_df.columns, 'get_scores: Whether or not the organism is Pyl+ needs to be specified.'
    n = sum([type(col) == int for col in alignment_df.columns]) # The alignment length, i.e. number of numeric columns. 
    
    alphabet = get_alphabet(order=1)
    f_dict = dict()
    for has_pyl, df in alignment_df.groupby('has_pyl'):
        f_dict[has_pyl] = np.array([get_probabilities(df[i].values, alphabet=alphabet) for i in range(n)])

    # Comparing Euclidean distances between column vectors. 
    scores = np.array([np.linalg.norm(f_dict[True][i] - f_dict[False][i]) for i in range(n)])
    return scores


# Want to support a similar computation for co-occurring pairs
def _get_scores_order_2(alignment_df):

    assert 'has_pyl' in alignment_df.columns, 'get_scores: Whether or not the organism is Pyl+ needs to be specified.'
    n = sum([type(col) == int for col in alignment_df.columns]) # The alignment length, i.e. number of numeric columns. 
    
    alphabet = get_alphabet(order=2)
    f_dict = dict()

    for has_pyl, df in alignment_df.groupby('has_pyl'):
        f = np.empty(shape=(n, n, len(alphabet)), dtype=np.float16)
        # Compute the co-ocurrence for each pair of tokens at each position. 
        for i in range(n):
            for j in range(n):
                col = np.char.add(df[i].values, df[j].values)
                f[i][j] = get_probabilities(col, alphabet=alphabet)
        f_dict[has_pyl] = f.copy()

    scores = np.empty((n, n))
    for i in range(n):
        for j in range(n):
            scores[i][j] = np.linalg.norm(f_dict[True][i][j] - f_dict[False][i][j])

    return scores


def get_scores(alignment_df, order:int=1):
    if order == 1:
        return _get_scores_order_1(alignment_df)
    elif order == 2:
        return _get_scores_order_2(alignment_df)