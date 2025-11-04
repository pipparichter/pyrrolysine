import pandas as pd 
import numpy as np 
from itertools import permutations, product


aas = 'ACDEFGHIKLMNPQRSTVWY' 
# aas = ''.join(list(dayhoff.values()))

def get_alphabet(order:int=1):
    alphabet = [''.join(pair) for pair in product(aas, repeat=order)]
    return sorted(alphabet)

def _get_frequencies(col, alphabet:list=None,):
    symbols, counts = np.unique(col, return_counts=True)
    counts = dict(zip(symbols, counts))
    counts = np.array([counts.get(a, 0) for a in alphabet])
    return counts / counts.sum() if (counts.sum() > 0) else np.array([np.nan] * len(alphabet))
    

def get_frequencies(alignment:np.ndarray, alphabet:list=None, order:int=1):
    alphabet = get_alphabet(order=order)
    n = alignment.shape[-1] # The alignment length, i.e. number of numeric columns. 
    f = np.empty(shape=[n]*order + [len(alphabet)], dtype=np.float16)

    for i in product(np.arange(n), repeat=order):
        col = np.array([''.join(residues) for residues in alignment.T[i]])
        # col = np.char.add.reduce(alignment.T[i], axis=0) # Concatenate over rows.
        f[i] = _get_frequencies(col, alphabet=alphabet)
    return f


def get_most_common_symbols(f:np.ndarray, order:int=1):
    alphabet = get_alphabet(order=order)
    return np.array(alphabet)[np.argmax(f, axis=-1)]

def get_scores(alignment_df, order:int=1, group_by:str='has_pyl', normalize:bool=True):

    max_score = np.sqrt(2)
    max_score = 2

    assert group_by in alignment_df.columns, f'get_scores: Specified grouping column, {group_by}, is not present.'
    groups = alignment_df[group_by].unique()
    assert len(groups) == 2, f'get_scores: Expected two groups, but got {', '.join(groups)}'

    n = sum([type(col) == int for col in alignment_df.columns]) # The alignment length, i.e. number of numeric columns. 

    f_dict = dict()
    for group, df in alignment_df.groupby(group_by):
        alignment = df[np.arange(n)].values # Extract the alignment array. 
        f_dict[group] = get_frequencies(alignment, order=order)

    # Comparing Euclidean distances between vectors. 
    scores = np.linalg.norm(f_dict[groups[1]] - f_dict[groups[0]], axis=-1, ord=1)
    scores = scores / max_score if normalize else scores
    return scores


def get_entropies(alignment_df):
    n = sum([type(col) == int for col in alignment_df.columns]) # The alignment length, i.e. number of numeric columns. 
    alignment = alignment_df[np.arange(n)].values # Extract the alignment array. 

    f = get_frequencies(alignment, order=1)
    symbols = get_most_common_symbols(f, order=1) 

    # Compute the entropy for each alignment column, which corresponds to a row in the frequency matrix. 
    # get_entropy = lambda col : 0 if (np.count_nonzero(col) == 1) else sum(-np.log(col[col != 0]) / np.log(np.count_nonzero(col)) * col[col != 0])
    get_entropy = lambda col : 0 if (np.count_nonzero(col) == 1) else sum(-np.log(col[col != 0]) / np.log(20) * col[col != 0])
    entropies = [get_entropy(col) for col in f] 
    return entropies, symbols 



# def get_mutual_information(alignment:np.ndarray):

#     length = alignment.shape[-1]
#     n = alignment.shape[0]

#     # First compute the frequency of each amino acid in the alphabet.
#     f = [{token.item():(col == token).mean() for token in np.unique(col)} for col in alignment.T]

#     # Then compute the co-ocurrence for each pair of tokens at each position. 
#     f_ab = [[None for _ in range(length)] for _ in range(length)]
#     for i in range(length):
#         for j in range(length):
#             pairs, counts = np.unique(np.char.add(alignment.T[i], alignment.T[j]), return_counts=True)
#             f_ab[i][j] = dict(zip(pairs, counts / n))

#     def h_i(i, f:np.ndarray):
#         return - np.sum([p * np.log2(p) for p in f[i].values()])

#     def h_ij(i, j, f_ab:np.ndarray=f_ab):
#         # k = len(f_ab[i][j]) # Get the alphabet size. 
#         return - np.sum([p * np.log2(p) for p in f_ab[i][j].values()])
    
#     h = np.empty(shape=(length, length))
#     for i in range(length):
#         for j in range(length):
#             if h_ij(i, j) == 0: # This will be zero for fully-conserved columns. 
#                 h[i, j] = np.nan 
#             else:
#                 h[i, j] = (h_i(i, f) + h_i(j, f) - h_ij(i, j)) / h_ij(i, j)
    
#     return h
