import matplotlib.pyplot as plt 
import pandas as pd
import numpy as np 
from tqdm import tqdm 
import glob
import os
import re 

blue = '#0000ff'
red = '#ff0000'
green = '#aaffaa'
gray = '#696969'

dayhoff = {'A':'A','G':'A','P':'A','S':'A','T':'A','D':'B','E':'B','N':'B','Q':'B','R':'C','H':'C','K':'C', 'M':'D','I':'D','L':'D','V':'D','F':'E','W':'E','Y':'E','C':'F' }

hp = {'A':'H','V':'H','L':'H','I':'H','P':'H','F':'H','W':'H','M':'H','Y':'H','G':'P','S':'P','T':'P','C':'P','N':'P','Q':'P','D':'P','E':'P','K':'P','R':'P','H':'P'}

ba = {'K':'B','R':'B','H':'B','D':'A','E':'A','A':'N','C':'N','F':'N','G':'N','I':'N','L':'N','M':'N','N':'N','P':'N','Q':'N','S':'N','T':'N','V':'N','W':'N','Y':'N'}


def make_itol_annotation_file(arf1_df:pd.DataFrame, palette=None, path:str=None, field:str=None, widths:dict=None):
    header_lines = ['TREE_COLORS', 'SEPARATOR COMMA', 'DATA']
    # NODE_ID TYPE COLOR LABEL_OR_STYLE SIZE_FACTOR
    lines = list()
    for row in arf1_df.itertuples():
        color = palette[getattr(row, field)]
        width = widths[getattr(row, field)]
        lines += [f'{row.Index},branch,{color},normal,{width}']
    lines = header_lines + lines 
    with open(path, 'w') as f:
        f.write('\n'.join(lines))



def mutual_information(msa_df:pd.DataFrame):

    length = len(msa_df.seq.iloc[0])

    # First compute the frequency of each amino acid in the alphabet.
    m = np.array([list(seq) for seq in msa_df.seq_a])
    f = [{token.item():(col == token).mean() for token in np.unique(col)} for col in m.T]

    # Then compute the co-ocurrence for each pair of tokens at each position. 

    # f_ab = [[None] * length_b] * length_a
    f_ab = [[None for _ in range(length)] for _ in range(length)]
    for i in range(length):
        for j in range(length):
            pairs = np.array([''.join(pair) for pair in (zip(m.T[i], m.T[j]))])
            f_ab[i][j] = {pair:(pairs == pair).mean() for pair in np.unique(pairs)}

    def h_i(i, f:np.ndarray):
        # k = len(f[i]) # Get the alphabet size. 
        return - np.sum([p * np.log2(p) for p in f[i].values()])

    def h_ij(i, j, f_ab:np.ndarray=f_ab, f_a:np.ndarray=f_a, f_b:np.ndarray=f_b):
        # k = len(f_ab[i][j]) # Get the alphabet size. 
        return - np.sum([p * np.log2(p) for p in f_ab[i][j].values()])
    
    h = np.empty(shape=(length, length))
    for i in range(length):
        for j in range(length):
            h[i, j] = (h_i(i, f) + h_i(j, f) - h_ij(i, j)) / h_ij(i, j)
    
    return h


def get_split_figure(bottom_range:tuple, top_range:tuple):

    height_ratio = (top_range[-1] - top_range[0]) / (bottom_range[-1] - bottom_range[0]) # Get ratio of top axes height to bottom axes height.

    fig_size = (5.4, 4.6)
    fig, (ax_top, ax_bottom) = plt.subplots(nrows=2, ncols=1, sharex=True, figsize=fig_size, height_ratios=[height_ratio, 1])
    
    ax_top.set_ylim(*top_range)
    ax_bottom.set_ylim(*bottom_range)

    # Hide the spines between the two plots. 
    ax_top.spines['bottom'].set_visible(False)
    ax_bottom.spines['top'].set_visible(False)
    ax_top.tick_params(labelbottom=False, bottom=False)  # Don't show tick labels on the top

    # Add diagonal lines to indicate the break. 
    d = 0.015 # * ((top_range[1] - top_range[0]) + (bottom_range[1] - bottom_range[0]))
    kwargs = dict(transform=ax_top.transAxes, color='k', clip_on=False)
    ax_top.plot((-d, +d), (-d, +d), **kwargs)
    ax_top.plot((1 - d, 1 + d), (-d, +d), **kwargs)

    kwargs.update(transform=ax_bottom.transAxes)
    ax_bottom.plot((-d, +d), (1 - d, 1 + d), **kwargs)
    ax_bottom.plot((1 - d, 1 + d), (1 - d, 1 + d), **kwargs)

    return fig, (ax_top, ax_bottom)


# The "ali coords" are the coordinates for the exact regions of the sequence that align to the profile HMM. 
# The "env coords" are the same as the alignment coordinates, but allowing for some uncertainty in the somain boundaries. 

# "acc" is the mean posterior probability of aligned residues in the maximum expected accuracy (MEA) alignment; a measure of how reliable 
# the overall alignment is (from 0 to 1, with 1.00 indicating a completely reliable alignment)

def _hmmer_merge(df:pd.DataFrame):
    '''When there are multiple hits for the same profile HMM on the same sequence, merge these into a single row; for example, if a sequence
    has separate hits for aRF-1 at the N-terminus and C-terminus, these are combined into a single entry.'''
    merged_df = list()
    for query_name, query_df in df.groupby('query_name'):
        for target_name, target_df in query_df.groupby('target_name'):
            row = dict()
            row['target_name'] = target_name 
            row['query_name'] = query_name
            row['e_value'] = target_df.e_value.min() # Use the best E-value. 
            row['query_coords'] = ','.join([f'{hit.query_from}..{hit.query_to}' for hit in target_df.itertuples()])
            row['target_coords'] = ','.join([f'{hit.target_from}..{hit.target_to}' for hit in target_df.itertuples()])
            row['envelope_coords'] = ','.join([f'{hit.envelope_from}..{hit.envelope_to}' for hit in target_df.itertuples()])
            row['target_length'] = target_df.target_length.iloc[0]
            row['query_length'] = target_df.query_length.iloc[0]
            row['n_hits'] = len(target_df) # Get the number of domain hits.
            merged_df.append(row)
    return pd.DataFrame(merged_df)


def _hmmer_load(path:str):
    '''Load a single HMM tabular output file, produced using --domtblout'''

    cols = ['target_name', 'accession', 'target_length', 'query_name', 'accession_', 'query_length', 'e_value', 'score', 'bias', 'domain_num', 'of', 'c_e_value', 'i_e_value',  'score_', 'bias_']
    cols += ['query_from', 'query_to', 'target_from', 'target_to', 'envelope_from', 'envelope_to',  'accuracy', 'description']

    # Omit the last two columns because parsing becomes difficult. 
    df = pd.read_csv(path, comment='#', header=None, usecols=np.arange(len(cols) - 2), names=cols[:-2], sep=r'\s+')
    df = df.drop(columns=['domain_num', 'of', 'accession', 'score_', 'bias_', 'accession_'])
    return _hmmer_merge(df)


# ['aRF1_eRF1', 'pelota', 'pyrrolys_PylC', 'pyrrolys_PylB', 'pyrrolys_PylD', 'PylS_Nterm']
def hmmer_load(data_dir='../data/hmmer', max_e_values:dict=dict(), query_names:list=['aRF1_eRF1', 'pelota', 'pyrrolys_PylC', 'pyrrolys_PylB', 'pyrrolys_PylD', 'PylS_Nterm', 'PylS_Cterm'], best_hit_only:bool=True):
    hmmer_df = list()
    for path in tqdm(glob.glob(os.path.join(data_dir, '*')), desc='hmmer_load'):
        genome_id = os.path.basename(path).replace('.tab', '')

        try:
            df_ = _hmmer_load(path).sort_values('e_value') # Sort so best E-values are first. 
            df_ = df_[df_.query_name.isin(query_names)]
            df_ = df_.drop_duplicates('target_name') if best_hit_only else df_ # Remove cases where multiple HMMs annotate the same sequence.
            df_ = df_.drop_duplicates('target_name') # This should resolve the cases where aPelota is being annotated as aRF-1.
            df_['genome_id'] = genome_id
            df_['max_e_value'] = df_.query_name.map(max_e_values)
            df_['max_e_value'] = df_.max_e_value.fillna(np.inf) # Set infinite upper bound it no E-value is specified. 
            df_ = df_[df_.e_value < df_.max_e_value]
            if len(df_) > 0:
                hmmer_df.append(df_)
        except:
            print(f'hmmer_load: Failed to load HMMer results for {path}.') # Fails on GCA_018302685.1, which has no hits for any HMM. 

    return pd.concat(hmmer_df)