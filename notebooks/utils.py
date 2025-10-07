import matplotlib.pyplot as plt 
import pandas as pd
import numpy as np 
from tqdm import tqdm 
import glob
import os


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
    d = 0.015 
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

def hmmer_merge_domain_hits(df:pd.DataFrame):
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

    cols = ['target_name', 'accession', 'target_length', 'query_name', 'accession_', 'query_length', 'e_value', 'score', 'bias', 'domain_num', 'of', 'c_e_value', 'i_e_value',  'score_', 'bias_']
    cols += ['query_from', 'query_to', 'target_from', 'target_to', 'envelope_from', 'envelope_to',  'accuracy', 'description']

    # Omit the last two columns because parsing becomes difficult. 
    df = pd.read_csv(path, comment='#', header=None, usecols=np.arange(len(cols) - 2), names=cols[:-2], sep=r'\s+')
    df = df.drop(columns=['domain_num', 'of', 'accession', 'score_', 'bias_', 'accession_'])
    return hmmer_merge_domain_hits(df)


def hmmer_load(data_dir='../data/hmmer', max_e_value:float=1e-5):
    df = list()
    for path in tqdm(glob.glob(os.path.join(data_dir, '*')), desc='hmmer_load'):
        genome_id = os.path.basename(path).replace('.tab', '')

        try:
            df_ = _hmmer_load(path)
            df_ = df_.sort_values('e_value')
            df_ = df_.drop_duplicates('target_name') # This should resolve the cases where aPelota is being annotated as aRF-1.
            df_['genome_id'] = genome_id

            if max_e_value is not None:
                # Only filter the E-values for the Pyl machinery genes.
                df_ = df_[(df_.e_value < max_e_value) | df_.query_name.str.contains('aRF')].copy()
            if len(df_) > 0:
                df.append(df_)
        except:
            print(f'hmmer_load: Failed to load HMMer results for {path}.') # Fails on GCA_018302685.1, which has no hits for any HMM. 

    return pd.concat(df)