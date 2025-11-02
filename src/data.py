import pandas as pd 
from src.files.fasta import FASTAFile
from tqdm import tqdm
import numpy as np 
import os
import re 



def get_arf1_data(hmmer_df:pd.DataFrame, genome_ids=None, data_dir:str='../data/prodigal'):
    if genome_ids is None:
        genome_ids = hmmer_df.genome_id.unique()
    hmmer_df = hmmer_df[hmmer_df.genome_id.isin(genome_ids) & hmmer_df.query_name.str.contains('aRF1')].copy()
    hmmer_df = hmmer_df.set_index('target_name')

    arf1_df = list()
    for genome_id, df in tqdm(list(hmmer_df.groupby('genome_id')), desc='get_arf1_data'):
        fasta_df = FASTAFile.from_fasta(os.path.join(data_dir, f'{genome_id}.fa')).to_df()
        fasta_df['length'] = fasta_df.seq.apply(len)
        # Remove all sequences which have long runs of unknown residues or run off the end of a contig.
        fasta_df['partial'] = ~fasta_df.description.str.contains('partial=00') | fasta_df.seq.str.contains('X{5,}', regex=True) 
        fasta_df = fasta_df[fasta_df.index.isin(df.index)]
        fasta_df = fasta_df.merge(df, left_index=True, right_index=True) # Add the HMMer data. 
        arf1_df.append(fasta_df)

    return pd.concat(arf1_df)

# def merge(coords:str, max_gap_size:int=0):
#     '''Some of the HMM hits overlap, which complicates analysis. Want to merge these hits to get runs of the 
#     target sequence which correspond to the HMM model. '''
#     coords = parse_coords(coords)
#     i = 0
#     while ((i + 1) < len(coords)):
#         gap_size = coords[i + 1][0] - coords[i][1]
#         if gap_size < max_gap_size:
#             merged_coord = _merge(coords.pop(i), coords.pop(i))
#             coords.insert(i, merged_coord)
#         else:
#             i += 1
#     coords = [(str(start), str(stop)) for start, stop in coords]
#     coords = ','.join(['..'.join(coord) for coord in coords])
#     return coords

# def _merge(coord_1, coord_2):
#     coords = list(coord_1) + list(coord_2)
#     return (min(coords), max(coords))