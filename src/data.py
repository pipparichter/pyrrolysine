import pandas as pd 
from src.files.fasta import FASTAFile
from tqdm import tqdm
import numpy as np 


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