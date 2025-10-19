from src.files.fasta import FASTAFile 
import pandas as pd 
import numpy as np 
from Bio.Seq import Seq
import os


# 5' end is N-terminus and 3' end is C-terminus. 
is_partial_at_3_prime_end = lambda row : ((row.partial[-1] == '1') & (row.strand == 1)) | ((row.partial[0] == '1') & (row.strand == -1))
get_contig_id = lambda id_ :  '_'.join(id_.split('_')[:-1])

def get_stop_codons(fa_file=None, fn_file=None):

    contig_id_to_contig_map = dict(zip(fn_file.ids, fn_file.seqs))
    
    fa_df = fa_file.to_df(parse_description=True) # Convert to a DataFrame, parsing the Prodigal headers to get gene coordinates. 
    fa_df['contig_id'] = fa_df.index
    fa_df['contig_id'] = fa_df.contig_id.apply(get_contig_id) # Extract contig ID. 
    fa_df['start'], fa_df['stop'], fa_df['strand'] = fa_df['start'].astype(int), fa_df['stop'].astype(int), fa_df['strand'].astype(int)

    def get_stop_codon(row):
        if is_partial_at_3_prime_end(row):
            return 'none'
        contig = contig_id_to_contig_map[row.contig_id]
        nt_seq = contig[row.start - 1:row.stop] # Make sure to shift start to be zero-indexed. 
        nt_seq = str(Seq(nt_seq).reverse_complement()) if (row.strand == -1) else nt_seq
        return nt_seq[-3:]
    
    fa_df['stop_codon'] = fa_df.apply(get_stop_codon, axis=1)
    return fa_df


def get_stop_codon_info(fa_path:str=None, fn_path:str=None):

    fn_file = FASTAFile().from_fasta(fn_path) # Load in the FASTA nucleotide file. 
    fa_file = FASTAFile().from_fasta(fa_path) # Load in the Prodigal protein predictions. 
    fa_df = get_stop_codons(fa_file=fa_file, fn_file=fn_file)

    info = dict()
    info.update(fa_df.stop_codon.value_counts().to_dict())
    # info.update(fa_df.stop_codon.value_counts(normalize=True))
    info['gc_content'] = fn_file.get_gc_content()

    return info


def build_stop_codon_dataset(genome_ids:list, fn_dir:str='../dta/ncbi/genoms', fa_dir:str='../data/prodigal', path:str='../data/stop_codon.csv'):

    if not os.path.exists(path):
        stop_codon_df = list() 
        for genome_id in tqdm(genome_ids, 'build_stop_codon_dataset'):
            row = {'genome_id':genome_id}
            fn_path = os.path.join(fn_dir, f'{genome_id}.fn')
            fa_path = os.path.join(fa_dir, f'{genome_id}.fa')
            row.update(get_stop_codon_info(fn_path=fn_path, fa_path=fa_path))
            stop_codon_df.append(row)

        stop_codon_df = pd.concat(stop_codon_df)
        stop_codon_df['total'] = stop_codon_df.TAG + stop_codon_df.TAA + stop_codon_df.TGA 
        stop_codon_df.to_csv(path)




