from src.files.fasta import FASTAFile, parse_prodigal_description, get_contig_id
import pandas as pd 
import numpy as np 
from Bio.Seq import Seq
import os
from tqdm import tqdm 



# 5' end is N-terminus and 3' end is C-terminus. 
# is_partial_at_3_prime_end = lambda row : ((row.partial[-1] == '1') & (row.strand == 1)) | ((row.partial[0] == '1') & (row.strand == -1))
is_partial_at_3_prime_end = lambda partial, strand : ((partial[-1] == '1') & (strand == 1)) | ((partial[0] == '1') & (strand == -1))


def _get_nt_seq(start:int=None, stop:int=None, strand:int=None, contig:str=None, partial:str=None, contig_id:str='', fn_path:str='', **kwargs):
    if is_partial_at_3_prime_end(partial, strand):
        return 'none'
    assert stop <= len(contig), f'_get_nt_seq: Specified stop position {stop} is greater than the length of contig {contig_id} {len(contig)} in file {fn_path}'
    assert start < len(contig), f'_get_nt_seq: Specified start position {start} is greater than the length of contig {contig_id} {len(contig)} in file {fn_path}'
    
    nt_seq = contig[start - 1:stop] # Make sure to shift start to be zero-indexed. 
    nt_seq = str(Seq(nt_seq).reverse_complement()) if (strand == -1) else nt_seq
    return nt_seq[-3:]


def _get_stop_codons_in_file(fa_file=None, fn_file=None):

    contig_id_to_contig_map = dict(zip(fn_file.ids, fn_file.seqs))
    fa_df = fa_file.to_df(parse_description=True) # Convert to a DataFrame, parsing the Prodigal headers to get gene coordinates. 
    fa_df['stop_codon'] = fa_df.apply(lambda row : _get_nt_seq(**row.to_dict(), contig=contig_id_to_contig_map[row.contig_id]), axis=1)
    return fa_df


def _get_stop_codons_in_dataframe(fa_df:pd.DataFrame):

    required_cols = np.array(['start', 'stop', 'strand', 'fn_path', 'partial', 'contig_id'])
    missing_cols = required_cols[~np.isin(required_cols, fa_df.columns)].tolist()
    assert len(missing_cols) == 0, f'_get_stop_codons_in_fa_dataframe: Input DataFrame is missing required columns {', '.join(missing_cols)}.'

    def get_contig(fn_path:str, contig_id:str):
        fn_file = FASTAFile().from_fasta(fn_path)
        assert contig_id in fn_file.ids, f'_get_stop_codons_in_dataframe: Contig {contig_id} is missing from the file {fn_path}'
        return fn_file.seqs[fn_file.ids == contig_id][0]

    # Accumulate required contigs across all FASTA nucleotide files specified in the DataFrame.
    contig_id_to_contig_map = {row.contig_id:get_contig(row.fn_path, row.contig_id) for row in fa_df.itertuples()}
    fa_df['stop_codon'] = fa_df.apply(lambda row : _get_nt_seq(**row.to_dict(), contig=contig_id_to_contig_map[row.contig_id]), axis=1)
    assert np.all(fa_df.stop_codon.apply(len) > 0), '_get_stop_codons_in_dataframe: Some of the stop codons are empty strings.'
    
    return fa_df


def _get_stop_codon_genome_metadata_info(fa_path:str=None, fn_path:str=None):

    fn_file = FASTAFile().from_fasta(fn_path) # Load in the FASTA nucleotide file. 
    fa_file = FASTAFile().from_fasta(fa_path) # Load in the Prodigal protein predictions. 
    fa_df = _get_stop_codons_in_file(fa_file=fa_file, fn_file=fn_file)

    info = dict()
    info.update(fa_df.stop_codon.value_counts().to_dict())
    # info.update(fa_df.stop_codon.value_counts(normalize=True))
    info['gc_content'] = fn_file.get_gc_content()

    return info


def build_stop_codon_genome_metadata_dataset(genome_ids:list, fn_dir:str='../dta/ncbi/genomes', fa_dir:str='../data/prodigal', path:str='../data/stop_codon_genome_metadata.csv'):
    
    if not os.path.exists(path):
        stop_codon_df = list() 
        for genome_id in tqdm(genome_ids, 'build_stop_codon_dataset'):
            row = {'genome_id':genome_id}
            fn_path = os.path.join(fn_dir, f'{genome_id}.fn')
            fa_path = os.path.join(fa_dir, f'{genome_id}.fa')
            row.update(_get_stop_codon_genome_metadata_info(fn_path=fn_path, fa_path=fa_path))
            stop_codon_df.append(row)

        stop_codon_df = pd.concat(stop_codon_df)
        stop_codon_df['total'] = stop_codon_df.TAG + stop_codon_df.TAA + stop_codon_df.TGA 
        stop_codon_df.to_csv(path)



def build_stop_codon_dataset(df:pd.DataFrame, fn_dir:str='../data/ncbi/genomes', path:str=None):
    
    fa_df = pd.DataFrame([parse_prodigal_description(description) for description in df.description])
    fa_df['contig_id'] = [get_contig_id(id_) for id_ in df.index]
    fa_df['fn_path'] = [os.path.join(fn_dir, f'{genome_id}.fn') for genome_id in df.genome_id]
    fa_df['start'], fa_df['stop'], fa_df['strand'] = fa_df['start'].astype(int), fa_df['stop'].astype(int), fa_df['strand'].astype(int)
    fa_df.index = df.index 
    fa_df = _get_stop_codons_in_dataframe(fa_df)

    assert np.all(fa_df.index == df.index), 'build_stop_codon_dataset: Indices do not align.'

    fa_df = fa_df.drop(columns=df.columns, errors='ignore')
    stop_codon_df = df.merge(fa_df, left_index=True, right_index=True)
    stop_codon_df.to_csv(path)



    


