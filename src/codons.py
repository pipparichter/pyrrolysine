from src.files.fasta import FASTAFile 
import pandas as pd 
import numpy as np 
from Bio.Seq import Seq


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



