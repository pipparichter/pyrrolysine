import os 
import re 
from typing import List, Dict, Tuple, NoReturn
from tqdm import tqdm 
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import pandas as pd 


# def parser_prodigal(description:str):
#     pattern = r'# ([\d]+) # ([\d]+) # ([-1]+) # ID=([^;]+);partial=([^;]+);start_type=([^;]+);rbs_motif=([^;]+);rbs_spacer=([^;]+);gc_cont=([\.\w]+)'
#     columns = ['start', 'stop', 'strand', 'ID', 'partial', 'start_type', 'rbs_motif', 'rbs_spacer', 'gc_content']
#     match = re.search(pattern, description)
#     parsed_header = {col:match.group(i + 1) for i, col in enumerate(columns)}
#     parsed_header['rbs_motif'] = 'none' if (parsed_header['rbs_motif'] == 'None') else parsed_header['rbs_motif']
#     return parsed_header


# def parser_default(description:str):
#     return {'description':description}


class FASTAFile():
    # prodigal_dtypes = {'start':int, 'stop':int, 'strand':int, 'gc_content':float, 'rbs_motif':str, 'rbs_spacer':str}

    def __init__(self):
        '''Initialize a FASTAFile object.'''

        return


    def __len__(self):
        return len(self.seqs)
    
    @classmethod
    def from_df(cls, df:pd.DataFrame):
        obj = cls()

        obj.seqs = df.seq.values
        obj.ids = df.index.values 
        obj.descriptions = df.description.values if ('description' in df.columns) else [''] * len(obj.ids)
        return obj
        
    @classmethod
    def from_fasta(cls, path:str):
        obj = cls()
        f = open(path, 'r')
        obj.seqs, obj.ids, obj.descriptions = [], [], []
        for record in SeqIO.parse(path, 'fasta'):
            obj.ids.append(record.id)
            obj.descriptions.append(record.description.replace(record.id, '').strip())
            obj.seqs.append(str(record.seq))
        f.close()
        return obj 


    # @staticmethod
    # def parse_coordinate(coord:str):
    #     strand = '-' if 'complement' in coord else '+'
    #     pattern = r'[<>]*(\d+)..(\d+)[<>]+'
    #     match_ = re.search(pattern, coord)

    #     if match_ is None:
    #         return None 
    #     else:
    #         return int(match_.group(1)), int(match_.group(2)), strand

    # def from_genbank(cls, path:str, ):

    #     with open(path, 'r') as f:
    #         content = f.read()

    #     contig_ids = [match_.group(1) for match_ in re.finditer('seqhdr="([^\s]+)')]
    #     contigs = re.split(r'^DEFINITION.+$', content, flags=re.MULTILINE)
    #     contigs = [re.sub('^FEATURES.+$', '', contig) for contig in contigs]

    #     for contig_id, contig in zip(contig_ids, contig):
    #         coords = 
        
            
    def to_df(self, prodigal_output:bool=True) -> pd.DataFrame:

        df = []
        for id_, seq, description in zip(self.ids, self.seqs, self.descriptions):
            row = {'description':description} 
            row['id'] = id_
            row['seq'] = seq
            df.append(row)
        df = pd.DataFrame(df).set_index('id')
        return df

    def write(self, path:str, mode:str='w') -> NoReturn:
        f = open(path, mode=mode)
        records = []
        for id_, seq, description in zip(self.ids, self.seqs, self.descriptions):
            record = SeqRecord(Seq(seq), id=id_, description=description)
            records.append(record)
        SeqIO.write(records, f, 'fasta')
        f.close()