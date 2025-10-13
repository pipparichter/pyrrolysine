import pandas as pd 
from src.files.fasta import FASTAFile
from tqdm import tqdm
import numpy as np 
import os
import re 

genus_to_order_map = dict()
genus_to_order_map['Methanosarcina'] = 'Methanosarcinales'
genus_to_order_map['Methanolobus'] = 'Methanosarcinales'
genus_to_order_map['Methanohalobium'] = 'Methanosarcinales'
genus_to_order_map['Methanosalsum'] = 'Methanosarcinales'
genus_to_order_map['Methanomethylophilus'] = 'Methanomassiliicoccales'
genus_to_order_map['Methanohalophilus'] = 'Methanosarcinales'
genus_to_order_map['Methanococcoides'] = 'Methanosarcinales'
genus_to_order_map['Methanomethylovorans'] = 'Methanosarcinales'
genus_to_order_map['Methanomassiliicoccus_A'] = 'Methanomassiliicoccales'
genus_to_order_map['Methanomassiliicoccus'] = 'Methanomassiliicoccales'
genus_to_order_map['Methanoprimaticola'] = 'Methanomassiliicoccales'
genus_to_order_map['Methanoplasma'] = 'Methanomassiliicoccales'
genus_to_order_map['Methermicoccus'] = 'Methanosarcinales'
genus_to_order_map['Hecatella'] = 'Hecatellales'
genus_to_order_map['Methanosalsum'] = 'Methanosarcinales'
genus_to_order_map['none'] = 'none'


get_genus = lambda taxonomy : re.search('g__([^;]+)', taxonomy).group(1) if (re.search('g__([^;]+)', taxonomy) is not None) else 'none'
get_species = lambda taxonomy : re.search('s__([^;]+)', taxonomy).group(1) if (re.search('s__([^;]+)', taxonomy) is not None) else 'none'
get_order = lambda taxonomy : re.search('o__([^;]+)', taxonomy).group(1) if (re.search('o__([^;]+)', taxonomy) is not None) else 'none'
get_phylum = lambda taxonomy : re.search('p__([^;]+)', taxonomy).group(1) if (re.search('p__([^;]+)', taxonomy) is not None) else 'none'

def load_kivenson_table_2_metadata(data_dir:str='../data'):


    taxa = ['Methanococcoides', 'Methanosarcina', 'Methanohalobium', 'Methanohalophilus', 'Methanimicrococcus', 'Methermicoccus']
    taxa += ['Methanolobus', 'Methanomethylovorans', 'Methanosalsum', 'Methanomicrobia', 'Methanonatronarchaeum', 'Methanohalarchaeum']
    taxa += ['MSBL_1_clade', 'Methanoplasma', 'Methanomethylophilus', 'Methanomassiliicoccus', 'Methanogranum', 'Thermoplasmatota']
    taxa += ['Hydrothermarchaeum', 'Borrarchaeum', 'Korarchaeota', 'Methylarchaceae', 'Hecatella', 'Bathyarchaeota']
    taxa = '|'.join(taxa)

    # gtdb_genome_metadata_df = pd.read_csv(f'{data_dir}/ar53_metadata_r226.tsv', sep='\t')
    # gtdb_genome_metadata_df.accession = [genome_id.replace('RS_','').replace('GB_', '') for genome_id in gtdb_genome_metadata_df.accession]
    # gtdb_genome_metadata_df = gtdb_genome_metadata_df.set_index('accession')

    # table_2_genome_ids = pd.read_csv(f'{data_dir}/kivenson_2023_table_2.csv', usecols=['Assembly Accession'])['Assembly Accession'].tolist()
    
    # print(f'load_kivenson_table_2_metadata: Loaded {len(table_2_genome_ids)} from Kivenson Table 2.')
    # table_2_genome_ids = np.intersect1d(table_2_genome_ids, gtdb_genome_metadata_df.index)
    # print(f'load_kivenson_table_2_metadata: {len(table_2_genome_ids)} from Kivenson Table 2 with GTDB metadata.')

    # table_2_metadata_df = gtdb_genome_metadata_df.loc[table_2_genome_ids]
    # table_2_metadata_df['order'] = table_2_metadata_df.gtdb_taxonomy.apply(get_order)
    # table_2_metadata_df['order'] = table_2_metadata_df.gtdb_taxonomy.apply(get_order)

    table_2_metadata_df = pd.read_csv(f'{data_dir}/kivenson_2023_table_2.csv', usecols=['Assembly Accession', 'Organism Name'])
    table_2_metadata_df = table_2_metadata_df.rename(columns={'Assembly Accession':'accession', 'Organism Name':'organism_name'})
    table_2_metadata_df = table_2_metadata_df.set_index('accession')
    table_2_metadata_df['genus'] = [re.search(taxa, name).group(0) if (re.search(taxa, name) is not None) else 'none' for name in table_2_metadata_df.organism_name]
    table_2_metadata_df['order'] = table_2_metadata_df.genus.map(genus_to_order_map)
    return table_2_metadata_df


def load_kivenson_table_5_metadata(data_dir:str='../data'):
    # For some reason, there are genome IDs in SI Table 5 which are not in the SI Table 2 list, so had to add them seperately. 
    table_5_metadata_df = pd.read_csv(f'{data_dir}/kivenson_2023_table_5.csv')
    table_5_metadata_df['accession'] = ['_'.join(file_name.split('_')[:2]) for file_name in table_5_metadata_df['Filename']]
    table_5_metadata_df = table_5_metadata_df.set_index('accession')
    table_5_metadata_df = table_5_metadata_df.drop(columns=['Euryarchaeota / Halobacteriota'])
    columns = {'Genus / taxonomy (from NCBI)':'ncbi_taxonomy', 'GTDB taxonomy':'gtdb_taxonomy', 'Genome':'genome', '%TAG':'tag_percent', 'Filename':'filename'}
    columns.update({f'Category {i}':f'category {i}' for i in range(1, 5)})
    table_5_metadata_df = table_5_metadata_df.rename(columns=columns)
    table_5_metadata_df['order'] = table_5_metadata_df.gtdb_taxonomy.apply(get_order)
    table_5_metadata_df['genus'] = table_5_metadata_df.gtdb_taxonomy.apply(get_genus)
    return table_5_metadata_df


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