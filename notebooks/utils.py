import matplotlib.pyplot as plt 
import pandas as pd
import numpy as np 
from tqdm import tqdm 
import glob
import os
from src.files.fasta import FASTAFile
import re 
import seaborn as sns
import src.files.fasta as fasta
import subprocess
from Bio import Phylo
import io
from src.tree import *


# ['#3182bd', '#6baed6', '#9ecae1', '#c6dbef',
# '#e6550d', '#fd8d3c', '#fdae6b', '#fdd0a2',
# '#31a354', '#74c476', '#a1d99b', '#c7e9c0',
# '#756bb1', '#9e9ac8', '#bcbddc', '#dadaeb',
# '#636363', '#969696', '#bdbdbd', '#d9d9d9']

def annotate_residues(annotations:dict, ax=None, lines_only:bool=False):
    for seq, (start, stop) in annotations.items():
        positions = np.arange(start, stop)
        for x, aa in zip(positions, list(seq)):
            if not lines_only:
                ax.text(x, ax.get_ylim()[-1], aa, ha='center', va='bottom')
            ax.axvline(x, ls='--', lw=0.5, color='gray')

GTDB_DATA_DIR = '../data/gtdb/'

darkblue = '#3182bd'
lightblue = '#6baed6'
gray = '#969696'
lightgreen = '#74c476'
darkgreen = '#31a354'
red = '#e6550d'
orange = '#fdae6b'
black= '#000000'


gtdb_get_genus = lambda taxonomy : re.search('g__([^;]+)', taxonomy).group(1) if (re.search('g__([^;]+)', taxonomy) is not None) else 'none'
gtdb_get_species = lambda taxonomy : re.search('s__([^;]+)', taxonomy).group(1) if (re.search('s__([^;]+)', taxonomy) is not None) else 'none'
gtdb_get_order = lambda taxonomy : re.search('o__([^;]+)', taxonomy).group(1) if (re.search('o__([^;]+)', taxonomy) is not None) else 'none'
gtdb_get_phylum = lambda taxonomy : re.search('p__([^;]+)', taxonomy).group(1) if (re.search('p__([^;]+)', taxonomy) is not None) else 'none'
gtdb_get_family = lambda taxonomy : re.search('f__([^;]+)', taxonomy).group(1) if (re.search('f__([^;]+)', taxonomy) is not None) else 'none'
gtdb_get_class = lambda taxonomy : re.search('c__([^;]+)', taxonomy).group(1) if (re.search('c__([^;]+)', taxonomy) is not None) else 'none'


def gtdb_load_ar53_metadata(data_dir:str=GTDB_DATA_DIR):
    gtdb_metadata_df = pd.read_csv(os.path.join(data_dir, 'ar53_metadata_r226.tsv'), sep='\t')
    gtdb_metadata_df.accession = [genome_id.replace('RS_','').replace('GB_', '') for genome_id in gtdb_metadata_df.accession]
    gtdb_metadata_df['genus'] = gtdb_metadata_df.gtdb_taxonomy.apply(gtdb_get_genus)
    gtdb_metadata_df['order'] = gtdb_metadata_df.gtdb_taxonomy.apply(gtdb_get_order)
    gtdb_metadata_df['species'] = gtdb_metadata_df.gtdb_taxonomy.apply(gtdb_get_species)
    gtdb_metadata_df['phylum'] = gtdb_metadata_df.gtdb_taxonomy.apply(gtdb_get_phylum)
    gtdb_metadata_df['family'] = gtdb_metadata_df.gtdb_taxonomy.apply(gtdb_get_family)
    gtdb_metadata_df['class'] = gtdb_metadata_df.gtdb_taxonomy.apply(gtdb_get_class)
    gtdb_metadata_df = gtdb_metadata_df.set_index('accession')
    return gtdb_metadata_df


def gtdb_load_ar53_tree(data_dir:str=GTDB_DATA_DIR, genome_ids:list=None):
    with open(os.path.join(data_dir, 'ar53.tree'), 'r') as f:
        tree = f.read()
    tree = tree.replace('RS_', '').replace('GB_', '') # Remove the prefixes to allow mapping. 
    tree = Phylo.read(io.StringIO(tree), format='newick')
    if genome_ids is not None:
        tree = tree_subset(tree, ids=genome_ids)
    return tree 


def load_arf1_dataset(path:str='../data/arf1_cleaned.csv', stop_codon_metadata_path:str='../data/arf1_stop_codon_metadata.csv', exclude_genome_ids=[]):

    gtdb_metadata_df = gtdb_load_ar53_metadata()
    gtdb_metadata_df['genome_id'] = gtdb_metadata_df.index

    arf1_df = pd.read_csv(path, index_col=0).drop(columns=['genus', 'order', 'Unnamed: 0'], errors='ignore')
    index = arf1_df.index
    arf1_df = arf1_df.merge(gtdb_metadata_df, left_on='genome_id', right_on='genome_id', how='left')
    arf1_df.index = index # Restore the index after the merge. 

    # Exclude the specified genome IDs, probably the Methanosarcina species that lost Pyl. 
    arf1_df = arf1_df[~arf1_df.genome_id.isin(exclude_genome_ids)].copy()

    # Add stop codon information to the data.
    stop_codon_metadata_df = pd.read_csv(stop_codon_metadata_path, index_col=0)
    arf1_df['tag_count'] = arf1_df.genome_id.map(stop_codon_metadata_df.groupby('genome_id').TAG.first())
    arf1_df['stop_codon_count'] = arf1_df.genome_id.map(stop_codon_metadata_df.groupby('genome_id')['total'].first())
    arf1_df['tag_percent'] = arf1_df.tag_count / arf1_df.stop_codon_count

    if 'has_pyl' in arf1_df.columns:
        # I think more granular categories could be helpful:
        # (1) Pyl+ and largely re-coded (TAG < 5%)
        # (2) Pyl+ which still use lots of TAG stops (TAG > 5%) 
        # (3) Pyl- (including the weird outliers

        masks = dict()
        masks['pyl+ recoded'] = (arf1_df.tag_percent < 0.05) & (arf1_df.has_pyl)
        masks['pyl+'] = (arf1_df.tag_percent >= 0.05) & (arf1_df.has_pyl)
        masks['pyl-'] = (~arf1_df.has_pyl)

        categories = list(masks.keys())
        arf1_df['category'] = np.select([masks[category] for category in categories], categories, default='none')

    return arf1_df


LEGEND_KWARGS = dict()
LEGEND_KWARGS['handletextpad'] = 0.5
LEGEND_KWARGS['handlelength'] = 1
LEGEND_KWARGS['labelspacing'] = 0.2
LEGEND_KWARGS['borderpad'] = 0.4
LEGEND_KWARGS['markerscale'] = 0.5
LEGEND_KWARGS['loc'] = 'lower right'
LEGEND_KWARGS['prop'] = {'size':6}

def shrink_legend_barplot(ax, **kwargs):
    legend = ax.get_legend()
    handles, labels = legend.get_patches(), [label.get_text() for label in legend.get_texts()] 
    kwargs_ = LEGEND_KWARGS.copy()
    kwargs_.update(kwargs)
    ax.legend(handles, labels, title='', **kwargs_)


def shrink_legend_scatterplot(ax, **kwargs):
    handles, labels = ax.get_legend_handles_labels()
    kwargs_ = LEGEND_KWARGS.copy()
    kwargs_.update(kwargs)
    ax.legend(handles, labels, title='', **kwargs_)



def run_muscle(input_path:str):
    muscle_output_path = re.sub(r'\.fa.*', '.afa', input_path) # input_path.replace('.fa', '.afa').replace('.fasta', '.afa')
    trimal_output_path = re.sub(r'\.fa.*', '_trimmed.afa', input_path)
    if not os.path.exists(muscle_output_path):
        print(f'run_muscle: Generating file {muscle_output_path}')
        subprocess.run(f'~/muscle5.1.linux_intel64 -align {input_path} -output {muscle_output_path}', shell=True, check=True)
    subprocess.run(f'trimal -in {muscle_output_path} -out {trimal_output_path}', shell=True, check=True)
    return trimal_output_path

    # if build_tree:
    #     tree_path = muscle_output_path.replace('.afa', '')
    #     # subprocess.run(f'fasttree {muscle_output_path} > {tree_path}', shell=True, check=True)
    #     # Use the trimmed output!
    #     subprocess.run(f'iqtree -s {trimal_output_path} -m MFP -bb 1000 -alrt 1000 -nt AUTO -pre {tree_path}', shell=True, check=True)


def plot_scores_1d(scores:np.ndarray, start=0, stop=150, y_label:str='', ax:plt.Axes=None, color:str='steelblue'):
    scores = scores[start:stop]
    sns.lineplot(x=np.arange(len(scores)), y=scores, ax=ax, color=color)
    ax.set_xticks(np.arange(len(scores)), labels=np.arange(start, stop), fontsize='x-small')
    ax.set_xlim(xmin=start, xmax=stop)
    ax.set_ylabel(y_label)

# def plot_scores_2d(scores:np.ndarray, start=0, stop=150, ax:plt.Axes=None):
#     cmap = matplotlib.colors.LinearSegmentedColormap.from_list('palette', ['white', 'steelblue'])
#     plt.imshow(scores[start:stop, start:stop], cmap=cmap)
#     ax.set_yticks(np.arange(stop - start), labels=np.arange(start, stop), fontsize='x-small')
#     ax.set_xticks(np.arange(start, stop), labels=np.arange(start, stop), fontsize='x-small')


def parse_prodigal_description(df:pd.DataFrame):
    n = len(df)
    assert 'description' in df.columns, 'parse_prodigal_description: No description column to parse.'
    df_ = pd.DataFrame([fasta._parse_prodigal_description(description) for description in df.description], index=df.index)
    df_['start'] = df_.start.astype(int)
    df_['stop'] = df_.stop.astype(int)
    df_['strand'] = df_.strand.astype(int)
    df_['contig_id'] = [fasta.get_contig_id(id_) for id_ in df_.index]
    df = df.drop(columns=df_.columns, errors='ignore').merge(df_, left_index=True, right_index=True)
    assert len(df) == n, 'parse_prodigal_description: Something went wrong while merging the DataFrames.'
    return df


def plot_residue_counts(figure_df:pd.DataFrame, position:int=None, palette=None, ax=None, hue:str=None, stat:str='probability'):
    ax_df = figure_df[[hue, position]].copy()
    sns.histplot(ax_df, x=position, hue=hue, multiple='dodge', ax=ax, palette=palette, stat=stat, common_norm=False, shrink=0.8)
    ax.set_ylabel(stat)
    ax.get_legend().set_title('')
    ax.set_xlabel('')
    ax.set_title(f'position {position}')


def get_domain_boundaries(seq:str):
    patterns = dict()
    patterns['NIKS'] = 'NIKS'
    patterns['YxCxxxF'] = 'Y.C'
    # patterns['GGQ'] = 'GGQ'
    patterns['GTS'] = 'GR.'
    patterns['GGQ'] = 'GGQ'

    boundaries = dict()
    for domain, pattern in patterns.items():
        start = re.search(pattern, seq).start()
        stop =  start + len(domain)
        print(f'get_domain_boundaries: {domain} boundaries {start}-{stop}.')
        boundaries[domain] = (start, stop)
    return boundaries 


dayhoff = {'A':'A','G':'A','P':'A','S':'A','T':'A','D':'B','E':'B','N':'B','Q':'B','R':'C','H':'C','K':'C', 'M':'D','I':'D','L':'D','V':'D','F':'E','W':'E','Y':'E','C':'F' }

hp = {'A':'H','V':'H','L':'H','I':'H','P':'H','F':'H','W':'H','M':'H','Y':'H','G':'P','S':'P','T':'P','C':'P','N':'P','Q':'P','D':'P','E':'P','K':'P','R':'P','H':'P'}

ba = {'K':'B','R':'B','H':'B','D':'A','E':'A','A':'N','C':'N','F':'N','G':'N','I':'N','L':'N','M':'N','N':'N','P':'N','Q':'N','S':'N','T':'N','V':'N','W':'N','Y':'N'}

dayhoff_descriptions = dict()
dayhoff_descriptions['A'] = 'small hydrophilic'
dayhoff_descriptions['B'] = 'acidic amide'
dayhoff_descriptions['C'] = 'basic'
dayhoff_descriptions['D'] = 'hydrophobic aliphatic'
dayhoff_descriptions['E'] = 'aromatic'
dayhoff_descriptions['F'] = 'sulfur containing'

def load_msa(path, ids:list=None, conservation_threshold:float=0.8):
    is_conserved = lambda col : (col != '-').astype(int).mean() > conservation_threshold # Flag conserved positions as those where at least 80 percent of the sequences do not have a gap. 

    alignment_df = FASTAFile().from_fasta(path).to_df()
    if ids is not None:
        alignment_df = alignment_df.loc[ids].copy()
    alignment = np.array([list(seq) for seq in alignment_df.seq])
    positions = np.where([is_conserved(col) for col in alignment.T])[0]
    alignment = alignment[:, positions].copy() # Ignore the extra junk in the alignment.
    return pd.DataFrame(alignment, index=ids)



def make_chimerax_attribute_file(positions:list, scores:list, path:str='../data/alphafold/arf1_alphafold_attributes.defattr'):
    content = ['attribute: score\nmatch mode: 1-to-1\nrecipient: residues\n']
    for i, score in zip(positions, scores):
        content += [f'\t:{i + 1}\t{score}\n'] # Residues are 1-indexed. 

    with open(path, 'w') as f:
        f.write(''.join(content))


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
def hmmer_load(data_dir='../data/hmmer', max_e_values:dict=dict(), genome_ids:list=None, query_names:list=['aRF1_eRF1', 'pelota', 'pyrrolys_PylC', 'pyrrolys_PylB', 'pyrrolys_PylD', 'PylS_Nterm', 'PylS_Cterm'], best_hit_only:bool=True):
    hmmer_df = list()

    paths = np.array(glob.glob(os.path.join(data_dir, '*')))
    
    if genome_ids is not None:
        path_genome_ids = np.array([os.path.basename(path).replace('.tab', '') for path in paths])
        paths = paths[np.isin(path_genome_ids, genome_ids)].copy()

    for path in tqdm(paths, desc='hmmer_load'):
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




# genus_to_order_map = dict()
# genus_to_order_map['Methanosarcina'] = 'Methanosarcinales'
# genus_to_order_map['Methanolobus'] = 'Methanosarcinales'
# genus_to_order_map['Methanohalobium'] = 'Methanosarcinales'
# genus_to_order_map['Methanosalsum'] = 'Methanosarcinales'
# genus_to_order_map['Methanomethylophilus'] = 'Methanomassiliicoccales'
# genus_to_order_map['Methanohalophilus'] = 'Methanosarcinales'
# genus_to_order_map['Methanococcoides'] = 'Methanosarcinales'
# genus_to_order_map['Methanomethylovorans'] = 'Methanosarcinales'
# genus_to_order_map['Methanomassiliicoccus_A'] = 'Methanomassiliicoccales'
# genus_to_order_map['Methanomassiliicoccus'] = 'Methanomassiliicoccales'
# genus_to_order_map['Methanoprimaticola'] = 'Methanomassiliicoccales'
# genus_to_order_map['Methanoplasma'] = 'Methanomassiliicoccales'
# genus_to_order_map['Methermicoccus'] = 'Methanosarcinales'
# genus_to_order_map['Hecatella'] = 'Hecatellales'
# genus_to_order_map['Methanosalsum'] = 'Methanosarcinales'
# genus_to_order_map['none'] = 'none'


# def load_kivenson_table_2_metadata(data_dir:str='../data'):
#     taxa = ['Methanococcoides', 'Methanosarcina', 'Methanohalobium', 'Methanohalophilus', 'Methanimicrococcus', 'Methermicoccus']
#     taxa += ['Methanolobus', 'Methanomethylovorans', 'Methanosalsum', 'Methanomicrobia', 'Methanonatronarchaeum', 'Methanohalarchaeum']
#     taxa += ['MSBL_1_clade', 'Methanoplasma', 'Methanomethylophilus', 'Methanomassiliicoccus', 'Methanogranum', 'Thermoplasmatota']
#     taxa += ['Hydrothermarchaeum', 'Borrarchaeum', 'Korarchaeota', 'Methylarchaceae', 'Hecatella', 'Bathyarchaeota']
#     taxa = '|'.join(taxa)
#     table_2_metadata_df = pd.read_csv(f'{data_dir}/kivenson_2023_table_2.csv', usecols=['Assembly Accession', 'Organism Name'])
#     table_2_metadata_df = table_2_metadata_df.rename(columns={'Assembly Accession':'accession', 'Organism Name':'organism_name'})
#     table_2_metadata_df = table_2_metadata_df.set_index('accession')
#     table_2_metadata_df['genus'] = [re.search(taxa, name).group(0) if (re.search(taxa, name) is not None) else 'none' for name in table_2_metadata_df.organism_name]
#     table_2_metadata_df['order'] = table_2_metadata_df.genus.map(genus_to_order_map)
#     return table_2_metadata_df


# def load_kivenson_table_5_metadata(data_dir:str='../data'):
#     # For some reason, there are genome IDs in SI Table 5 which are not in the SI Table 2 list, so had to add them seperately. 
#     table_5_metadata_df = pd.read_csv(f'{data_dir}/kivenson_2023_table_5.csv')
#     table_5_metadata_df['accession'] = ['_'.join(file_name.split('_')[:2]) for file_name in table_5_metadata_df['Filename']]
#     table_5_metadata_df = table_5_metadata_df.set_index('accession')
#     table_5_metadata_df = table_5_metadata_df.drop(columns=['Euryarchaeota / Halobacteriota'])
#     columns = {'Genus / taxonomy (from NCBI)':'ncbi_taxonomy', 'GTDB taxonomy':'gtdb_taxonomy', 'Genome':'genome', '%TAG':'tag_percent', 'Filename':'filename'}
#     columns.update({f'Category {i}':f'category {i}' for i in range(1, 5)})
#     table_5_metadata_df = table_5_metadata_df.rename(columns=columns)
#     table_5_metadata_df['order'] = table_5_metadata_df.gtdb_taxonomy.apply(get_order)
#     table_5_metadata_df['genus'] = table_5_metadata_df.gtdb_taxonomy.apply(get_genus)
#     return table_5_metadata_df