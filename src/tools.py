import os 
import subprocess 
from tqdm import tqdm 
import glob 
import matplotlib.cm as cm
import matplotlib.colors
import numpy as np
import pandas as pd 
import re 

def get_palette(name=None, values=None):
    name = 'coolwarm' if (name is None) else name
    cmap = cm.get_cmap(name)
    # Want to normalize between 0 and 1. Use this to get the color string, but preserve the original value for mapping. 
    normalized_values = (np.array(values) + min(values)) / max(values)
    palette = {value:matplotlib.colors.to_hex(cmap(normalized_value)) for value, normalized_value in zip(values, normalized_values)}
    return palette


def make_itol_annotation_file(df:pd.DataFrame, palette='coolwarm', path:str=None, field:str=None, styles:dict=dict(), sizes=dict()):
    header_lines = ['TREE_COLORS', 'SEPARATOR COMMA', 'DATA']     # NODE_ID TYPE COLOR LABEL_OR_STYLE SIZE_FACTOR

    lines = list()

    palette = get_palette(name=palette, values=df[field].values) if (type(palette) != dict) else palette

    if type(palette) == str:
        values = df[field].unique()
        palette = get_palette(name=palette, values=values)

    for row in df.itertuples():
        color = palette[getattr(row, field)]
        style = styles.get(row.Index, 'normal') 
        size = sizes.get(row.Index, 1)
        lines += [f'{row.Index},label,{color},{style},{size}']
        # lines += [f'{row.Index},label,{color},normal,{width}']
    lines = header_lines + lines 
    with open(path, 'w') as f:
        f.write('\n'.join(lines))



# def make_itol_annotation_file(df:pd.DataFrame, palette='coolwarm', path:str=None, field:str=None, size:int=10):
#     header_lines = ['DATASET_SYMBOL', 'SEPARATOR SPACE', 'DATASET_LABEL leaves', 'COLOR #000000', 'DATA']

#     if type(palette) == str:
#         values = df[field].unique()
#         palette = get_palette(name=palette, values=values)

#     lines = list()
#     for row in df.itertuples():
#         color = palette[getattr(row, field)]
#         lines += [f'{row.Index} 2 {color} 1 {size} 0 {row.Index}']
#         # lines += [f'{row.Index} {color} 1 {size} #000000 1']
#     lines = header_lines + lines 
#     with open(path, 'w', encoding='utf-8') as f:
#         f.write('\n'.join(lines))


def run_prodigal(input_dir:str='../data/ncbi/genomes/', output_dir='../data/prodigal/'):
    for path in tqdm(glob.glob(os.path.join(input_dir, '*')), 'run_prodigal'):
        output_path = re.sub('.fn|.fna', '.fa', os.path.basename(path))
        output_path = os.path.join(output_dir, output_path)
        try:
            if not os.path.exists(output_path):
                subprocess.run(f'prodigal -i {path} -a {output_path}', shell=True, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        except:
            print(f'run_prodigal: Could not run Prodigal on {path}')

def run_hmmer(input_dir:str='../data/prodigal/', query_path='../data/hmms/query.hmm', output_dir:str='../data/hmmer', overwrite:bool=False):
    for path in tqdm(glob.glob(os.path.join(input_dir, '*')), 'run_hmmer'):
        output_path = re.sub('.fa|.faa', '.tab', os.path.basename(path))
        output_path = os.path.join(output_dir, output_path)
        if not os.path.exists(output_path) or overwrite:
            try:
                subprocess.run(f'hmmsearch --domtblout {output_path} {query_path} {path}', shell=True, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
            except:
                print(f'run_hmmer: HMM search failed for {path}')