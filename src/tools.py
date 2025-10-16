import os 
import subprocess 
from tqdm import tqdm 
import glob 
import matplotlib.cm as cm
import matplotlib.colors
import numpy as np
import pandas as pd 


def get_palette(name=None, values=None):
    name = 'coolwarm' if (name is None) else name
    cmap = cm.get_cmap(name)
    # Want to normalize between 0 and 1. Use this to get the color string, but preserve the original value for mapping. 
    normalized_values = (np.array(values) + min(values)) / max(values)
    palette = {value:matplotlib.colors.to_hex(cmap(normalized_value)) for value, normalized_value in zip(values, normalized_values)}
    return palette


def make_itol_annotation_file(df:pd.DataFrame, palette='coolwarm', path:str=None, field:str=None, widths:dict=dict()):
    header_lines = ['TREE_COLORS', 'SEPARATOR COMMA', 'DATA']     # NODE_ID TYPE COLOR LABEL_OR_STYLE SIZE_FACTOR

    lines = list()

    palette = get_palette(name=palette, values=df[field].values) if (type(palette) != dict) else palette

    if type(palette) == str:
        values = df[field].unique()
        palette = get_palette(name=palette, values=values)

    for row in df.itertuples():
        color = palette[getattr(row, field)]

        width = widths.get(getattr(row, field), 2)
        lines += [f'{row.Index},branch,{color},normal,{width}']
    lines = header_lines + lines 
    with open(path, 'w') as f:
        f.write('\n'.join(lines))


def run_prodigal(input_dir:str='../data/ncbi/genomes/', output_dir='../data/prodigal/'):
    for path in tqdm(glob.glob(os.path.join(input_dir, '*')), 'run_prodigal'):
        output_path = re.sub('.fn|.fna', '.fa', os.path.basename(path))
        output_path = os.path.join(output_dir, output_path)
        if not os.path.exists(output_path):
            subprocess.run(f'prodigal -i {path} -a {output_path}', shell=True, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

def run_hmmer(input_dir:str='../data/prodigal/', query_path='../data/query.hmm', output_dir:str='../data/hmmer'):
    for path in tqdm(glob.glob(os.path.join(input_dir, '*')), 'run_hmmer'):
        output_path = re.sub('.fa|.faa', '.fa', os.path.basename(path))
        output_path = os.path.join(output_dir, output_path)
        if not os.path.exists(output_path):
            try:
                subprocess.run(f'hmmsearch --domtblout {output_path} {query_path} {path}', shell=True, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
            except:
                print(f'run_hmmer: HMM search failed for {path}')