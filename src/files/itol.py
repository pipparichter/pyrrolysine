import os 
import subprocess 
from tqdm import tqdm 
import glob 
import matplotlib.cm as cm
import matplotlib.colors
import numpy as np
import pandas as pd 
import re 


darkblue = '#3182bd'
lightblue = '#6baed6'
gray = '#969696'
lightgreen = '#74c476'
darkgreen = '#31a354'
red = '#e6550d'
orange = '#fdae6b'


def get_palette(name=None, values=None):
    name = 'coolwarm' if (name is None) else name
    cmap = cm.get_cmap(name)
    # Want to normalize between 0 and 1. Use this to get the color string, but preserve the original value for mapping. 
    normalized_values = (np.array(values) + min(values)) / max(values)
    palette = {value:matplotlib.colors.to_hex(cmap(normalized_value)) for value, normalized_value in zip(values, normalized_values)}
    return palette


def make_itol_annotation(df:pd.DataFrame, palette='coolwarm', path:str=None, field:str=None, styles:dict=dict(), sizes=dict()):
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


def make_itol_connection(connections, path:str=None, dataset_label:str='dataset'):

    lines = ['DATASET_CONNECTION']
    lines += ['SEPARATOR COMMA'] # Cannot use comma for this dataset. 
    lines += [f'DATASET_LABEL,{dataset_label}']
    lines += ['COLOR,#00ff00']
    lines += ['CURVE_ANGLE,60']
    lines += ['DRAW_ARROWS,0']
    lines += ['MAXIMUM_LINE_WIDTH,5']
    lines += ['ALIGN_TO_LABELS,1']
    lines += ['DATA']

    for (node_1, node_2) in connections:
        lines += [f'{node_1},{node_2},1,{gray},normal']
    with open(path, 'w') as f:
        f.write('\n'.join(lines))


def make_itol_tanglegram(connections:list, tree_path='../data/arf1_cleaned_rep_seq.tree', dataset_label:str='dataset', path:str=None):

    lines = ['DATASET_TANGLEGRAM']
    lines += ['SEPARATOR TAB'] # Cannot use comma for this dataset. 
    lines += [f'DATASET_LABEL\t{dataset_label}']
    lines += ['COLOR\t#00ff00']
    lines += ['TANGLEGRAM_TREE']

    with open(tree_path, 'r') as f:
        tree = f.read().strip()
        lines += [tree]
    lines += ['END_TANGLEGRAM_TREE']
    lines += ['DEFAULT_CONNECTION_STYLE\tnormal']
    lines += ['DEFAULT_CONNECTION_COLOR\t#ffffff'] # Make default color white. 
    lines += ['DEFAULT_CONNECTION_WIDTH\t1']
    lines += ['DATA']
    for (node_1, node_2) in connections:
        lines += [f'connect\t{node_1}\t{node_2}\t1\t{gray}']
    with open(path, 'w') as f:
        f.write('\n'.join(lines))