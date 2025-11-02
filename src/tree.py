import pandas as pd 
import numpy as np 
import os 
from Bio import Phylo
import io 
import subprocess 
from tqdm import tqdm 
import glob 
import matplotlib.cm as cm
import matplotlib.colors
import numpy as np
import pandas as pd 
import re 


def tree_relabel(input_path, output_path, label_map:dict=None):

    assert len(list(label_map.keys())) == len(set(label_map.keys())), 'tree_relabel: New label IDs are not unique.'

    with open(input_path, 'r') as f:
        tree = f.read()
    for old_label, new_label in label_map.items():
        tree = tree.replace(old_label, new_label)
    with open(output_path, 'w') as f:
        f.write(tree)


def tree_remove_nonterminals(tree):

    for clade in tree.find_clades():
        if not clade.is_terminal():
            clade.name = None
    return tree 


def tree_subset(tree, ids:list=[]):

    nodes = [tree.find_any(id_) for id_ in ids]
    tree = tree.common_ancestor(nodes)

    for leaf in tree.get_terminals():
        if leaf.name not in ids:
            tree.prune(leaf)

    tree = tree_remove_nonterminals(tree)
    return tree 


def tree_write(tree, path:str=None):
    Phylo.write(tree, path, format='newick')


def tree_write(tree, path:str, fmt:str='newick'):
    Phylo.write(tree, path, format=fmt)


def tree_get_palette_continuous(name:str='coolwarm', values=None):
    cmap = cm.get_cmap(name)
    # Want to normalize between 0 and 1. Use this to get the color string, but preserve the original value for mapping. 
    normalized_values = (np.array(values) + min(values)) / max(values)
    palette = {value:matplotlib.colors.to_hex(cmap(normalized_value)) for value, normalized_value in zip(values, normalized_values)}
    return palette

def tree_get_palette_discrete(name:str='coolwarm', values=None):
     
    n = len(np.unique(values)) # Number of discrete colors you want
    cmap = cm.get_cmap(name, n) # Get the continuous coolwarm colormap
    cmap = matplotlib.colors.ListedColormap(cmap(np.arange(n))) # Convert to a ListedColormap for better control
    palette = {value: matplotlib.colors.to_hex(cmap(i)) for i, value in enumerate(values)}
    return palette


def _tree_get_legend(palette:dict, title:str=''):

    labels = list(palette.keys())
    lines = list()
    lines += [f'LEGEND_TITLE,{title}']
    lines += ['LEGEND_SHAPES,' + ','.join(['1'] * len(labels))]
    lines += ['LEGEND_COLORS,' + ','.join([palette[label] for label in labels])]
    lines += ['LEGEND_LABELS,' + ','.join(labels)]
    return lines 


def tree_make_annotation_file(df:pd.DataFrame, palette=None, path:str=None, field:str=None, styles:dict=dict(), sizes=dict(), legend:bool=True, dataset_label:str='dataset'):
    lines = ['DATASET_STYLE']
    lines += ['SEPARATOR COMMA'] # Cannot use comma for this dataset. 
    lines += [f'DATASET_LABEL,{dataset_label}']
    lines += ['COLOR,#000000']
    if legend:
        lines += _tree_get_legend(palette, title=field)
    # 2097|1502,label,clade,#ff0000,1,normal
    lines += ['DATA']

    for row in df.itertuples():
        color = palette[getattr(row, str(field))]
        style = styles.get(row.Index, 'normal') 
        size = sizes.get(row.Index, 1)
        lines += [f'{row.Index},label,node,{color},{size},{style}']
        # lines += [f'{row.Index},label,{color},normal,{width}']
    with open(path, 'w') as f:
        f.write('\n'.join(lines))


# def tree_make_annotation_file(df:pd.DataFrame, palette=None, path:str=None, field:str=None, styles:dict=dict(), sizes=dict(), legend:bool=True):
#     lines = ['TREE_COLORS', 'SEPARATOR COMMA']     # NODE_ID TYPE COLOR LABEL_OR_STYLE SIZE_FACTOR
#     if legend:
#         lines += _tree_get_legend(palette, title=field)
#     lines += ['DATA']
#     for row in df.itertuples():
#         color = palette[getattr(row, str(field))]
#         style = styles.get(row.Index, 'normal') 
#         size = sizes.get(row.Index, 1)
#         lines += [f'{row.Index},label,{color},{style},{size}']
#         # lines += [f'{row.Index},label,{color},normal,{width}']
#     with open(path, 'w') as f:
#         f.write('\n'.join(lines))





# def make_itol_connection(connections, path:str=None, dataset_label:str='dataset'):

#     lines = ['DATASET_CONNECTION']
#     lines += ['SEPARATOR COMMA'] # Cannot use comma for this dataset. 
#     lines += [f'DATASET_LABEL,{dataset_label}']
#     lines += ['COLOR,#00ff00']
#     lines += ['CURVE_ANGLE,60']
#     lines += ['DRAW_ARROWS,0']
#     lines += ['MAXIMUM_LINE_WIDTH,5']
#     lines += ['ALIGN_TO_LABELS,1']
#     lines += ['DATA']

#     for (node_1, node_2) in connections:
#         lines += [f'{node_1},{node_2},1,{gray},normal']
#     with open(path, 'w') as f:
#         f.write('\n'.join(lines))


# def make_itol_tanglegram(connections:list, tree_path='../data/arf1_cleaned_rep_seq.tree', dataset_label:str='dataset', path:str=None):

#     lines = ['DATASET_TANGLEGRAM']
#     lines += ['SEPARATOR TAB'] # Cannot use comma for this dataset. 
#     lines += [f'DATASET_LABEL\t{dataset_label}']
#     lines += ['COLOR\t#00ff00']
#     lines += ['TANGLEGRAM_TREE']

#     with open(tree_path, 'r') as f:
#         tree = f.read().strip()
#         lines += [tree]
#     lines += ['END_TANGLEGRAM_TREE']
#     lines += ['DEFAULT_CONNECTION_STYLE\tnormal']
#     lines += ['DEFAULT_CONNECTION_COLOR\t#ffffff'] # Make default color white. 
#     lines += ['DEFAULT_CONNECTION_WIDTH\t1']
#     lines += ['DATA']
#     for (node_1, node_2) in connections:
#         lines += [f'connect\t{node_1}\t{node_2}\t1\t{gray}']
#     with open(path, 'w') as f:
#         f.write('\n'.join(lines))