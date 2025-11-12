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
from itertools import combinations

def tree_get_distances(tree, query_ids=None, target_ids=None):
    query_leaves = [tree.find_any(genome_id) for genome_id in query_ids]
    if target_ids is not None:
        target_leaves = [tree.find_any(genome_id) for genome_id in target_ids]
        # leaf_pairs = [(query_leaf, target_leaf) for query_leaf in query_leaves for target_leaf in target_leaves]
        leaf_pairs = list(zip(query_leaves, target_leaves)) # Compute pairwise distances. 
        leaves = query_leaves + target_leaves
    else:
        leaves = query_leaves
        leaf_pairs = list(combinations(query_leaves, 2))
    print(f'tree_get_distances: Computing distances between {len(leaf_pairs)} pairs of leaves.')

    # Precompute the distances to the roots. 
    distances_to_root = {leaf:tree.distance(tree.root, leaf) for leaf in tqdm(leaves, desc='tree_get_distances: Precomputing distances to leaves.')}
    distance_df = list()
    
    for query_leaf, target_leaf in tqdm(leaf_pairs, desc='tree_get_distances'): # Using combinations, so there should not be any duplicates or zeros. 
    # for query_leaf, target_leaf in list(combinations(leaves, 2)): # Using combinations, so there should not be any duplicates or zeros. 
        row = {'query_genome_id':query_leaf.name, 'target_genome_id':target_leaf.name}

        lca = tree.common_ancestor(query_leaf, target_leaf) # Using the LCA to get distance is *much* faster (thanks ChatGPT)... should look into why. 
        lca_distance_to_root = distances_to_root.get(lca, tree.distance(tree.root, lca)) # Only compute the distance if it's not already in the dictionary. 
        distances_to_root[lca] = lca_distance_to_root # Store in the lookup table. 
        row['distance'] = distances_to_root[query_leaf] + distances_to_root[target_leaf] - (2 * lca_distance_to_root) 
        distance_df.append(row) # Only add if the distance is not too far. 

    distance_df = pd.DataFrame(distance_df)
    return distance_df

# def _tree_get_distances_preprocess(df:pd.DataFrame):
#     df = df.copy()
#     df['protein_id'] = df.index 
#     df = df.set_index('genome_id') # Need to use the genome ID as an index to work with the tree. 
#     return df

# def tree_get_distances(query_df, tree, target_df=None, max_distance=2.5):

#     query_df = _tree_get_distances_preprocess(query_df)
#     query_leaves = [tree.find_any(genome_id) for genome_id in query_df.index]

#     if target_df is not None:
#         target_df = _tree_get_distances_preprocess(target_df)
#         target_leaves = [tree.find_any(genome_id) for genome_id in target_df.index]
#         leaf_pairs = [(query_leaf, target_leaf) for query_leaf in query_leaves for target_leaf in target_leaves]
#         leaves = query_leaves + target_leaves
#     else:
#         target_df = query_df 
#         leaves = query_leaves
#         leaf_pairs = list(combinations(query_leaves, 2))

#     # Precompute the distances to the roots. 
#     distances_to_root = {leaf:tree.distance(tree.root, leaf) for leaf in leaves} # Can use the clade objects themselves as keys. 
#     distance_df = list()
    
#     for query_leaf, target_leaf in tqdm(leaf_pairs, desc='tree_get_distances'): # Using combinations, so there should not be any duplicates or zeros. 
#     # for query_leaf, target_leaf in list(combinations(leaves, 2)): # Using combinations, so there should not be any duplicates or zeros. 
#         query_genome_id, target_genome_id = query_leaf.name, target_leaf.name 
#         row = {'query_genome_id':query_leaf.name, 'target_genome_id':target_leaf.name}

#         lca = tree.common_ancestor(query_leaf, target_leaf) # Using the LCA to get distance is *much* faster (thanks ChatGPT)... should look into why. 
#         lca_distance_to_root = distances_to_root.get(lca, tree.distance(tree.root, lca)) # Only compute the distance if it's not already in the dictionary. 
#         distances_to_root[lca] = lca_distance_to_root # Store in the lookup table. 
#         distance = distances_to_root[query_leaf] + distances_to_root[target_leaf] - (2 * lca_distance_to_root)
        
#         if (max_distance is None) or (distance < max_distance):
#             row['distance'] = distance 
#             row['query_category'], row['target_category'] = query_df.loc[query_genome_id].category, target_df.loc[target_genome_id].category
#             row['query_protein_id'], row['target_protein_id'] = query_df.loc[query_genome_id].protein_id, target_df.loc[target_genome_id].protein_id
#             row['category'] = ' vs. '.join(sorted([row['query_category'], row['target_category']]))
#             distance_df.append(row) # Only add if the distance is not too far. 

#     distance_df = pd.DataFrame(distance_df)
#     return distance_df


def tree_relabel_file(input_path, output_path, label_map:dict=None):

    assert len(list(label_map.keys())) == len(set(label_map.keys())), 'tree_relabel: New label IDs are not unique.'

    with open(input_path, 'r') as f:
        tree = f.read()
    for old_label, new_label in label_map.items():
        tree = tree.replace(old_label, new_label)
    with open(output_path, 'w') as f:
        f.write(tree)


def tree_relabel(tree, label_map:dict=None):

    assert len(list(label_map.keys())) == len(set(label_map.keys())), 'tree_relabel: New label IDs are not unique.'
    for node in tree.get_terminals():
        node.name = label_map.get(node.name, node.name)
    return tree 

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