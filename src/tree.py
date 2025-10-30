import pandas as pd 
import numpy as np 
import os 
from Bio import Phylo
import io 


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


def tree_load_ar53(path:str='../data/ar53.tree', genome_ids:list=None):
    with open(path, 'r') as f:
        tree = f.read()
    tree = tree.replace('RS_', '').replace('GB_', '') # Remove the prefixes to allow mapping. 
    tree = Phylo.read(io.StringIO(tree), format='newick')
    
    if genome_ids is not None:
        tree = tree_subset(tree, ids=genome_ids)

    return tree 


def tree_write(tree, path:str, fmt:str='newick'):
    Phylo.write(tree, path, format=fmt)