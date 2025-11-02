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


