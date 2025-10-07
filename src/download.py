import os 
import pandas as pd 
import subprocess
from tqdm import tqdm 
import shutil
import numpy as np 
import io 
import requests 
import re 
import glob
import zipfile 

class NCBI():
    cleanup_files = ['README.md', 'md5sum.txt', 'ncbi.zip']
    cleanup_dirs = ['ncbi_dataset']

    src_dir = 'ncbi_dataset/data'

    def __init__(self):
        pass

    def get_genomes(self, genome_ids:list, include:list=['protein', 'genome'], data_dir:str='../data/ncbi/'):

        pbar = tqdm(genome_ids)
        for genome_id in pbar:
            pbar.set_description(f'NCBI.get_genomes: Downloading data for {genome_id}.')

            output_path = os.path.join(data_dir, f'{genome_id}.zip')
            if os.path.exists(output_path):
                continue

            cmd = f"datasets download genome accession {genome_id} --filename {output_path} --include {','.join(include)} --no-progressbar"
            try:
                subprocess.run(cmd, shell=True, check=True, stdout=subprocess.DEVNULL)
            except Exception as err:
                print(f'NCBI.get_genomes: Failed to download data for {genome_id}. Returned error message "{err}"')

    
def extract(data_dir:str='../data/ncbi'):
    protein_file_pattern = r'faa'
    genome_file_pattern = r'fna'

    def extract_file(zf, member:str, output_path:str=None):
        if not os.path.exists(output_path):
            with zf.open(member, 'r') as src, open(output_path, 'wb') as dst:
                dst.write(src.read())

    for path in tqdm(glob.glob(os.path.join(data_dir, '*.zip')), desc=f'extract: Extracting zipfiles in {data_dir}'):
        genome_id = os.path.basename(path).replace('.zip', '')
        
        try:
            zf = zipfile.ZipFile(path, 'r')
            for member in zf.namelist():
                if re.search(protein_file_pattern, str(member)) is not None:
                    extract_file(zf, member, output_path=os.path.join(data_dir, 'proteins', f'{genome_id}.fa'))
                if re.search(genome_file_pattern, str(member)) is not None:
                    extract_file(zf, member, output_path=os.path.join(data_dir, 'genomes', f'{genome_id}.fn'))
            zf.close()
        except:
            print(f'extract: Could not extract file at {path}')