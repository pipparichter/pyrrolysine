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
from Bio import Entrez 

Entrez.email = 'prichter@berkeley.edu'


def download_genomes(genome_ids:list, include:list=['protein', 'genome'], data_dir:str='../data/ncbi/'):

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

    
def extract_genomes(data_dir:str='../data/ncbi'):
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

patterns = dict()
patterns['id'] = r'<Id>(\d+)</Id>'
patterns['coordinates'] = r'<GBQualifier_name>coded_by</GBQualifier_name>\s+<GBQualifier_value>(.+)</GBQualifier_value>' # From a protein fetch result. 
patterns['taxonomy'] = r'<GBSeq_taxonomy>(.+)</GBSeq_taxonomy>' # From a nuccore fetch result.
patterns['organism'] = r'<GBSeq_organism>(.+)</GBSeq_organism>' # From a nuccore fetch result. 
patterns['nuccore_id'] = r'<GBSeq_source-db>accession (.+)</GBSeq_source-db' # From a protein fetch result. 
patterns['assembly_id'] = r'<AssemblyAccession>(.+)</AssemblyAccession>' # From a protein fetch result. 
patterns['seq'] = r'<GBSeq_sequence>(.+)</GBSeq_sequence>'

def parse_entrez(result, pattern=None, multiple:bool=False, flags=re.DOTALL):
    matches = [match.group(1) for match in re.finditer(patterns[pattern], result, flags=flags)]
    return matches if multiple else ('none' if (len(matches) == 0) else matches[0])


def download_nuccore(nuccore_id:str, fn_dir=None):
    fn_path = os.path.join(fn_dir, f'{nuccore_id}.fn')
    if not os.path.exists(fn_path):
        result = Entrez.efetch(db='nuccore', id=nuccore_id, rettype='fasta').read()
        with open(fn_path, 'w') as f:
            f.write(result)
    return fn_path


def download_protein_info(protein_id):
    # TODO: The non-redundant protein sequences need to be handled differently.  
    info = {'id':protein_id}
    
    result = Entrez.efetch(db='protein', id=protein_id, rettype='html').read().decode('utf-8')
    nuccore_id = parse_entrez(result, pattern='nuccore_id')
    if nuccore_id == 'none':
        return info

    info['coordinates'] = parse_entrez(result, pattern='coordinates', flags=re.MULTILINE|re.DOTALL)
    info['nuccore_id'] = nuccore_id

    result = Entrez.efetch(db='nuccore', id=nuccore_id, rettype='html').read().decode('utf-8')

    info['organism'] = parse_entrez(result, pattern='organism')
    info['taxonomy'] = parse_entrez(result, pattern='taxonomy')

    return info

def download_protein(protein_id:str):
    
    result = Entrez.efetch(db='protein', id=protein_id, rettype='html').read().decode('utf-8')
    return parse_entrez(result, 'seq')


def download_nr_protein_info(protein_id):
    assert protein_id.startswith('WP'), f'download_nr_protein_info: Expected a protein accession beginning with WP, but got {id_}'
    result = Entrez.efetch(db='ipg', id=protein_id, rettype='html').read().decode('utf-8')
    
    info = list()
    pattern = r'<CDS  accver="(.+)" start="(\d+)" stop="(\d+)" strand="(.+)" taxid="(.+)" org="([^"]+)".+</CDSList></Protein>'
    for match in re.finditer(pattern, result):
        info_ = {'id':protein_id}
        info_['nuccore_id'] = match.group(1)
        info_['start'] = match.group(2)
        info_['stop'] = match.group(3)
        info_['strand'] = match.group(4)
        info_['organism'] = match.group(6)
        info_['html'] = match.group(0)
        info.append(info_)

    return info


def download_assembly_id(nuccore_id):
    '''Get the assembly ID associated with the specified nuccore ID.'''
    get_alternate_nuccore_id = lambda nuccore_id : nuccore_id[:8] + '0' * 7 # For some reason, sometimes you need to manualy adjust the IDs?
    
    try:
        result = Entrez.elink(dbfrom='nuccore', db='assembly', id=nuccore_id, rettype='html').read().decode('utf-8')
        assembly_link_id = parse_entrez(result, pattern='id', multiple=True)[-1] # The second ID is the actual link. 
    except:
        try:
            result = Entrez.elink(dbfrom='nuccore', db='assembly', id=get_alternate_nuccore_id(nuccore_id), rettype='html').read().decode('utf-8')
            assembly_link_id = parse_entrez(result, pattern='id', multiple=True)[-1] # The second ID is the actual link. 
        except:
            print(f'download_assembly_id: Error raised when finding assembly for nuccore ID {nuccore_id}')
            return 'none'

    result = Entrez.esummary(db='assembly', id=assembly_link_id, rettype='html').read().decode('utf-8')
    assembly_id = parse_entrez(result, pattern='assembly_id')
    return assembly_id 