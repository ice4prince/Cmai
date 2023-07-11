#!/usr/bin/env python
# coding: utf-8

# In[7]:


import pandas as pd
import numpy as np
from Bio import SeqIO
import os


# In[ ]:


parser = argparse.ArgumentParser(description='Parameters for pair model.')

# Add a optional argument
parser.add_argument('--code', type=str, help='the CLAnO directory',default = '/project/DPDS/Wang_lab/shared/BCR_antigen/code/CLAnO')
parser.add_argument('--input',type = str, help = 'the input files in csv which should include Antigen_id,BCR_Vh,BCR_CDR3h',default = 'data/example/binary_example.csv')
parser.add_argument('--fasta',type = str, help = 'if no sequence included in the input, the seperate fasta file of antigens is required',default =None)
parser.add_argument('--continuous', action='store_true', help='swtich the mode from binary to continuous, default mode is binary.')
parser.add_argument('--pre_dir',type = str, help='the directory to save the preprocessed data.',default = 'data/intermediates')

args = parser.parse_args()

CODE_DIR = args.code
INPUT = args.input
FASTA = args.fasta
CONT = args.continuous
PRE_DIR = args.pre_dir


# In[103]:


# CODE_DIR = '/project/DPDS/Wang_lab/shared/BCR_antigen/code/CLAnO'
# INPUT = "data/example/binary_example.csv"
# FASTA = "data/example/binary_antigens.fasta"
# CONT = False
# PRE_DIR = 'data/intermediates'


# In[30]:


os.chdir(CODE_DIR)
if CONT:
    MODE = 'continuous'
    print('mode switching ...')
else:
    MODE = 'binary'
print('Entering',MODE,'mode now...')


# In[69]:


def write_fasta(df,filename):
    df = df.drop_duplicates(subset='Antigen_id')
    with open(filename, 'w') as f:
        for index, row in df.iterrows():
            # write each row to the file in FASTA format
            f.write('>' + str(row['Antigen_id']) + '\n')
            f.write(str(row['Antigen_seq']) + '\n')


# In[65]:


def parse_fasta(filename):
    with open(filename, 'r') as file:
        id_seq_dict = {}
        for line in file:
            line = line.rstrip()  # remove newline character
            if line.startswith('>'):
                id_ = line[1:]  # remove '>'
                id_seq_dict[id_] = ''
            else:
                id_seq_dict[id_] += line
        return id_seq_dict


# In[84]:


def add_BCR_id(df):
    df_copy = df.copy()
    df_copy['BCR'] = df_copy['BCR_Vh']+df_copy['BCR_CDR3h']
    unique_B=df_copy['BCR'].unique()
# assign names to each unique colB value
    dict_BCR = {unique_B[i]: f'BCR_{i}' for i in range(len(unique_B))}
    df['BCR_id'] = df_copy['BCR'].map(dict_BCR)
    return df


# In[98]:


def check_input(df,fastaname=None):
    global MODE
    print('Checking the columns...')    
    if MODE == 'continuous':
        required_cols = ['Antigen_id', 'BetterBCR_Vh', 'BetterBCR_CDR3h','WorseBCR_Vh','WorseBCR_CDR3h']
    else:
        required_cols = ['Antigen_id','BCR_Vh','BCR_CDR3h']
    for col in required_cols:
        if col not in df.columns:
            print(f"{col} is missing.")
            exit() 
    print('All required columns are included.')# This exits the function. If you want to exit the script completely, use `exit()`
    if 'Antigen_seq' in df.columns:
        write_fasta(df,PRE_DIR+'/antigens.fasta')
        print('Write antigen sequence to intermediates folder.')
    else:
        fasta_dict = parse_fasta(fastaname)
        df['Antigen_seq'] = df['Antigen_id'].map(fasta_dict)
        print('read in the antigen sequence from the fasta file.')
    if MODE == 'binary':
        if 'BCR_id' not in df.columns:
            df = add_BCR_id(df)
            print('BCR_id is added.')
    if MODE == 'continuous':
        if 'BetterBCR_id' not in df.columns:
            df_better = df[['BetterBCR_Vh','BetterBCR_CDR3h']]
            df_worse = df[['WorseBCR_Vh','WorseBCR_CDR3h']]
            df_better.columns= df_worse.columns = ['BCR_Vh','BCR_CDR3h']
            df_better = add_BCR_id(df_better)
            df_worse = add_BCR_id(df_worse)
            df['BetterBCR_id'] = df_better['BCR_id']
            df['WorseBCR_id'] = df_worse['BCR_id']
            print('BetterBCR_id and worseBCR_id are both added.')
    print(df.head())
    return df


# In[99]:


input_file = pd.read_csv(INPUT)
input_checked = check_input(input_file)


# In[107]:


input_checked.to_csv(PRE_DIR+'/processed_input.csv')
