#!/usr/bin/env python
# coding: utf-8

# In[ ]:





# In[2]:






import os
import sys
from sys import exit
import pandas as pd
import numpy as np
import random
from torch.utils.data import Dataset, DataLoader, TensorDataset, SubsetRandomSampler,Sampler
import torch.multiprocessing as mp
import torch
import torch.nn as nn
import torch.nn.functional as F
import math
import time
import colorama
from colorama import Fore,Back,Style
import matplotlib.pyplot as plt
import csv
import gc
import pickle
from Bio import SeqIO
module_path = os.path.abspath(os.path.join('..'))
if module_path not in sys.path:
    sys.path.append(module_path)
from datetime import date
from datetime import datetime
import torch.optim as optim
from torch.optim.lr_scheduler import ReduceLROnPlateau
import argparse
from glob import glob
from sklearn.model_selection import train_test_split
print("Today's date:",date.today())
print(str(datetime.now()))


# In[3]:






current = os.path.dirname(os.path.realpath(__file__))
parent = os.path.dirname(current)


parser = argparse.ArgumentParser(description='Parameters for the binding model.')

# Add a optional argument
parser.add_argument('--code', type=str, help='the Cmai directory',default = parent)
parser.add_argument('--input',type = str, help = 'the input folder for the preprocessed input',default = os.path.join(parent,'data/example/output'))
parser.add_argument('--npy_dir',type = str, help = 'the npy folder if different with input directory/NPY',default = None)
parser.add_argument('--background_npy',type = str, help = 'the npy folder for background antigens if different with the input directory/NPY',default = None)
parser.add_argument('--out',type = str, help = 'the directory for output files',default = os.path.join(parent,'data/example/output'))
# parser.add_argument('--species', action='store_true', help='match the species of background BCR to the target BCR. NOTE: the species MUST BE specified and unique in the target BCR input.')
parser.add_argument('--seed', type=int, help='the seed for the first 100 background BCRs. To use the prepared embeded 100 BCRs, keep the seed to default 1',default = 1)
parser.add_argument('--subsample1', type=int, help='the initial sample size of background BCRs. The default is 100',default = 100)
parser.add_argument('--subsample2', type=float, help='the initial sample size ratio of background Antigens. The default is 0.03',default = 0.03)
parser.add_argument('--bottomline1', type=int, help='the maximum size for subsample of background BCRs, which should no more than 1000000. The deafult is 10000',default = 10000)
parser.add_argument('--bottomline2', type=float, help='the maximum size ratio for subsample of background antigens, which should no more than 100%. The deafult is 0.3',default = 0.3)
parser.add_argument('--no_rank', action='store_true', help='Only export the predicted score but no rank in background BCRs, default is False.')
parser.add_argument('--export_background',action='store_true', help='Only export the score dict for background BCRs of amount of the bottomline number, default is False.')
parser.add_argument('--verbose', action='store_true', help='Enable verbose output, default is False.')
parser.add_argument('--no_merge', action='store_true', help='Unable merging output to input, default is False.')
parser.add_argument('--debug', action='store_true', help='Enable debug mode and print intermediates output every step.')
parser.add_argument('--Antigen_only', action='store_true', help='Only get the rank% of anchored antigen in background BCRs. Default is False')
parser.add_argument('--BCR_only', action='store_true', help='Only get the rank% of anchored BCRs in background antigens. Default is False')


args = parser.parse_args()

CODE_DIR = args.code
INPUT_DIR = args.input

OUT_DIR = args.out
BACK_BATCH_SIZE = 1
BATCH_SIZE = 1
# MATCHING_SPECIES = args.species
SEED = args.seed
SUBSAMPLE1 = args.subsample1
SUBSAMPLE2 = args.subsample2
BOTTOMLINE1 = args.bottomline1
BOTTOMLINE2 = args.bottomline2
ANTIGEN_ONLY = args.Antigen_only #anchor antigen
BCR_ONLY = args.BCR_only #anchor BCR
VERBOSE = args.verbose
DEBUG = args.debug
if args.npy_dir is not None:
    NPY_DIR = args.npy_dir
else:
    NPY_DIR = INPUT_DIR+'/NPY' ###need to add a command to move the pair.npy under results/pred/ to the intermediates/
# BACK_DIR = args.background_npy
if args.background_npy is not None:
    BACK_DIR = args.background_npy
else:
    BACK_DIR = NPY_DIR
if ANTIGEN_ONLY and not BCR_ONLY:
    print(Fore.GREEN+'Only Ranking the antigens...'+Style.RESET_ALL)
elif BCR_ONLY and not ANTIGEN_ONLY:
    print(Fore.GREEN+'Only Ranking the BCRs...'+Style.RESET_ALL)
else:
    print(Fore.GREEN+'Ranking both antigens and BCRs...'+Style.RESET_ALL)

if VERBOSE:
    print('Verbose Mode is: ON!')


# In[124]:







# current = '/project/DPDS/Wang_lab/shared/BCR_antigen/code/Cmai/scripts'
# parent = '/project/DPDS/Wang_lab/shared/BCR_antigen/code/Cmai'
# CODE_DIR = parent
# INPUT_DIR = os.path.join(parent,'data/example/output1_2')

# OUT_DIR = os.path.join(parent,'data/example/output')
# BACK_BATCH_SIZE = 1
# BATCH_SIZE = 1
# # MATCHING_SPECIES = args.species
# SEED = 1
# SUBSAMPLE1 = 100
# SUBSAMPLE2 = 0.03
# BOTTOMLINE1 = 10000
# BOTTOMLINE2 = 0.3
# BOSE = False
# DEBUG = False
# BCR_ONLY = False
# ANTIGEN_ONLY = False
# VERBOSE = True
# DEBUG = False
# # NPY_DIR = '/project/DPDS/Wang_lab/shared/BCR_antigen/data/rfscript/Antigen_embed/CLEANED/NPY'
# NPY_DIR = '/project/DPDS/Wang_lab/shared/BCR_antigen/data/output/Jun/NPY'
# BACK_DIR = '/project/DPDS/Wang_lab/shared/BCR_antigen/data/rfscript/Antigen_embed/CLEANED/NPY'


# In[5]:




# INPUT_DIR = '/project/DPDS/Wang_lab/shared/BCR_antigen/data/output/output_Catherin1'
# NPY_DIR = INPUT_DIR+'/NPY'
# OUT_DIR = '/project/DPDS/Wang_lab/shared/BCR_antigen/data/output/output_Catherin1'


# In[125]:






np.random.seed(SEED)
# torch.use_deterministic_algorithms(True, warn_only=True)
torch.manual_seed(SEED)


# CODE_DIR = '/project/DPDS/Wang_lab/shared/BCR_antigen/code/Cmai'
# #CUTOFF = 1800
# SEED = 1
# VERBOSE = False
# MATCHING_SPECIES = False
# SUBSAMPLE = 100
# BOTTOMLINE = 10000
# OUT_DIR = '/project/DPDS/Wang_lab/shared/BCR_antigen/code/Cmai/data/example'


# In[126]:






if DEBUG:
    print('Entering DEBUG mode.')

if BOTTOMLINE1 >1000000:
    print('The bottomline cannot be larger than 1,000,000.')
    exit()


# In[127]:




os.chdir(CODE_DIR)

# os.chdir(current+'/wrapV')
# glob('wrapV/*')
# from Vwrap import embedV
# os.chdir(current+'/wrapCDR3')
# from CDR3wrap import embedCDR3
from wrapV.Vwrap import embedV ##input needs to be list of strings
from wrapCDR3.CDR3wrap import embedCDR3 ##input needs to be list of strings


# In[128]:



# os.chdir(parent)



# In[129]:






BACKGROUND = 'data/background/backgroundBCR.csv.gz'
MODEL = 'models/antigenModel.pth'
MODEL2 = 'models/bcrModel.pth'
# NPY_DIR = INPUT_DIR+'/NPY'
INPUT = INPUT_DIR+'/processed_input.csv'
##MARK HERE

CHANNEL_ANTIGEN = 600
CLIP = None
LR = 0.005
T = 0.005
LAMBDA = 0
w = 100


# In[130]:






print('system version:',sys.version)
device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")
print('device:',device)
torch.set_printoptions(precision=10)


# In[131]:






# def get_now():
#     now = str(datetime.now())
#     tosec = now.split('.')[0]
#     out_now = tosec.split(' ')[0]+'_'+tosec.split(' ')[1].replace(':',"-")
#     return out_now

# def filter_big_antigens(dataset,cutoff):
#     dataset['aalens']=list(map(len,dataset['Antigen'].tolist()))
#     data_filtered = dataset.loc[dataset['aalens']< cutoff]
#     print('After removing antigens larger than '+str(cutoff)+', '+str(100*data_filtered.shape[0]/dataset.shape[0])+'% antigens remained.')
#     return data_filtered
#
# def check_bad_bcr(seq):
#     allowed_letters = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
#     uppercase_string = seq.upper()
#     return not any(char not in allowed_letters for char in uppercase_string)
#
# def filter_bad_bcr(df):
#     substrings = ['Vh', 'CDR3h', 'Antigen_seq']
#     some_columns = [col for col in df.columns if any(sub in col for sub in substrings)]
#     for col in some_columns:
#         df[col] = df[col].apply(lambda x: x.replace(' ', ''))
#     mask = df[some_columns].applymap(check_bad_bcr)
#     filtered_df = df[mask.all(axis=1)]
#     return filtered_df
#
# def preprocess(df):
#     df = filter_bad_bcr(df)
# #    df = filter_big_antigens(df,cutoff)
#     df = df.sort_values('Antigen_id')
#     df = df.assign(record_id = ['record_' + str(s) for s in range(df.shape[0])])
#     return df

# def has_space(string):
#     return " " in string
# def check_bad_letters(df,binary = True):
#     substrings = ['Vh', 'CDR3h', 'Antigen']
#     some_columns = [col for col in df_copy.columns if any(sub in col for sub in substrings)]
#     for col in some_columns:
#         if not all(df[col].apply(check_bad_bcr)):
#             print(Fore.RED +str(col)+' contains uncommon aa or symbols!')
#         else:
#             print(Fore.GREEN +str(col)+' PASS!')
#     print(Style.RESET_ALL)

def check_npy(df):
    for antigen in df['Antigen_id'].unique():
        if not os.path.exists(NPY_DIR+'/'+antigen+'.pair.npy'):
            print('The embedding for antigen: '+antigen+' is not found in NPY directory, Skipping...')
            deleted_rows = df[df['Antigen_id']==antigen]
            with open(INPUT_DIR+'/Skipped_entry.txt','a') as report:
                for _,row in deleted_rows.iterrows():
                    report.write(','.join(map(str, row.values)) + '\n')
            df = df[df['Antigen_id']!=antigen]
    return df

def build_BCR_dict(dataset,colname,precise = False):
    cols = dataset.filter(like = colname)
    uniq_keys = pd.unique(cols.values.ravel()).tolist()
    if colname == 'CDR3h':
        uniq_embedding,_,input_keys = embedCDR3(uniq_keys,precise = precise)
    elif colname == 'Vh':
        uniq_embedding,_,input_keys = embedV(uniq_keys,precise = precise)
    i = 0
    mydict = {}
    for key in uniq_keys:
        mydict[key] = uniq_embedding[i]
        i += 1
    return(mydict)

def predict_size(len,datatype = 'float32'):
    # Define the shape of the tensor
    shape = [1, len, len, 318]
    element_size = np.dtype(datatype).itemsize
    num_elements = math.prod(shape)
    tensor_size = num_elements * element_size
#    size_gb = tensor_size/(1024**3)
#    print(f"The tensor takes up {size_gb:.2f} GB of memory.")
#    size_gb = tensor_size/(1024*1024*1024)
#    print(f"Tensor size: {size_gb:.2f} GB")
    return tensor_size

def get_antigen_dict(df,npy_dir=NPY_DIR):
    antigen_pool = df['Antigen_id'].unique()
    antigen_dict = {}
    for antigen in antigen_pool:
        antigen_dict[antigen] = np.load(npy_dir+'/'+antigen+'.pair.npy').astype(np.float32)/w
    return antigen_dict

def get_bcr_dict(df):
    BCR_dict = {}
    for _, row in df.iterrows():
        key = row['BCR_id']
        value = (row['BCR_Vh'], row['BCR_CDR3h'])
        if key not in BCR_dict:
            BCR_dict[key] = set()
        BCR_dict[key].add(value)
    return BCR_dict


# In[132]:






class checkDataset(Dataset):
    def __init__(self, dataframe, antigen_dict, antigen_fpath_dict,len_dict):
        self.dataframe = dataframe
#        self.group_ids = self.generate_group_index(group_size=GROUP_SIZE)  # Adjust the group_size as needed
        self.your_data_list = self.get_my_data_list()
        self.antigen_fpath_dict = antigen_fpath_dict
#        self.bcr_pool = self.data[['BCR_id','BCR_Vh','BCR_CDR3h']].to_dict(orient='records')
        self.antigen_dict = antigen_dict
        self.cdr3_dict = {}
        self.v_dict = {}
        self.lens_dict = len_dict
        self.antigen_in = {}

    def __getitem__(self, idx):
        your_dict = self.your_data_list[idx]
#         self.antigen_feat = self.extract_antigen(antigen)[0].to(device)
#         self.lengthen = len(self.antigen_dict[antigen])
#         bcr_feat = self.__embedding_BCR(cdr3_key,v_key,precise = True)
#         pair_feat = self.__comb_embed_gpu(bcr_feat)
#         return pair_feat,index_key
        index_key = your_dict['record_id']
        if DEBUG:
            print('EntryID:'+str(index_key))
        antigen_key = your_dict['Antigen_id']
        if DEBUG:
            print('Antigen:'+str(antigen_key))
        aalens_key = self.lens_dict[antigen_key]
        bcr_key = your_dict['BCR_id']
        if DEBUG:
            print('BCR:'+str(bcr_key))
        cdr3_key = your_dict['BCR_CDR3h']
        if DEBUG:
            print('BCR_CDR3h:'+str(cdr3_key))
        v_key = your_dict['BCR_Vh']
        if DEBUG:
            print('BCR_Vh:'+str(v_key))
        bcr_feat = self.__embedding_BCR(cdr3_key,v_key,precise = True) #MARK HERE: EDIT the embedding_BCR function using a dictionary of bcr_id:bcr_embedding
        if antigen_key not in self.antigen_in.keys():
            self.__get_antigen_in(antigen_key)
#            print(self.antigen_in)
#        print('antigen shape',self.antigen_in[antigen_key].shape)
        pair_feat = self.__comb_embed_gpu(self.antigen_in[antigen_key],bcr_feat)
        return pair_feat, index_key, antigen_key, bcr_key#better_feat, worse_feat,

    def get_my_data_list(self,selected_cols = 'BCR|BCR|Antigen|record_id|'):
        ds_to_dict = self.dataframe.filter(regex=selected_cols)
    #ds_to_dict = dataset[selected_cols].set_index('record_id')
        my_data_list = ds_to_dict.to_dict(orient='records')
        return my_data_list

    def __comb_embed_gpu(self,antigen_feat,BCR_feat):
#        print('antigen to comb:',antigen_name)
        lengthen = CHANNEL_ANTIGEN
        #rint('The current antigen is: ',antigen_name)
#        print('length get from dict_len:',lengthen)
        single_antigen_g = antigen_feat[0]
#        single_antigen_g = F.normalize(single_antigen_g, p=2, dim=3)
#        print('shape of antigen from antigen dict:',single_antigen_g.shape)
#        single_antigen_g = torch.from_numpy(single_antigen).to(device)
        single_BCR_g = torch.from_numpy(BCR_feat).to(device)
#        print('single BCR shape:',single_BCR_g.shape)
#        single_BCR_g = torch.from_numpy(BCR_feat).half().to(device)
        BCR_t = torch.tile(single_BCR_g,(lengthen,lengthen,1))
#        print('tiled bcr shape',BCR_t.shape)
#        print('shape of BCR_tiled:',BCR_t.shape)#empty
        pair_feat_g = torch.cat((single_antigen_g,BCR_t),dim=2)
        del single_BCR_g,BCR_t
        torch.cuda.empty_cache()
        return pair_feat_g#.half()


    def __get_antigen_in(self,antigen_name):
        #print('Next antigen:',antigen_name)
#        if not antigen_name in self.antigen_in:
        self.antigen_in.clear()
        if not antigen_name in self.antigen_dict:
            antigen_to_in = self.extract_antigen(antigen_name)
        else:
            antigen_to_in = self.antigen_dict[antigen_name]
        try:
            antigen_tensor = torch.from_numpy(antigen_to_in)
            self.antigen_in[antigen_name] = self.pool_antigen(antigen_tensor,CHANNEL_ANTIGEN).to(device) ###ON CPU

        except RuntimeError as e:
            if "CUDA out of memory" in str(e):
                print("CUDA out of memory. Clearing cache and trying again...")
                torch.cuda.empty_cache()
                try:
                    self.antigen_in[antigen_name] = self.pool_antigen(torch.from_numpy(antigen_to_in),CHANNEL_ANTIGEN).to(device)
                except RuntimeError as e:
                    if "CUDA out of memory" in str(e):
                        print("Still out of memory after clearing cache.")
                    else:
                        raise
            else:
                raise

    def pool_antigen(self,antigen_input,out_n_channel):
#        lengthen = antigen_input.shape[1]
        pooling_layer = nn.AdaptiveAvgPool2d((out_n_channel,out_n_channel))
        output = pooling_layer(antigen_input.permute(0,3,1,2)).permute(0,2,3,1)
        return output

    def extract_antigen(self,antigen_name,verbose=False):
        try:
            antigen_import = np.load(self.antigen_fpath_dict+'/'+antigen_name+'.pair.npy').astype(np.float32)/w
            if not antigen_import.shape[1] == self.lens_dict[antigen_name]:
                print(Fore.RED + 'antigen ' +str(antigen_name)+' embedding '+str(antigen_import.shape[1])+' is NOT in the correct shape '+str(self.lens_dict[antigen_name])+'!'+ Style.RESET_ALL)
                exit()
#             single_antigen = antigen_import
            self.antigen_dict[antigen_name] = antigen_import

        except ValueError:
            print('The embedding of antigen %s cannot be found!' % antigen_name)

        return single_antigen

    def __embedding_BCR(self,cdr3_seq,v_seq,precise = False):
        if cdr3_seq not in self.cdr3_dict:
#            print('CDR3 not in dictionary!!')
#            df1 = pd.DataFrame()
            cdr3_feat,*_ = embedCDR3([cdr3_seq],precise = precise)
            cdr3_feat = cdr3_feat[0]
            self.cdr3_dict[cdr3_seq]=cdr3_feat
        else:
#            print('CDR3 in dictionary!!')
            cdr3_feat = self.cdr3_dict[cdr3_seq]
        if DEBUG:
            print('cdr3:'+str(cdr3_feat))
        if v_seq not in self.v_dict:
#            print('V not in dictionary!!')
#            df2 = pd.DataFrame([v_seq])
            v_feat,*_ = embedV([v_seq],precise = precise)
            v_feat = v_feat[0]
            self.v_dict[v_seq]=v_feat
        else:
#            print('V in dictionary!!')
            v_feat = self.v_dict[v_seq]
        if DEBUG:
            print('v:'+str(v_feat))
        bcr_feat = np.concatenate((cdr3_feat,v_feat))
        return bcr_feat

    def __len__(self):
        return len(self.your_data_list)


# In[133]:





# to geth rank_antigen among background BCRs
class rankDataset(Dataset):
    def __init__(self, dataframe, antigen, cdr3_dict, v_dict,subsample_ratio=1/10000,seed=SEED):
        self.seed = seed
#         self.antigen_fpath_dict = antigen_fpath_dict
        self.antigen_feat = antigen
        self.cdr3_dict = cdr3_dict
        self.v_dict = v_dict
        self.background = self.subsample_data(dataframe, subsample_ratio=subsample_ratio)
        self.bcr_pool = self.background[['Vh','CDR3h']].to_dict(orient='records')
#        print(len(self.bcr_pool))
#         self.lengthen = len(self.antigen_dict[antigen])
#        print(self.lengthen)

    def subsample_data(self, dataframe, subsample_ratio):
        if not (0 < subsample_ratio <= 1):
            raise ValueError("subsample_ratio must be between 0 and 1")
        n = int(subsample_ratio * len(dataframe))  # Compute number of rows to select
        return dataframe.iloc[:n]
    def __getitem__(self, idx):
        bcr_dict = self.bcr_pool[idx]
#         print(idx)
#        index_key = bcr_dict['ID']
        v_key = bcr_dict['Vh']
        cdr3_key = bcr_dict['CDR3h']
        bcr_feat = self.__embedding_BCR(cdr3_key,v_key,precise = True)
        pair_feat = self.__comb_embed_gpu(bcr_feat)
        return pair_feat

    def __comb_embed_gpu(self,bcr_feat):
#         if MODE == 'binary':
        lengthen = CHANNEL_ANTIGEN
#         else:
#             lengthen = self.lengthen
        antigen_tensor =  torch.from_numpy(self.antigen_feat).float()
        single_antigen_g = self.pool_antigen(antigen_tensor,CHANNEL_ANTIGEN)[0].to(device)
        single_BCR_g = torch.from_numpy(bcr_feat).to(device)
        BCR_t = torch.tile(single_BCR_g,(lengthen,lengthen,1))
        pair_feat_g = torch.cat((single_antigen_g,BCR_t),dim=2)
        del single_BCR_g,BCR_t
        torch.cuda.empty_cache()
        return pair_feat_g

#     def extract_antigen(self,antigen_name):
#         try:
#             antigen_import = np.load(str(self.antigen_fpath_dict+'/'+antigen_name+'.pair.npy'))/w
# #             if not antigen_import.shape[1] == self.lengthen:
# #                 print(Fore.RED + 'antigen ' +str(antigen_name)+' embedding'+str(antigen_import.shape[1])+' is NOT in the correct shape '+str(self.lens_dict[antigen_name])+'!'+ Style.RESET_ALL)
# #                 exit()
#             single_antigen = torch.from_numpy(antigen_import).float()
# #             if MODE == 'binary':
#             single_antigen = self.pool_antigen(single_antigen,CHANNEL_ANTIGEN) ###ON CPU
#         except ValueError:
#             print('The embedding of antigen %s cannot be found!' % antigen_name)
#         return single_antigen

    def __embedding_BCR(self,cdr3_seq,v_seq,precise = True):
        if cdr3_seq not in self.cdr3_dict:
            cdr3_feat,*_ = embedCDR3([cdr3_seq],precise = precise)
            cdr3_feat = cdr3_feat[0]
            self.cdr3_dict[cdr3_seq]=cdr3_feat
        else:
            cdr3_feat = self.cdr3_dict[cdr3_seq]
#         print(cdr3_feat)
        if v_seq not in self.v_dict:
            v_feat,*_ = embedV([v_seq],precise = precise)
            v_feat = v_feat[0]
            self.v_dict[v_seq]=v_feat
        else:
            v_feat = self.v_dict[v_seq]
#         print(v_feat)
        bcr_feat = np.concatenate((cdr3_feat,v_feat))
        return bcr_feat

    def pool_antigen(self,antigen_tensor,out_n_channel):
#        lengthen = antigen_input.shape[1]
        pooling_layer = nn.AdaptiveAvgPool2d((out_n_channel,out_n_channel))
        output = pooling_layer(antigen_tensor.permute(0,3,1,2)).permute(0,2,3,1)
        return output

    def __len__(self):
        return len(self.bcr_pool)


# In[134]:





###to get rank_bcr with background antigens
class rankDataset2(Dataset):
    def __init__(self,bcr_v,bcr_cdr3,antigen_ids,antigen_fpath_dict,subsample_ratio=0.1,seed = SEED):
        gc.collect()
        torch.cuda.empty_cache()
        self.seed = seed
        random.seed(self.seed)
        if not (0 < subsample_ratio <= 1):
            if (subsample_ratio >1):
                subsample_ratio = 1
                print('Resetting subsample_ratio to 1')
            else:
                raise ValueError("subsample_ratio must be between 0 and 1")
        n = int(subsample_ratio * len(antigen_ids))  # Compute number of rows to select
        self.antigen_ids = antigen_ids[:n]
        self.bcr_v=bcr_v
        self.bcr_cdr3=bcr_cdr3
        self.cdr3_dict = {}
        self.v_dict = {}
        self.antigen_fpath_dict = antigen_fpath_dict
        self.antigen_dict={}
        self.antigen_in = {}
        # print('The directory for background antigens:',self.antigen_fpath_dict)
        # self.lens_dict = {}
    def __getitem__(self,idx):
        bcr_feat = self.__embedding_BCR(self.bcr_cdr3,self.bcr_v,precise = True)
        antigen_key = self.antigen_ids[idx]
        # print('The randomly picked antigen is:'+str(antigen_key))
        self.__get_antigen_in(antigen_key)
        # antigen_feat = self.antigen_in[antigen_key]
#        print('antigen shape',self.antigen_in[antigen_key].shape)
        pair_feat = self.__comb_embed_gpu(antigen_key,bcr_feat)
        return pair_feat

    def __comb_embed_gpu(self,antigen_name,BCR_feat):
#        print('antigen to comb:',antigen_name)
        lengthen = CHANNEL_ANTIGEN
        #rint('The current antigen is: ',antigen_name)
#        print('length get from dict_len:',lengthen)
        single_antigen_g = self.antigen_in[antigen_name][0]
        self.antigen_in[antigen_name].cpu()
        del self.antigen_in[antigen_name]
#        single_antigen_g = F.normalize(single_antigen_g, p=2, dim=3)
#        print('shape of antigen from antigen dict:',single_antigen_g.shape)
#        single_antigen_g = torch.from_numpy(single_antigen).to(device)
        single_BCR_g = torch.from_numpy(BCR_feat).to(device)
#        print('single BCR shape:',single_BCR_g.shape)
#        single_BCR_g = torch.from_numpy(BCR_feat).half().to(device)
        BCR_t = torch.tile(single_BCR_g,(lengthen,lengthen,1))
#        print('tiled bcr shape',BCR_t.shape)
#        print('shape of BCR_tiled:',BCR_t.shape)#empty
        pair_feat_g = torch.cat((single_antigen_g,BCR_t),dim=2)
        del single_BCR_g,BCR_t
        torch.cuda.empty_cache()
        return pair_feat_g#.half()

    # def __import_antigen(self,antigen_name):
    #     if not antigen_name in self.antigen_dict:
    #         antigen_to_in = self.extract_antigen(antigen_name)
    #     else:
    #         antigen_to_in = self.antigen_dict[antigen_name]
    #     return antigen_to_in

    def __toGPU_and_pool(self,antigen_name):
        # print('pooling and importing the antigen',antigen_name,'to GPU...')
        antigen_tensor = torch.from_numpy(self.extract_antigen(antigen_name))
        # print('the size of antigen:',antigen_name,'is:',antigen_tensor.shape)
        self.antigen_in[antigen_name] = self.pool_antigen(antigen_tensor,CHANNEL_ANTIGEN).to(device)

    def __get_antigen_in(self,antigen_name):
        if antigen_name not in self.antigen_in:
            try:
                self.__toGPU_and_pool(antigen_name)
            except RuntimeError as e:
                if "CUDA out of memory" in str(e):
                    print("Out of memory when importing the anchor antigen. Clear gpu and try again.")
                    for key in list(self.antigen_in.keys()):
                        self.antigen_in[key].cpu()
                        del self.antigen_in[key]
                    self.antigen_in.clear()
                    gc.collect()
                    torch.cuda.empty_cache()
                    try:
                        self.__toGPU_and_pool(antigen_name)
                    except RuntimeError as e:
                        if "CUDA out of memory" in str(e):
                            print("Still out of memory after clearing GPU when importing the anchor antigen.")
                        else:
                            print("Runtime error during importing the anchor antigen after clearing GPU:", str(e))
                            raise
                else:
                    print("Runtime error during importing the anchor antigen before clearing GPU:", str(e))
                    raise


    def pool_antigen(self,antigen_input,out_n_channel):
#        lengthen = antigen_input.shape[1]
        pooling_layer = nn.AdaptiveAvgPool2d((out_n_channel,out_n_channel))
        output = pooling_layer(antigen_input.permute(0,3,1,2)).permute(0,2,3,1)
        return output


    def extract_antigen(self,antigen_name,verbose=False):
        if antigen_name in self.antigen_dict:
            single_antigen = self.antigen_dict[antigen_name]
        else:
            try:
                # print('The directory for background antigens:',self.antigen_fpath_dict)
                # print('The antigen name:',antigen_name)
                antigen_import = np.load(str(self.antigen_fpath_dict+'/'+antigen_name+'.pair.npy')).astype(np.float32)/w
                # if not antigen_import.shape[1] == self.lens_dict[antigen_name]:
                #     print(Fore.RED + 'antigen ' +str(antigen_name)+' embedding '+str(antigen_import.shape[1])+' is NOT in the correct shape '+str(self.lens_dict[antigen_name])+'!'+ Style.RESET_ALL)
                #     exit()
                single_antigen = antigen_import
#                single_antigen = torch.from_numpy(antigen_import).to(device) ###ON GPU
#            print(npy.shape)
                self.antigen_dict[antigen_name] = single_antigen
                if verbose:
                    print(Fore.RED + 'New antigen added to dictionary:'+antigen_name+Fore.RESET)
            except ValueError:
                print('The embedding of antigen %s cannot be found!' % antigen_name)
#             if verbose:
#                 print(Fore.RED + 'New antigen added to dictionary:',antigen_name)
#                 print(Fore.RED + 'Number of antigens included:',len(self.antigen_dict))
#                 print(Style.RESET_ALL)
        return single_antigen

    def __embedding_BCR(self,cdr3_seq,v_seq,precise = False):
        if cdr3_seq not in self.cdr3_dict:
#            print('CDR3 not in dictionary!!')
#            df1 = pd.DataFrame()
            cdr3_feat,*_ = embedCDR3([cdr3_seq],precise = precise)
            cdr3_feat = cdr3_feat[0]
            self.cdr3_dict[cdr3_seq]=cdr3_feat
        else:
#            print('CDR3 in dictionary!!')
            cdr3_feat = self.cdr3_dict[cdr3_seq]
        if v_seq not in self.v_dict:
#            print('V not in dictionary!!')
#            df2 = pd.DataFrame([v_seq])
            v_feat,*_ = embedV([v_seq],precise = precise)
            v_feat = v_feat[0]
            self.v_dict[v_seq]=v_feat
        else:
#            print('V in dictionary!!')
            v_feat = self.v_dict[v_seq]
        bcr_feat = np.concatenate((cdr3_feat,v_feat))
        return bcr_feat

    def __len__(self):
        return len(self.antigen_ids)


# In[135]:






class SelfAttentionPooling(nn.Module):
    """
    Implementation of SelfAttentionPooling
    Original Paper: Self-Attention Encoding and Pooling for Speaker Recognition
    https://arxiv.org/pdf/2008.01077v1.pdf
    """
    def __init__(self, input_dim,hidden_dim):
        super(SelfAttentionPooling, self).__init__()
#        hidden_dim=10
        self.W1 = nn.Linear(input_dim, hidden_dim)
        self.relu = nn.LeakyReLU(0.1)
        self.W2 = nn.Linear(hidden_dim, 1)

    def forward(self, batch_rep):
        """
        input:
            batch_rep : size (N, T, H), N: batch size, T: sequence length, H: Hidden dimension

        attention_weight:
            att_w : size (N, T, 1)

        return:
            utter_rep: size (N, H)
        """
        att_w = F.softmax(self.W2(self.relu(self.W1(batch_rep))).squeeze(-1),dim=1).unsqueeze(-1)

        utter_rep = torch.sum(batch_rep* att_w, dim=1)

        return utter_rep


# In[136]:






class mix_model(nn.Module):
    def __init__(self):
        super(mix_model,self).__init__()
#         if mode =='binary':
        self.model1 = nn.Sequential(
            nn.Linear(318,40),#.to(torch.float64),
            # in (1,len,len,318)
            # out (1,len,len.50)
            nn.LeakyReLU(0.1),
            nn.Linear(40,30),#.to(torch.float64),
            nn.LeakyReLU(0.1),
            nn.Linear(30,20),#.to(torch.float64),
            # out (1,len,len,20)
            nn.LeakyReLU(0.1)
        )
#         else:
#             self.model10 = nn.Sequential(
#                 nn.Linear(318,40),#.to(torch.float64),
#                 # in (1,len,len,318)
#                 # out (1,len,len.50)
#                 nn.LeakyReLU(0.1),
#             )

#             self.model11 = nn.Sequential(
#                 nn.AdaptiveAvgPool2d((CHANNEL_ANTIGEN,CHANNEL_ANTIGEN))
#             )

#             self.model12 = nn.Sequential(
#                 nn.Linear(40,30),#.to(torch.float64),
#                 nn.LeakyReLU(0.1),
#                 nn.Linear(30,20),#.to(torch.float64),
#                 # out (1,len,len,20)
#                 nn.LeakyReLU(0.1)
#             )
        self.model2 = SelfAttentionPooling(input_dim=20,hidden_dim=30)
        self.model2_1 = SelfAttentionPooling(input_dim=20,hidden_dim=30)
        # input_dim = hidden size (number of channels)
#        self.model2 = nn.MultiheadAttention(embed_dim=20,num_heads=N_HEADS,batch_first=True)
        self.model3 = nn.Sequential(
            nn.Linear(20,15),
            nn.LeakyReLU(0.1),
            nn.Linear(15,1)
        )
#         self.alpha1 = nn.Parameter(torch.randn(1))
#         self.beta1 = nn.Parameter(torch.randn(1))
#         self.alpha2 = nn.Parameter(torch.randn(1))
#         self.beta2 = nn.Parameter(torch.randn(1))
    def forward(self,x): ###because in getitem, return is .cuda(), now input is on gpu
#         x = torch.empty(0)
#         x = x.to(device)
#        x = x.permute(0,2,1,3)
#         print('after permute',x.shape)
        x0 = torch.empty(0)
        x0 = x0.to(device)
#         if mode=='binary':
        x = self.model1(x)
        for i in range(len(x)):
            k = x[i]
            k = self.model2(k).unsqueeze(0)
            k = self.model2_1(k)
            x0 = torch.cat((x0, k), dim=0)
#         else:
#             x = self.model10(x)
#             for i in range(len(x)):
#                 k = x[i]
#                 k=k.permute(2,0,1)
#                 k=self.model11(k)
#                 k=k.permute(1,2,0)
#                 k=self.model12(k)
#                 k = self.model2(k).unsqueeze(0)
#                 k = self.model2_1(k)
#                 x0 = torch.cat((x0, k), dim=0)
        x0 = F.normalize(x0)
        x0 = self.model3(x0).squeeze()
#        if binary:
        out  = x0
        return(out)
#
# class mix_model2(nn.Module):
#     def __init__(self):
#         super(mix_model2,self).__init__()
#         self.model1 = nn.Sequential(
#             nn.Linear(318,40),#.to(torch.float64),
#             # in (1,len,len,318)
#             # out (1,len,len.50)
#             nn.LeakyReLU(0.1),
#             nn.Linear(40,30),#.to(torch.float64),
#             nn.LeakyReLU(0.1),
#             nn.Linear(30,20),#.to(torch.float64),
#             # out (1,len,len,20)
#             nn.LeakyReLU(0.1)
#         )
#         self.model2 = SelfAttentionPooling(input_dim=20,hidden_dim=30)
#         self.model2_1 = SelfAttentionPooling(input_dim=20,hidden_dim=30)
#         # input_dim = hidden size (number of channels)
# #        self.model2 = nn.MultiheadAttention(embed_dim=20,num_heads=N_HEADS,batch_first=True)
#         self.model3 = nn.Sequential(
#             nn.Linear(20,15),
#             nn.LeakyReLU(0.1),
#             nn.Linear(15,1)
#         )
# #     def __init__(self):
# #         super(mix_model,self).__init__()
# #         self.model1 = nn.Sequential(
# #             nn.Linear(318,40),#.to(torch.float64),
# #             # in (1,len,len,318)
# #             # out (1,len,len.50)
# #             nn.LeakyReLU(0.1),
# #             nn.Linear(40,30),#.to(torch.float64),
# #             nn.LeakyReLU(0.1),
# #             nn.Linear(30,20),#.to(torch.float64),
# #             # out (1,len,len,20)
# #             nn.LeakyReLU(0.1)
# #         )
# #         self.model2 = SelfAttentionPooling(input_dim=20,hidden_dim=30)
# #         self.model2_1 = SelfAttentionPooling(input_dim=20,hidden_dim=30)
# #         # input_dim = hidden size (number of channels)
# # #        self.model2 = nn.MultiheadAttention(embed_dim=20,num_heads=N_HEADS,batch_first=True)
# #         self.model3 = nn.Sequential(
# #             nn.Linear(20,15),
# #             nn.LeakyReLU(0.1),
# #             nn.Linear(15,1)
# #         )
#         self.alpha1 = nn.Parameter(torch.randn(1))
#         self.beta1 = nn.Parameter(torch.randn(1))
#         self.alpha2 = nn.Parameter(torch.randn(1))
#         self.beta2 = nn.Parameter(torch.randn(1))
# #        self.sigmoid = nn.Sigmoid()
#     def forward(self,x):#,binary = True, is_10X = True): ###because in getitem, return is .cuda(), now input is on gpu
# #         x = torch.empty(0)
# #         x = x.to(device)
# #        x = x.permute(0,2,1,3)
# #         print('after permute',x.shape)
#         x = self.model1(x)
# #         print('after model1',x.shape)
#         x0 = torch.empty(0)
#         x0 = x0.to(device)
#         for i in range(len(x)):
#             k = x[i]
#             k = self.model2(k).unsqueeze(0)
# #             print('after model2',k.shape)
#             k = self.model2_1(k)
# #             print('after model2_1',k.shape)
#             x0 = torch.cat((x0, k), dim=0)
# #         print('after loop:',x0.shape)
#         x0 = F.normalize(x0)
#         x0 = self.model3(x0).squeeze()
#         # if binary:
#         out  = x0
#         # else:
#             # if is_10X:
#                 # out = x0*(-math.exp(self.alpha1))+self.beta1
#             # else:
#                 # out = x0*(-math.exp(self.alpha2))+self.beta2
#         return(out)


# In[137]:






def background_scores(dataloader,model):
    score_ls = []
    model.eval()
    for batch in dataloader:
        score = model(batch.unsqueeze(0))
        score_ls.append(score.item())
    return score_ls


# In[138]:






def check_score(dataloader,model1,model2):
    res = pd.DataFrame(columns=['record_id', 'Antigen', 'BCR_id','Score_antigen','Score_bcr'])
#    print(zip(index_idx,antigen_index))
#    print(b_pair.shape,w_pair.shape)
    model1.eval()
    model2.eval()
    for batch in dataloader:
        pair, index_idx, antigen_key, bcr_key =batch ##Change the InLoader, when binary is False, w_pair = score
        score1 = model1(pair.unsqueeze(0)).item()
        score2 = model2(pair.unsqueeze(0)).item()
        res_dict = {'record_id':index_idx,'Antigen':antigen_key,'BCR_id':bcr_key,'Score_antigen':score1,'Score_bcr':score2}
        df = pd.DataFrame({key: [value] for key, value in res_dict.items()})
        res = pd.concat([res,df],axis=0,ignore_index=True)
        # pd.concat([res,df],axis=0,ignore_index=True)
    return(res)


# In[139]:






def locate_rank(number, my_list):
    # Step 1: Sort the list in ascending order
    sorted_list = sorted(my_list)

    # Step 2: Insert the number into the sorted list while maintaining the sorted order
    from bisect import insort_left
    insort_left(sorted_list, number)

    # Step 3: Find the index of the inserted number in the list
    index = sorted_list.index(number)

    # Step 4: Calculate the percentage rank using the index and the length of the list
    percentage_rank = (index / len(sorted_list))

    return percentage_rank


# In[ ]:






# In[140]:


# for bcr_id, bcr_set in bcr_dict.items():
#     for bcr in bcr_set:
#         bcr_v,bcr_cdr3 = bcr
#         bcr_loader = rankDataset2(bcr_v,bcr_cdr3,antigen_ids,BACK_DIR,subsample_ratio=subsample2,seed = SEED)
#         bcr_score_background = background_scores(bcr_loader,model_mix2)
#         bcr_score_dict[bcr_id]=bcr_score_background


# In[141]:






def generate_score_dict(df,antigen_ids,antigen_dict,bcr_dict,cdr3_dict,v_dict,antigen_fpath_dict,model1,model2,
                        antigen_score_dict={},bcr_score_dict={},subsample1=1/10000,subsample2=0.1,seed =SEED,antigen_only = False,bcr_only=False):
    # print('The npy directory for backgroudn antigens:',antigen_fpath_dict)
    if antigen_only and not bcr_only:###antigen anchor, rank BCRs in background BCRs
        i = 0
        for antigen_id, antigen in antigen_dict.items():
            if VERBOSE:
                print(f"the {i}/{len(antigen_dict)} Antigen: {antigen_id}")
            background_loader = rankDataset(df, antigen, cdr3_dict, v_dict,subsample_ratio=subsample1,seed=SEED)
            antigen_score_background = background_scores(background_loader,model1)
            antigen_score_dict[antigen_id]=antigen_score_background
            i += 1
    elif bcr_only and not antigen_only: ###BCR_anchor, rank antigens in background antigens
        j= 0
        for bcr_id, bcr_set in bcr_dict.items():
            if VERBOSE:
                print(f"the {j}/{len(bcr_dict)} BCR: {bcr_id}")
            for bcr in bcr_set:
                bcr_v,bcr_cdr3 = bcr
                bcr_loader = rankDataset2(bcr_v,bcr_cdr3,antigen_ids,antigen_fpath_dict,subsample_ratio=subsample2,seed = SEED)
                bcr_score_background = background_scores(bcr_loader,model2)
                bcr_score_dict[bcr_id]=bcr_score_background
            j += 1

    else:
        i = 0
        j = 0
        for antigen_id, antigen in antigen_dict.items():
            if VERBOSE:
                print(f"the {i}/{len(antigen_dict)} Antigen: {antigen_id}")
            background_loader = rankDataset(df, antigen, cdr3_dict, v_dict,subsample_ratio=subsample1,seed=seed)
            antigen_score_background = background_scores(background_loader,model1)
            antigen_score_dict[antigen_id]=antigen_score_background
            i += 1
        for bcr_id, bcr_set in bcr_dict.items():
            if VERBOSE:
                print(f"the {j}/{len(bcr_dict)} BCR: {bcr_id}")
            for bcr in bcr_set:
                bcr_v,bcr_cdr3 = bcr
                bcr_loader = rankDataset2(bcr_v,bcr_cdr3,antigen_ids,antigen_fpath_dict,subsample_ratio=subsample2,seed = SEED)
                bcr_score_background = background_scores(bcr_loader,model2)
                bcr_score_dict[bcr_id]=bcr_score_background
            j += 1
    return(antigen_score_dict,bcr_score_dict)


# In[142]:


# antigen_score_dict,bcr_score_dict = generate_score_dict(background,antigen_ids,antigen_dict,bcr_dict,CDR3h_dict,Vh_dict,NPY_DIR,model_mix,model_mix2,
#                                                             antigen_score_dict=antigen_score_dict,bcr_score_dict=bcr_score_dict,
#                                                             subsample1=subsample1/1000000,subsample2=subsample2,seed = SEED,antigen_only = False,bcr_only=True)


# In[143]:


def get_subset_dict(df,antigen_dict,bcr_dict):
    sub_antigen_dict = {key: antigen_dict[key] for key in list(df['Antigen'].unique()) if key in antigen_dict}
    sub_bcr_dict  = {key: bcr_dict[key] for key in list(df['BCR_id'].unique()) if key in bcr_dict}
    return(sub_antigen_dict,sub_bcr_dict)


# In[144]:


# df = df_both
# antigen_score_dict = {}
# bcr_score_dict = {}
# sub_antigen_dict,sub_bcr_dict = get_subset_dict(df,antigen_dict,bcr_dict)
# print(sub_antigen_dict.keys(),sub_bcr_dict.keys())
# antigen_score_dict,bcr_score_dict = generate_score_dict(background,antigen_ids,sub_antigen_dict,sub_bcr_dict,CDR3h_dict,Vh_dict,NPY_DIR,model_mix,model_mix2,
#                                                             antigen_score_dict=antigen_score_dict,bcr_score_dict=bcr_score_dict,
#                                                             subsample1=subsample1/1000000,subsample2=subsample2,seed = SEED,antigen_only = False,bcr_only=False)
# df_updated = calculate_rank(df, antigen_score_dict, bcr_score_dict, antigen_only=False, bcr_only=False)
# # df_res = rank_check(df_updated,0.1,0.1)


# In[145]:


# print(df_updated)


# In[146]:


# res =rank_check(df_updated,0.1,0.1)
# print(res)


# In[147]:


# print(df_updated)
# df_updated=rank_check(df_updated,0.1,0.1)
# print(df_updated)


# In[148]:


# # calculate_rank(target,antigen_score_dict,bcr_score_dict,antigen_dict,len_dict,model_mix,model_mix2)
# check_loader = checkDataset(target, antigen_dict, NPY_DIR,len_dict)
# res_check = check_score(check_loader,model_mix,model_mix2)


# In[221]:


def get_check_score(df, antigen_dict,len_dict,model1,model2):
    check_loader = checkDataset(df, antigen_dict, NPY_DIR,len_dict)
    res_check = check_score(check_loader,model1,model2)
    return res_check

def determine_check(row):
    # Check which columns exist in the row
    has_check_antigen = 'check_antigen' in row.index
    has_check_bcr = 'check_bcr' in row.index

    # Case 1: Both 'check_antigen' and 'check_bcr' exist
    if has_check_antigen and has_check_bcr:
        if row['check_antigen'] == 'no' and row['check_bcr'] == 'no':
            return 'both'
        elif row['check_antigen'] == 'no':
            return 'antigen'
        elif row['check_bcr'] == 'no':
            return 'bcr'
        else:
            return 'pass'

    # Case 2: Only 'check_antigen' exists
    elif has_check_antigen:
        return 'antigen' if row['check_antigen'] == 'no' else 'pass'

    # Case 3: Only 'check_bcr' exists
    elif has_check_bcr:
        return 'bcr' if row['check_bcr'] == 'no' else 'pass'

    # Case 4: Neither column exists (shouldn't happen in practice)
    else:
        return 'pass'  # Default to 'pass' if both columns are missing

def rank_check(df, cutoff_antigen,cutoff_bcr,antigen_only=False,bcr_only=False):
    df = df.copy()  # Avoid modifying the original DataFrame directly

    has_rank_antigen = 'Rank_antigen' in df.columns
    has_rank_bcr = 'Rank_bcr' in df.columns

    # Remove 'check_antigen' and 'check_bcr' if they already exist
    df = df.drop(columns=['check_antigen', 'check_bcr', 'check'], errors='ignore')

    # Initialize check columns as 'yes' by default
    if not bcr_only:
        df['check_antigen'] = 'yes'
    if not antigen_only:
        df['check_bcr'] = 'yes'

    # Case 1: Neither Rank_antigen nor Rank_bcr exists
    if not has_rank_antigen and not has_rank_bcr:
        if not bcr_only:
            df['check_antigen'] = 'no'
        if not antigen_only:
            df['check_bcr'] = 'no'

    # Case 2: Rank_antigen missing ? Check Rank_bcr based on Score_bcr
    elif not has_rank_antigen:
        if not bcr_only:
            df['check_antigen'] = 'no'
        df['check_bcr'] = df['Rank_bcr'].apply(lambda x: 'no' if x <= cutoff_bcr else 'yes')

    # Case 3: Rank_bcr missing ? Check Rank_antigen based on Score_antigen
    elif not has_rank_bcr:
        if not antigen_only:
            df['check_bcr'] = 'no'
        df['check_antigen'] = df['Rank_antigen'].apply(lambda x: 'no' if x <= cutoff_antigen else 'yes')

    # Case 4: Both Rank_antigen and Rank_bcr exist ? Apply cutoff logic
    else:
        df['check_antigen'] = df['Rank_antigen'].apply(lambda x: 'no' if x <= cutoff_antigen else 'yes')
        df['check_bcr'] = df['Rank_bcr'].apply(lambda x: 'no' if x <= cutoff_bcr else 'yes')

    df['check'] = df.apply(determine_check, axis=1)

    return df  # Ensure the modified DataFrame is returned


# In[150]:


def split_result(df):
    if 'check' not in df.columns:
        raise ValueError("The input DataFrame must contain a 'check' column.")
    df = df.copy()
    if 'Rank_antigen' not in df.columns:
        df['Rank_antigen'] = pd.NA  # Use NA for missing values
    if 'Rank_bcr' not in df.columns:
        df['Rank_bcr'] = pd.NA
    df_both = df[df['check'] == 'both'].copy()
    df_antigen = df[df['check'] == 'antigen'].copy()
    df_bcr = df[df['check'] == 'bcr'].copy()
    df_pass = df[df['check'] == 'pass'].copy()
    return df_both,df_antigen,df_bcr,df_pass


# In[223]:


def calculate_rank(res_check,antigen_score_dict,bcr_score_dict,antigen_only=False,bcr_only=False):


    if antigen_only and not bcr_only:
        if 'Rank_antigen' not in res_check.columns:
            res_check['Rank_antigen'] = pd.NA  # Keeps missing values consistent
        res_check['Rank_antigen'] = res_check.apply(
            lambda row: locate_rank(row['Score_antigen'], antigen_score_dict[row['Antigen']])
            if pd.notna(row['Score_antigen']) and row['Antigen'] in antigen_score_dict else row['Rank_antigen'],
            axis=1
        )
        return res_check

    elif bcr_only and not antigen_only:
        if 'Rank_bcr' not in res_check.columns:
            res_check['Rank_bcr'] = pd.NA
        res_check['Rank_bcr'] = res_check.apply(
            lambda row: locate_rank(row['Score_bcr'], bcr_score_dict[row['BCR_id']])
            if pd.notna(row['Score_bcr']) and row['BCR_id'] in bcr_score_dict else row['Rank_bcr'],
            axis=1
        )
        return res_check
    elif not antigen_only and not bcr_only:
        if 'Rank_antigen' not in res_check.columns:
            res_check['Rank_antigen'] = pd.NA
        if 'Rank_bcr' not in res_check.columns:
            res_check['Rank_bcr'] = pd.NA
        res_check['Rank_antigen'] = res_check.apply(
            lambda row: locate_rank(row['Score_antigen'], antigen_score_dict[row['Antigen']])
            if pd.notna(row['Score_antigen']) and row['Antigen'] in antigen_score_dict else row['Rank_antigen'],
            axis=1
        )

        res_check['Rank_bcr'] = res_check.apply(
            lambda row: locate_rank(row['Score_bcr'], bcr_score_dict[row['BCR_id']])
            if pd.notna(row['Score_bcr']) and row['BCR_id'] in bcr_score_dict else row['Rank_bcr'],
            axis=1
        )
#         res_check['max_Rank'] = res_check[['Rank_antigen', 'Rank_bcr']].max(axis=1)
#         res_check['ave_Rank'] = res_check[['Rank_antigen', 'Rank_bcr']].mean(axis=1)
        return res_check


# In[224]:


def pip_score_rank(df,antigen_dict,bcr_dict,df_background,antigen_ids,CDR3h_dict,Vh_dict,model1,model2,
                  antigen_score_dict={},bcr_score_dict={},subsample1=100,subsample2=0.03,
                   seed = SEED,antigen_only = False,bcr_only=False):
    sub_antigen_dict,sub_bcr_dict = get_subset_dict(df,antigen_dict,bcr_dict)
#     print(sub_antigen_dict.keys(),sub_bcr_dict.keys())
#     print(antigen_only)
    antigen_score_dict,bcr_score_dict = generate_score_dict(df_background,antigen_ids,sub_antigen_dict,sub_bcr_dict,CDR3h_dict,Vh_dict,BACK_DIR,model1,model2,
                                                                antigen_score_dict=antigen_score_dict,bcr_score_dict=bcr_score_dict,
                                                                subsample1=subsample1/1000000,subsample2=subsample2,seed = SEED,antigen_only = antigen_only,bcr_only=bcr_only)
    df_updated = calculate_rank(df, antigen_score_dict, bcr_score_dict, antigen_only=antigen_only, bcr_only=bcr_only)
    return df_updated#,sub_antigen_dict,antigen_score_dict


# In[152]:


# subsample2=0.03
# subsample1=100


# In[153]:


# cutoff_antigen=cutoffs_dict[subsample1]
# cutoff_bcr=cutoffs_dict2[subsample2*1000]
# res = rank_check(res_check,cutoff_antigen,cutoff_bcr)
# df_both, df_antigen, df_bcr, df_pass = split_result(res)
# dfs_to_concat = [df_pass]
# if not df_both.empty:
#     both_updated=pip_score_rank(df_both,antigen_dict,bcr_dict,background,antigen_ids,CDR3h_dict,Vh_dict,model_mix,model_mix2,
#                   antigen_score_dict=antigen_score_dict,bcr_score_dict=bcr_score_dict,subsample1=subsample1,subsample2=subsample2,
#                    seed = SEED,antigen_only = False,bcr_only=False)
#     dfs_to_concat.append(both_updated)
# if not df_antigen.empty:
#     antigen_updated=pip_score_rank(df_antigen,antigen_dict,bcr_dict,background,antigen_ids,CDR3h_dict,Vh_dict,model_mix,model_mix2,
#                   antigen_score_dict=antigen_score_dict,bcr_score_dict=bcr_score_dict,subsample1=subsample1,subsample2=subsample2,
#                    seed = SEED,antigen_only = True,bcr_only=False)
#     dfs_to_concat.append(antigen_updated)
# if not df_bcr.empty:
#     bcr_updated=pip_score_rank(df_bcr,antigen_dict,bcr_dict,background,antigen_ids,CDR3h_dict,Vh_dict,model_mix,model_mix2,
#                   antigen_score_dict=antigen_score_dict,bcr_score_dict=bcr_score_dict,subsample1=subsample1,subsample2=subsample2,
#                    seed = SEED,antigen_only = False,bcr_only=True)
#     dfs_to_concat.append(bcr_updated)
# if dfs_to_concat:
#     f_res = pd.concat(dfs_to_concat, ignore_index=True)
# else:
#     f_res = pd.DataFrame()
# # subsample1=subsample1*10
# # subsample2=subsample2*10
# print(f_res)


# In[154]:


# antigen_score_dict={}
# bcr_score_dict={}
# cutoff_antigen=cutoffs_dict[subsample1]
# cutoff_bcr=cutoffs_dict2[subsample2*1000]
# f_res = rank_check(res_check,cutoff_antigen,cutoff_bcr)
# while ###MARK HERE
# f_res,antigen_score_dict,bcr_score_dict=one_round_rank(f_res,antigen_score_dict,bcr_score_dict,subsample1,subsample2,seed = SEED,antigen_only = False,bcr_only=False)
# subsample1=subsample1*10
# subsample2=subsample2*10
# f_res=rank_check(f_res,cutoff_antigen,cutoff_bcr)


# In[156]:


# mask_antigen=res_check['Rank_antigen']>=cutoffs_dict[100]
# if mask_antigen.any():
#     output.loc[mask_antigen, ["Rank_antigen"]]=res_check.loc[mask_antigen, ["Rank_antigen"]]
# output


# In[157]:




# def one_round_rank(s_target,f_antigens,output,antigen_score_dict,bcr_score_dict,
#                    subsample1=100,subsample2=0.01,seed =SEED,antigen_only=False,bcr_only=False):
#     columns = ["record_id", "Antigen", "BCR_id", "Score_antigen", "Score_bcr"]
#     output=pd.DataFrame(columns=columns)

#     antigen_score_dict,bcr_score_dict = generate_score_dict(background,antigen_ids,antigen_dict,bcr_dict,CDR3h_dict,Vh_dict,BACK_DIR,model_mix,model_mix2,
#                                                             antigen_score_dict=antigen_score_dict,bcr_score_dict=bcr_score_dict,
#                                                             subsample1=subsample1/1000000,subsample2=subsample2,seed = SEED,
#                                                             antigen_only = antigen_only,bcr_only=bcr_only)
#     res = calculate_rank(s_target,antigen_score_dict,bcr_score_dict,f_antigens,len_dict,model_mix,model_mix2,antigen_only = antigen_only,bcr_only=bcr_only)
#     output[columns]=res[columns]
#     if VERBOSE:
#         print('res after ranking is:')
#         print(res)
#     if antigen_only and not bcr_only:
#         mask_antigen = res['Rank_antigen']>=cutoffs_dict[subsample1]
#         if mask_antigen.any():
#             output.loc[mask_antigen, ["Rank_antigen"]]=res.loc[mask_antigen,["Rank_antigen"]]
#             s_target = s_target.loc[~mask_antigen,]
# #         output = pd.concat([output,res[res['Rank_antigen']>=cutoffs_dict[subsample1]]],axis =0)
# #         f_res = res[res['Rank_antigen']<cutoffs_dict[subsample1]]
#     elif bcr_only and not antigen_only:
#         mask_bcr = res['Rank_bcr']>=cutoffs_dict2[subsample2*1000]
#         if mask_bcr.any():
#             output.loc[mask_bcr, ["Rank_bcr"]]=res.loc[mask_bcr,["Rank_bcr"]]
#             s_target = s_target.loc[~mask_bcr,]
# #         output = pd.concat([output,res[res['Rank_bcr']>=cutoffs_dict2[subsample2*1000]]],axis =0)
# #         f_res = res[res['Rank_bcr']<cutoffs_dict2[subsample2*1000]]
#     else:
#         mask_antigen = res['Rank_antigen']>=cutoffs_dict[subsample1]
#         mask_bcr = res['Rank_bcr']>=cutoffs_dict2[subsample2*1000]
#         if mask_antigen.any():
#             output.loc[mask_antigen, ["Rank_antigen"]]=res.loc[mask_antigen,["Rank_antigen"]]
#         if mask_bcr.any():
#             output.loc[mask_bcr, ["Rank_bcr"]]=res.loc[mask_bcr,["Rank_bcr"]]

# #         output = pd.concat([output,res[(res['Rank_antigen']>=cutoffs_dict[subsample1])&(res['Rank_bcr']>=cutoffs_dict2[subsample2*1000])]],axis =0)
# #         f_res = res[(res['Rank_antigen']<cutoffs_dict[subsample1])|(res['Rank_bcr']<cutoffs_dict2[subsample2*1000])]
#     if VERBOSE:
#         print('The results after this round with',subsample1,'BCRs and',round(subsample2*1115),'Antigens:')
#         print(output)
#     percentage_completed = output.shape[0] / target.shape[0] * 100
#     print(f'Completed {percentage_completed:.2f}% of entries...')

#     f_antigens = {k: v for k, v in f_antigens.items() if k in f_res['Antigen'].values}
#     s_target = s_target[s_target['record_id'].isin(f_res['record_id'])]

#     return(output,f_antigens,s_target,antigen_score_dict,bcr_score_dict)


# In[88]:


# # list(res_check['BCR_id'].unique())
# print(f"{Fore.GREEN}Ranking in {subsample1} background BCRs and/or {round(subsample2 * 1000)} background antigens...{Style.RESET_ALL}")

# print(f"{Fore.LIGHTMAGENTA_EX}Entering BOTH channel:\n Antigen: {list(res_check['Antigen'].unique())} \n BCR: {list(res_check['BCR_id'].unique())}{Style.RESET_ALL}")
# print(f"{Fore.LIGHTCYAN_EX}Entering BOTH channel:\n Antigen: {list(res_check['Antigen'].unique())}{Style.RESET_ALL}")
# print(f"{Fore.LIGHTYELLOW_EX}Entering BCR channel:\n BCR: {list(res_check['BCR_id'].unique())}{Style.RESET_ALL}")


# In[249]:


def one_round_rank(df,antigen_score_dict,bcr_score_dict,subsample1,subsample2,seed = SEED,antigen_only = False,bcr_only=False):
    df_both, df_antigen, df_bcr, df_pass = split_result(df)
    dfs_to_concat = [df_pass]
    if not antigen_only and not bcr_only:
        print(f"{Fore.GREEN}Ranking in {subsample1} background BCRs and/or {round(subsample2 * 1000)} background antigens...{Style.RESET_ALL}")
        if not df_both.empty:
            print(f"{Fore.LIGHTMAGENTA_EX}Entering BOTH channel:\n Antigen: {list(df_both['Antigen'].unique())} \n BCR: {list(df_both['BCR_id'].unique())}{Style.RESET_ALL}")
            both_updated=pip_score_rank(df_both,antigen_dict,bcr_dict,background,antigen_ids,CDR3h_dict,Vh_dict,model_mix,model_mix2,
                          antigen_score_dict=antigen_score_dict,bcr_score_dict=bcr_score_dict,subsample1=subsample1,subsample2=subsample2,
                           seed = SEED,antigen_only = False,bcr_only=False)
            dfs_to_concat.append(both_updated)
        if not df_antigen.empty:
            print(f"{Fore.LIGHTCYAN_EX}Entering ANTIGEN channel:\n Antigen: {list(df_antigen['Antigen'].unique())}{Style.RESET_ALL}")
            antigen_updated=pip_score_rank(df_antigen,antigen_dict,bcr_dict,background,antigen_ids,CDR3h_dict,Vh_dict,model_mix,model_mix2,
                          antigen_score_dict=antigen_score_dict,bcr_score_dict=bcr_score_dict,subsample1=subsample1,subsample2=subsample2,
                           seed = SEED,antigen_only = True,bcr_only=False)
            dfs_to_concat.append(antigen_updated)
        if not df_bcr.empty:
            print(f"{Fore.LIGHTYELLOW_EX}Entering BCR channel:\n BCR: {list(df_bcr['BCR_id'].unique())}{Style.RESET_ALL}")
            bcr_updated=pip_score_rank(df_bcr,antigen_dict,bcr_dict,background,antigen_ids,CDR3h_dict,Vh_dict,model_mix,model_mix2,
                          antigen_score_dict=antigen_score_dict,bcr_score_dict=bcr_score_dict,subsample1=subsample1,subsample2=subsample2,
                           seed = SEED,antigen_only = False,bcr_only=True)
            dfs_to_concat.append(bcr_updated)
    elif antigen_only:
        print(f"{Fore.GREEN}Ranking in {subsample1} background BCRs...{Style.RESET_ALL}")
        if not df_both.empty:
            print(f"{Fore.LIGHTMAGENTA_EX}Entering BOTH channel:\n Antigen: {list(df_both['Antigen'].unique())} \n BCR: {list(df_both['BCR_id'].unique())}{Style.RESET_ALL}")
            both_updated=pip_score_rank(df_both,antigen_dict,bcr_dict,background,antigen_ids,CDR3h_dict,Vh_dict,model_mix,model_mix2,
                          antigen_score_dict=antigen_score_dict,bcr_score_dict=bcr_score_dict,subsample1=subsample1,subsample2=subsample2,
                           seed = SEED,antigen_only = True,bcr_only=False)
            dfs_to_concat.append(both_updated)
        if not df_antigen.empty:
            print(f"{Fore.LIGHTCYAN_EX}Entering ANTIGEN channel:\n Antigen: {list(df_antigen['Antigen'].unique())}{Style.RESET_ALL}")
            antigen_updated=pip_score_rank(df_antigen,antigen_dict,bcr_dict,background,antigen_ids,CDR3h_dict,Vh_dict,model_mix,model_mix2,
                          antigen_score_dict=antigen_score_dict,bcr_score_dict=bcr_score_dict,subsample1=subsample1,subsample2=subsample2,
                           seed = SEED,antigen_only = True,bcr_only=False)
            dfs_to_concat.append(antigen_updated)
    elif bcr_only:
        print(f"{Fore.GREEN}Ranking in {round(subsample2 * 1000)} background antigens...{Style.RESET_ALL}")
        if not df_both.empty:
            print(f"{Fore.LIGHTMAGENTA_EX}Entering BOTH channel:\n Antigen: {list(df_both['Antigen'].unique())} \n BCR: {list(df_both['BCR_id'].unique())}{Style.RESET_ALL}")
            both_updated=pip_score_rank(df_both,antigen_dict,bcr_dict,background,antigen_ids,CDR3h_dict,Vh_dict,model_mix,model_mix2,
                          antigen_score_dict=antigen_score_dict,bcr_score_dict=bcr_score_dict,subsample1=subsample1,subsample2=subsample2,
                           seed = SEED,antigen_only = False,bcr_only=True)
            dfs_to_concat.append(both_updated)
        if not df_bcr.empty:
            print(f"{Fore.LIGHTYELLOW_EX}Entering BCR channel:\n BCR: {list(df_bcr['BCR_id'].unique())}{Style.RESET_ALL}")
            bcr_updated=pip_score_rank(df_bcr,antigen_dict,bcr_dict,background,antigen_ids,CDR3h_dict,Vh_dict,model_mix,model_mix2,
                          antigen_score_dict=antigen_score_dict,bcr_score_dict=bcr_score_dict,subsample1=subsample1,subsample2=subsample2,
                           seed = SEED,antigen_only = False,bcr_only=True)
            dfs_to_concat.append(bcr_updated)

    if dfs_to_concat:
        f_res = pd.concat(dfs_to_concat, ignore_index=True)
        f_res = f_res.dropna(axis=1,how='all')
    else:
        f_res = pd.DataFrame()
    return f_res,antigen_score_dict,bcr_score_dict


# In[246]:





# In[159]:





target_file = pd.read_csv(INPUT) # required columns 'Vh','CDR3h', optional 'species'
background = pd.read_csv(BACKGROUND,compression='gzip') # with columns 'Vh','CDR3h','file','species'




# In[ ]:






# ###MARK HERE TO EDIT ####stea
# file_path = '/project/DPDS/Wang_lab/shared/BCR_antigen/data/rfscript/Antigen_embed/CLEANED/NPY/antigen_list2.txt'

# # Open and read the file
# with open(file_path, 'r') as file:
#     lines = file.readlines()

# # Strip newline characters from each line
# antigen_ids = [line.strip() for line in lines]
# target_file['Antigen_id']=['Hie/H4Hubei','Hie/H4Hubei','Chen/ova','Chen/ova']


# In[ ]:






# if 'BCR_species' in target_file.columns:
#     unique_values = target_file['BCR_species'].dropna().unique()
#     if len(unique_values) == 1 and unique_values[0] in {'human', 'mouse'}:
#         species = unique_values[0]
#         print('The species is:',species)
#         if MATCHING_SPECIES:
#             print('matching the background species...')
#             background = background[background['species']==species]
#     else:
#         print('The target species are not unique or not in the background species.')
# else:
#     print('No species are specified!')


# In[160]:






target = check_npy(target_file)
# target = preprocess(target_file)
# target.to_csv(INPUT_DIR+'/filtered_input.csv')


# In[161]:




bcr_dict = get_bcr_dict(target)
for key,value in bcr_dict.items():
    print(key)
antigen_dict = get_antigen_dict(target,NPY_DIR)
for key,value in antigen_dict.items():
    print(key,value.shape)


# In[162]:





cutoffs_dict = pd.read_csv('paras/cutoff_table.txt', sep='\t').set_index('backsample')['cutoffs'].to_dict()
print('threshold to enter the next level of ranking for background bcrs:',cutoffs_dict)
cutoffs_dict2 = pd.read_csv('paras/cutoff_table2.txt', sep='\t').set_index('backsample')['cutoffs'].to_dict()
print('threshold to enter the next level of ranking for background antigens:',cutoffs_dict2)



# In[163]:






len_dict = target.drop_duplicates(subset='Antigen_id').set_index('Antigen_id')['Antigen_seq'].map(len).to_dict()
# len_dict={}
# for key, value in antigen_dict.items():
#     len_dict[key]=value.shape[1]


# In[164]:






model_mix = mix_model()
checkpoint = torch.load(MODEL)
model_mix.load_state_dict(checkpoint)
model_mix.to(device)
optimizer = torch.optim.Adam(model_mix.parameters(),lr=LR)
scheduler = ReduceLROnPlateau(optimizer,mode = 'min',factor=0.1,patience =10,verbose=False,threshold_mode='rel',threshold=1e-5)


# In[165]:






model_mix2 = mix_model()
checkpoint2 = torch.load(MODEL2)
model_mix2.load_state_dict(checkpoint2)
model_mix2.to(device)
optimizer2 = torch.optim.Adam(model_mix2.parameters(),lr=LR)
scheduler2 = ReduceLROnPlateau(optimizer2,mode = 'min',factor=0.1,patience =10,verbose=False,threshold_mode='rel',threshold=1e-5)


# In[ ]:






# score_dict = generate_score_dict(background,score_dict,antigen_dict,CDR3h_dict,Vh_dict,model_mix)
# res = calculate_rank(target,score_dict,antigen_dict,len_dict,model_mix)
# for antigen_id, antigen in antigen_dict.items():
#     background_loader = DataLoader(rankDataset(background, antigen, CDR3h_dict, Vh_dict,subsample_ratio=1/10000,seed=SEED),1)
#     score_background = background_scores(background_loader,model_mix)
#     score_dict[antigen_id]=score_background


# In[166]:






antigen_score_dict = {}
bcr_score_dict = {}

if BOTTOMLINE1==10000 and SEED == 1:
    with open('data/background/back1w_V_dict.pkl','rb') as f:
        Vh_dict = pickle.load(f)
    with open('data/background/back1w_CDR3_dict.pkl','rb') as f:
        CDR3h_dict = pickle.load(f)
else:
    back100 = background.sample(frac=1/10000, random_state=SEED)
    Vh_dict = build_BCR_dict(back100,'Vh',precise = True)
    CDR3h_dict = build_BCR_dict(back100,'CDR3h',precise = True)


if BOTTOMLINE2==0.3 and SEED == 1:
    with open("data/background/default300.txt", "r") as file:
        antigen_ids = [line.strip() for line in file.readlines()]
else:
    antigen_files = glob(BACK_DIR+'/*/*pair.npy')
    antigen_npys = [antigen.split('/')[-2]+'/'+antigen.split('/')[-1] for antigen in antigen_files]
    # antigen_files[0].split('/')[-2]+'/'+antigen_files[0].split('/')[-1]
    antigen_ids = [antigen_id.replace('.pair.npy','') for antigen_id in antigen_npys]
# antigen_score_dict,bcr_score_dict = generate_score_dict(background,antigen_ids,antigen_dict,CDR3h_dict,Vh_dict,NPY_DIR,model_mix,model_mix2,
#                                  # antigen_score_dict=antigen_score_dict,bcr_score_dict=bcr_score_dict,subsample1=1/10000,subsample2=0.1,seed =SEED)
# if DEBUG:
#     with open(OUT_DIR+'/first100_V_dict.pkl','wb') as f:
#         pickle.dump(Vh_dict, f)
#     with open(OUT_DIR+'/first100_CDR3_dict.pkl','wb') as f:
#         pickle.dump(CDR3h_dict, f)


# In[ ]:




# subsample1 = 100
# subsample2 = 0.01
# antigen_score_dict,bcr_score_dict = generate_score_dict(background,antigen_ids,antigen_dict,bcr_dict,CDR3h_dict,Vh_dict,NPY_DIR,model_mix,model_mix2,
#                                                             antigen_score_dict=antigen_score_dict,bcr_score_dict=bcr_score_dict,
#                                                             subsample1=subsample1/1000000,subsample2=subsample2,seed = SEED)


# In[50]:






if args.no_rank:
    check_loader = checkDataset(target, antigen_dict, NPY_DIR,len_dict)
    res_check = check_score(check_loader,model_mix,model_mix2)
    if not args.no_merge:
        merged = target.merge(res_check, on=['record_id','BCR_id'], how='left')
        merged.to_csv(OUT_DIR+'/merged_results_no_rank.csv',index=False)
        exit()
    else:
        res_check.to_csv(OUT_DIR+'/binding_results_no_rank.csv',index=False)
        exit()


if args.export_background:
    antigen_score_dict,bcr_score_dict = {},{}
    # if BOTTOMLINE1==10000 & SEED == 1:
    #     with open('data/background/back1w_V_dict.pkl','rb') as f:
    #         Vh_dict = pickle.load(f)
    #     with open('data/background/back1w_CDR3_dict.pkl','rb') as f:
    #         CDR3h_dict = pickle.load(f)
    if not ANTIGEN_ONLY and not BCR_ONLY:
        antigen_score_dict,bcr_score_dict = generate_score_dict(background,antigen_ids,antigen_dict,bcr_dict,CDR3h_dict,Vh_dict,BACK_DIR,model_mix,model_mix2,
                                                                antigen_score_dict=antigen_score_dict,bcr_score_dict=bcr_score_dict,
                                                                subsample1=subsample1/1000000,subsample2=subsample2,seed = SEED,antigen_only = False,bcr_only=False)
        print('Exporting dict for '+str(BOTTOMLINE1)+' background BCRs...')
        with open(OUT_DIR+'/background_antigen_score_dict.pkl', 'wb') as outfile:
            pickle.dump(antigen_score_dict, outfile)
        print('Exporting dict for '+str(BOTTOMLINE2*1148)+' background antigens...')
        with open(OUT_DIR+'/background_bcr_score_dict.pkl', 'wb') as outfile:
            pickle.dump(bcr_score_dict, outfile)
        exit()
    if ANTIGEN_ONLY and not BCR_ONLY:
        antigen_score_dict,bcr_score_dict = generate_score_dict(background,antigen_ids,antigen_dict,bcr_dict,CDR3h_dict,Vh_dict,BACK_DIR,model_mix,model_mix2,
                                                                antigen_score_dict=antigen_score_dict,bcr_score_dict=bcr_score_dict,
                                                                subsample1=subsample1/1000000,subsample2=subsample2,seed = SEED,antigen_only = True,bcr_only=False)
        print('Exporting dict for '+str(BOTTOMLINE1)+' background BCRs...')
        with open(OUT_DIR+'/background_antigen_score_dict.pkl', 'wb') as outfile:
            pickle.dump(antigen_score_dict, outfile)
        exit()
    elif BCR_ONLY and not ANTIGEN_ONLY:
        antigen_score_dict,bcr_score_dict = generate_score_dict(background,antigen_ids,antigen_dict,bcr_dict,CDR3h_dict,Vh_dict,BACK_DIR,model_mix,model_mix2,
                                                                antigen_score_dict=antigen_score_dict,bcr_score_dict=bcr_score_dict,
                                                                subsample1=subsample1/1000000,subsample2=subsample2,seed = SEED,antigen_only = False,bcr_only=True)
        print('Exporting dict for '+str(BOTTOMLINE2*1148)+' background antigens...')
        with open(OUT_DIR+'/background_bcr_score_dict.pkl', 'wb') as outfile:
            pickle.dump(bcr_score_dict, outfile)
        exit()


    else:
        print('Exporting dict for '+str(BOTTOMLINE1)+' background BCRs and '+str(BOTTOMLINE2*1148)+' background antigens...')
        with open(OUT_DIR+'/background_antigen_score_dict.pkl', 'wb') as outfile:
            pickle.dump(antigen_score_dict, outfile)
        with open(OUT_DIR+'/background_bcr_score_dict.pkl', 'wb') as outfile:
            pickle.dump(bcr_score_dict, outfile)
        exit()





# In[169]:




# BOTTOMLINE1=1000
# BOTTOMLINE2=0.3
# # DEBUG = True
# # print(BOTTOMLINE1,BOTTOMLINE2)


# In[ ]:




# antigen_score_dict,bcr_score_dict = {},{}
# if BOTTOMLINE1==10000 & SEED == 1:
#     with open('data/background/back1w_V_dict.pkl','rb') as f:
#         Vh_dict = pickle.load(f)
#     with open('data/background/back1w_CDR3_dict.pkl','rb') as f:
#         CDR3h_dict = pickle.load(f)
# antigen_score_dict,bcr_score_dict = generate_score_dict(background,antigen_ids,antigen_dict,bcr_dict,CDR3h_dict,Vh_dict,NPY_DIR,model_mix,model_mix2,
#                              antigen_score_dict={},bcr_score_dict={},subsample1=BOTTOMLINE1/1000000,subsample2=BOTTOMLINE2,seed =SEED,
#                              antigen_only = ANTIGEN_ONLY,bcr_only=BCR_ONLY)

# print('Exporting dict for '+str(BOTTOMLINE1)+' background BCRs and '+str(BOTTOMLINE2*1150)+' background antigens...')
# with open(OUT_DIR+'/background_antigen_score_dict.pkl', 'wb') as outfile:
#     pickle.dump(antigen_score_dict, outfile)
# with open(OUT_DIR+'/background_bcr_score_dict.pkl', 'wb') as outfile:
#     pickle.dump(bcr_score_dict, outfile)


# In[ ]:





# antigen_score_dict
# bcr_score_dict


# In[ ]:






# subsample1 = SUBSAMPLE1
# subsample2 = SUBSAMPLE2
# s_target =  target
# f_antigens = antigen_dict
# output = pd.DataFrame()
# antigen_score_dict = {}
# bcr_score_dict = {}
# output,f_res,f_antigens,s_target,antigen_score_dict,bcr_score_dict=one_round_rank(s_target,f_antigens,output,antigen_score_dict,bcr_score_dict,subsample1=subsample1,subsample2=subsample2,seed=SEED)


# In[ ]:




# BOTTOMLINE1 = 1000
# BOTTOMLINE2 = 0.1


# In[170]:




print('subsample1:',SUBSAMPLE1,'subsample2:',SUBSAMPLE2,'BOTTOMLINE1:',BOTTOMLINE1,'BOTTOMLINE2:',BOTTOMLINE2)


# In[ ]:




# subsample1 = SUBSAMPLE1
# subsample2 = SUBSAMPLE2
# print('the initial subsample for background BCRs is:',subsample1)
# print('the initial subsample for background antigens is:',subsample2)
# s_target =  target
# f_antigens = antigen_dict
# output = pd.DataFrame()
# antigen_score_dict = {}
# bcr_score_dict = {}


# In[ ]:




# start_time=time.time()
# print('Ranking in',subsample1,'background BCRs and/or',str(round(subsample2*1000)),'background antigens...')
# # antigen_score_dict,bcr_score_dict = generate_score_dict(background,antigen_ids,antigen_dict,CDR3h_dict,Vh_dict,NPY_DIR,model_mix,model_mix2,
# #                              antigen_score_dict=antigen_score_dict,bcr_score_dict=bcr_score_dict,subsample1=subsample1/1000000,subsample2=subsample2,seed =SEED)
# # res = calculate_rank(s_target,antigen_score_dict,bcr_score_dict,f_antigens,len_dict,model_mix,model_mix2)
# # output = pd.concat([output,res[(res['Rank_antigen']>=cutoffs_dict[subsample])&(res['Rank_bcr']>=cutoffs_dict[subsample])]],axis =0)
# # percentage_completed = output.shape[0] / target.shape[0] * 100
# # print(f'Completed {percentage_completed:.2f}% of entries...')
# # f_res = res[(res['Rank_antigen']<cutoffs_dict[subsample])|(res['Rank_bcr']<cutoffs_dict[subsample])]
# # f_antigens = {k: v for k, v in f_antigens.items() if k in f_res['Antigen'].values}
# # s_target = s_target[s_target['record_id'].isin(f_res['record_id'])]
# # output,f_res,f_antigens,s_target,antigen_score_dict,bcr_score_dict=one_round_rank(s_target,f_antigens,output,antigen_score_dict,bcr_score_dict,subsample1=subsample1,subsample2=subsample2,seed=SEED)
# output,f_res,f_antigens,s_target,antigen_score_dict,bcr_score_dict=one_round_rank(s_target,f_antigens,output,antigen_score_dict,bcr_score_dict,
#                                                                                   subsample1=subsample1,subsample2=subsample2,seed=SEED,
#                                                                                   antigen_only=ANTIGEN_ONLY,bcr_only=BCR_ONLY)
# print('Processing time for',str(subsample1),'BCRs and/or',str(subsample2*1150),'Antigens is:',str(time.time()-start_time))
# subsample1 = subsample1*10
# subsample2 = subsample2*10


# In[ ]:




# antigen_score_dict,bcr_score_dict = generate_score_dict(background,antigen_ids,antigen_dict,bcr_dict,CDR3h_dict,Vh_dict,BACK_DIR,model_mix,model_mix2,
#                                                         antigen_score_dict=antigen_score_dict,bcr_score_dict=bcr_score_dict,
#                                                         subsample1=subsample1/1000000,subsample2=subsample2,seed = SEED,
#                                                         antigen_only = False,bcr_only=True)


# In[227]:


# f_res


# In[255]:


# BCR_ONLY = True
# ANTIGEN_ONLY = False
# f_res
# # f_res =  rank_check(res_check,0,0,antigen_only=ANTIGEN_ONLY)
# # i_res,antigen_score_dict,bcr_score_dict=one_round_rank(f_res,antigen_score_dict,bcr_score_dict,subsample1,subsample2,seed = SEED,antigen_only = True)


# In[229]:


res_check=get_check_score(target,antigen_dict,len_dict,model_mix,model_mix2)


# In[256]:


subsample1 = SUBSAMPLE1
subsample2 = SUBSAMPLE2
print('the initial subsample for background BCRs is:',subsample1)
print('the initial subsample for background antigens is:',subsample2)
f_res =  rank_check(res_check,0,0,antigen_only=ANTIGEN_ONLY,bcr_only =BCR_ONLY)
print(f_res.head())
if not ANTIGEN_ONLY and not BCR_ONLY:
    while not (f_res['check']=='pass').all():
        if subsample1<=BOTTOMLINE1 and subsample2 <= BOTTOMLINE2:
            start_time=time.time()
            i_res,antigen_score_dict,bcr_score_dict=one_round_rank(f_res,antigen_score_dict,bcr_score_dict,subsample1,subsample2,seed = SEED,antigen_only = False,bcr_only=False)
            print(f"{Fore.BLUE}The initial ranks:\n{i_res.loc[:, ~i_res.columns.str.contains('check', case=False)]}{Style.RESET_ALL}")
            cutoff_antigen=cutoffs_dict[subsample1]
            cutoff_bcr=cutoffs_dict2[subsample2*1000]
            print(f"{Fore.LIGHTWHITE_EX}After checking with the cutoff for antigen {Fore.RED}{cutoff_antigen}{Fore.LIGHTWHITE_EX} "
                  f"and the cutoff for BCR: {Fore.RED}{cutoff_bcr}{Style.RESET_ALL}")
            f_res=rank_check(i_res,cutoff_antigen,cutoff_bcr)
            print(f"{Fore.BLUE}{f_res}{Style.RESET_ALL}")
            print(f"{Fore.LIGHTWHITE_EX}Processing time for {round(subsample1)} BCRs and/or {round(subsample2 * 1148)} Antigens is: {time.time() - start_time:.2f} seconds{Style.RESET_ALL}")
            percentage_completed = (f_res['check'] == 'pass').mean() * 100
            print(f'Completed {percentage_completed:.2f}% of entries...')
            # Print the sample size update message
            print(f"{Fore.YELLOW}Sample size x10...{Style.RESET_ALL}")
            subsample1 *= 10
            subsample2 *= 10
        elif subsample1<=BOTTOMLINE1 and subsample2 > BOTTOMLINE2:
            start_time=time.time()
            i_res,antigen_score_dict,bcr_score_dict=one_round_rank(f_res,antigen_score_dict,bcr_score_dict,subsample1,subsample2,seed = SEED,antigen_only = True,bcr_only=False)
            print(f"{Fore.BLUE}The initial ranks:\n{i_res.loc[:, ~i_res.columns.str.contains('check', case=False)]}{Style.RESET_ALL}")
            cutoff_antigen=cutoffs_dict[subsample1]
            cutoff_bcr=cutoffs_dict2[subsample2*1000]
            print(f"{Fore.LIGHTWHITE_EX}After checking with the cutoff for antigen {Fore.RED}{cutoff_antigen}{Fore.LIGHTWHITE_EX} "
                  f"and the cutoff for BCR: {Fore.RED}{cutoff_bcr}{Style.RESET_ALL}")
            f_res=rank_check(i_res,cutoff_antigen,cutoff_bcr)
            print(f"{Fore.BLUE}{f_res}{Style.RESET_ALL}")
            print(f"{Fore.LIGHTWHITE_EX}Processing time for {round(subsample1)} BCRs is: {time.time() - start_time:.2f} seconds{Style.RESET_ALL}")
            percentage_completed = (f_res['check'] == 'pass').mean() * 100
            print(f'Completed {percentage_completed:.2f}% of entries...')
            # Print the sample size update message
            print(f"{Fore.YELLOW}Sample size x10...{Style.RESET_ALL}")
            subsample1 *= 10
        elif subsample1>BOTTOMLINE1 and subsample2 <= BOTTOMLINE2:
            start_time=time.time()
            i_res,antigen_score_dict,bcr_score_dict=one_round_rank(f_res,antigen_score_dict,bcr_score_dict,subsample1,subsample2,seed = SEED,antigen_only = False,bcr_only=True)
            print(f"{Fore.BLUE}The initial ranks:\n{i_res.loc[:, ~i_res.columns.str.contains('check', case=False)]}{Style.RESET_ALL}")
            cutoff_antigen=cutoffs_dict[subsample1]
            cutoff_bcr=cutoffs_dict2[subsample2*1000]
            print(f"{Fore.LIGHTWHITE_EX}After checking with the cutoff for antigen {Fore.RED}{cutoff_antigen}{Fore.LIGHTWHITE_EX} "
                  f"and the cutoff for BCR: {Fore.RED}{cutoff_bcr}{Style.RESET_ALL}")
            f_res=rank_check(i_res,cutoff_antigen,cutoff_bcr)
            print(f"{Fore.BLUE}{f_res}{Style.RESET_ALL}")
            print(f"{Fore.LIGHTWHITE_EX}Processing time for {round(subsample2 * 1148)} Antigens is: {time.time() - start_time:.2f} seconds{Style.RESET_ALL}")
            percentage_completed = (f_res['check'] == 'pass').mean() * 100
            print(f'Completed {percentage_completed:.2f}% of entries...')
            # Print the sample size update message
            print(f"{Fore.YELLOW}Sample size x10...{Style.RESET_ALL}")
            subsample2 *= 10
        else:
            print(f"{Fore.RED}Both subsample1 and subsample2 are above their bottom limits. Exiting loop.{Style.RESET_ALL}")
            break
    f_res['max_Rank'] = f_res[['Rank_antigen', 'Rank_bcr']].max(axis=1)
    f_res['ave_Rank'] = f_res[['Rank_antigen', 'Rank_bcr']].mean(axis=1)
    output=f_res.loc[:, ~f_res.columns.str.contains('check', case=False)]
    print(output.head())
    with open(OUT_DIR+'/background_antigen_score_dict_'+str(BOTTOMLINE1)+'BCRs.pkl', 'wb') as outfile:
        pickle.dump(antigen_score_dict, outfile)
    with open(OUT_DIR+'/background_bcr_score_dict_'+str(BOTTOMLINE2)+'Antigens.pkl', 'wb') as outfile:
        pickle.dump(bcr_score_dict, outfile)
elif ANTIGEN_ONLY and not BCR_ONLY:
    while not (f_res['check']=='pass').all():
        if subsample1<=BOTTOMLINE1:
            start_time=time.time()
            i_res,antigen_score_dict,bcr_score_dict=one_round_rank(f_res,antigen_score_dict,bcr_score_dict,subsample1,subsample2,seed = SEED,antigen_only = True)
            print(f"{Fore.BLUE}The initial ranks:\n{i_res.loc[:, ~i_res.columns.str.contains('check', case=False)]}{Style.RESET_ALL}")
            cutoff_antigen=cutoffs_dict[subsample1]
#             cutoff_bcr=cutoffs_dict2[subsample2*1000]
            print(f"{Fore.LIGHTWHITE_EX}After checking with the cutoff for antigen {Fore.RED}{cutoff_antigen}{Style.RESET_ALL}")
            f_res=rank_check(i_res,cutoff_antigen,cutoff_bcr,antigen_only=True)
            print(f"{Fore.BLUE}{f_res}{Style.RESET_ALL}")
            print(f"{Fore.LIGHTWHITE_EX}Processing time for {round(subsample1)} BCRs is: {time.time() - start_time:.2f} seconds{Style.RESET_ALL}")
            percentage_completed = (f_res['check'] == 'pass').mean() * 100
            print(f'Completed {percentage_completed:.2f}% of entries...')
            # Print the sample size update message
            print(f"{Fore.YELLOW}Sample size x10...{Style.RESET_ALL}")
            subsample1 *= 10
        else:
            print(f"{Fore.RED}Subsample1{subsample1} is above the bottom limit{BOTTOMLINE1}. Exiting loop.{Style.RESET_ALL}")
            break
    output=f_res.loc[:, ~f_res.columns.str.contains('check', case=False)]
    print(output.head())
    with open(OUT_DIR+'/background_antigen_score_dict_'+str(BOTTOMLINE1)+'BCRs.pkl', 'wb') as outfile:
        pickle.dump(antigen_score_dict, outfile)
elif not ANTIGEN_ONLY and BCR_ONLY:
    while not (f_res['check']=='pass').all():
        if subsample2 <= BOTTOMLINE2:
            start_time=time.time()
            i_res,antigen_score_dict,bcr_score_dict=one_round_rank(f_res,antigen_score_dict,bcr_score_dict,subsample1,subsample2,seed = SEED,bcr_only = True)
            print(f"{Fore.BLUE}The initial ranks:\n{i_res.loc[:, ~i_res.columns.str.contains('check', case=False)]}{Style.RESET_ALL}")
#             cutoff_antigen=cutoffs_dict[subsample1]
            cutoff_bcr=cutoffs_dict2[subsample2*1000]
            print(f"{Fore.LIGHTWHITE_EX}After checking with the cutoff for BCR: {Fore.RED}{cutoff_bcr}{Style.RESET_ALL}")
            f_res=rank_check(i_res,cutoff_antigen,cutoff_bcr,bcr_only=True)
            print(f"{Fore.BLUE}{f_res}{Style.RESET_ALL}")
            print(f"{Fore.LIGHTWHITE_EX}Processing time for {round(subsample2 * 1148)} Antigens is: {time.time() - start_time:.2f} seconds{Style.RESET_ALL}")
            percentage_completed = (f_res['check'] == 'pass').mean() * 100
            print(f'Completed {percentage_completed:.2f}% of entries...')
            # Print the sample size update message
            print(f"{Fore.YELLOW}Sample size x10...{Style.RESET_ALL}")
            subsample2 *= 10
        else:
            print(f"{Fore.RED}Subsample2 {subsample2} is above the bottom limit {BOTTOMLINE2}. Exiting loop.{Style.RESET_ALL}")
            break
    output=f_res.loc[:, ~f_res.columns.str.contains('check', case=False)]
    print(output.head())
    with open(OUT_DIR+'/background_bcr_score_dict_'+str(BOTTOMLINE2)+'Antigens.pkl', 'wb') as outfile:
        pickle.dump(bcr_score_dict, outfile)


# In[ ]:





# In[ ]:





# In[123]:


if args.no_merge:
    output.to_csv(OUT_DIR+'/binding_results.csv',index=False)
else:
    merged = target.merge(output, on=['record_id','BCR_id'], how='left')
    merged.to_csv(OUT_DIR+'/merged_results.csv',index=False)


# In[59]:


# subsample1 = SUBSAMPLE1
# subsample2 = SUBSAMPLE2
# print('the initial subsample rer background BCRs is:',subsample1)
# print('the initial subsample for background antigens is:',subsample2)
# f_res =  res_check
# # f_antigens = antigen_dict
# # output = pd.DataFrame()
# antigen_score_dict = {}
# bcr_score_dict = {}
# while not (f_res['check']=='pass').all():
#     if subsample1<=BOTTOMLINE1 and subsample2 <= BOTTOMLINE2:
#         start_time=time.time()
#         print('Ranking in',subsample1,'background BCRs and/or',str(round(subsample2*1000)),'background antigens...')
#         # antigen_score_dict,bcr_score_dict = generate_score_dict(background,antigen_ids,antigen_dict,CDR3h_dict,Vh_dict,NPY_DIR,model_mix,model_mix2,
#         #                              antigen_score_dict=antigen_score_dict,bcr_score_dict=bcr_score_dict,subsample1=subsample1/1000000,subsample2=subsample2,seed =SEED)
#         # res = calculate_rank(s_target,antigen_score_dict,bcr_score_dict,f_antigens,len_dict,model_mix,model_mix2)
#         # output = pd.concat([output,res[(res['Rank_antigen']>=cutoffs_dict[subsample])&(res['Rank_bcr']>=cutoffs_dict[subsample])]],axis =0)
#         # percentage_completed = output.shape[0] / target.shape[0] * 100
#         # print(f'Completed {percentage_completed:.2f}% of entries...')
#         # f_res = res[(res['Rank_antigen']<cutoffs_dict[subsample])|(res['Rank_bcr']<cutoffs_dict[subsample])]
#         # f_antigens = {k: v for k, v in f_antigens.items() if k in f_res['Antigen'].values}
#         # s_target = s_target[s_target['record_id'].isin(f_res['record_id'])]
#         # output,f_res,f_antigens,s_target,antigen_score_dict,bcr_score_dict=one_round_rank(s_target,f_antigens,output,antigen_score_dict,bcr_score_dict,subsample1=subsample1,subsample2=subsample2,seed=SEED)
#         output,f_res,f_antigens,s_target,antigen_score_dict,bcr_score_dict=one_round_rank(s_target,f_antigens,output,antigen_score_dict,bcr_score_dict,
#                                                                                           subsample1=subsample1,subsample2=subsample2,seed=SEED,
#                                                                                           antigen_only=ANTIGEN_ONLY,bcr_only=BCR_ONLY)
#         print('Processing time for',str(subsample1),'BCRs and/or',str(subsample2*1148),'Antigens is:',str(time.time()-start_time))
#         subsample1 = subsample1*10
#         subsample2 = subsample2*10
# #         print(subsample1,subsample2)
#     elif subsample1 > BOTTOMLINE1 and (subsample2*1000) <= BOTTOMLINE2:
#         while (subsample2*1148) <= BOTTOMLINE2 and len(s_target) > 0:
#             start_time=time.time()

#             print('Ranking in',str(round(subsample2*1000)),'background antigens...')
#             output,f_res,f_antigens,s_target,antigen_score_dict,bcr_score_dict=one_round_rank(s_target,f_antigens,output,antigen_score_dict,bcr_score_dict,
#                                                                                               subsample1=subsample1,subsample2=subsample2,seed=SEED,bcr_only=True)
#             # antigen_score_dict,bcr_score_dict = generate_score_dict(background,antigen_ids,antigen_dict,CDR3h_dict,Vh_dict,NPY_DIR,model_mix,model_mix2,
#             #                              antigen_score_dict=antigen_score_dict,bcr_score_dict=bcr_score_dict,subsample1=subsample1/1000000,subsample2=subsample2,seed =SEED,bcr_only=True)
#             # res = calculate_rank(s_target,antigen_score_dict,bcr_score_dict,f_antigens,len_dict,model_mix,model_mix2)
#             # output = pd.concat([output,res[(res['Rank_antigen']>=cutoffs_dict[subsample])&(res['Rank_bcr']>=cutoffs_dict[subsample])]],axis =0)
#             # percentage_completed = output.shape[0] / target.shape[0] * 100
#             # print(f'Completed {percentage_completed:.2f}% of entries...')
#             # f_res = res[(res['Rank_antigen']<cutoffs_dict[subsample])|(res['Rank_bcr']<cutoffs_dict[subsample])]
#             # f_antigens = {k: v for k, v in f_antigens.items() if k in f_res['Antigen'].values}
#             # s_target = s_target[s_target['record_id'].isin(f_res['record_id'])]
#             print('Processing time for',str(subsample2*1148),'Antigens is:',str(time.time()-start_time))
#             subsample2 = subsample2*10
#     elif subsample1 <= BOTTOMLINE1 and (subsample2*1000) > BOTTOMLINE2:
#         while subsample1 <= BOTTOMLINE1 and len(s_target) > 0:
#             print('Ranking in',subsample1,'background BCRs...')
#             start_time=time.time()

#             # antigen_score_dict,bcr_score_dict = generate_score_dict(background,antigen_ids,antigen_dict,CDR3h_dict,Vh_dict,NPY_DIR,model_mix,model_mix2,
#             #                              antigen_score_dict=antigen_score_dict,bcr_score_dict=bcr_score_dict,subsample1=subsample1/1000000,subsample2=subsample2,seed =SEED,antigen_only=True)
#             # res = calculate_rank(s_target,antigen_score_dict,bcr_score_dict,f_antigens,len_dict,model_mix,model_mix2)
#             # output = pd.concat([output,res[(res['Rank_antigen']>=cutoffs_dict[subsample])&(res['Rank_bcr']>=cutoffs_dict[subsample])]],axis =0)
#             # percentage_completed = output.shape[0] / target.shape[0] * 100
#             # print(f'Completed {percentage_completed:.2f}% of entries...')
#             # f_res = res[(res['Rank_antigen']<cutoffs_dict[subsample])|(res['Rank_bcr']<cutoffs_dict[subsample])]
#             # f_antigens = {k: v for k, v in f_antigens.items() if k in f_res['Antigen'].values}
#             # s_target = s_target[s_target['record_id'].isin(f_res['record_id'])]
#             output,f_res,f_antigens,s_target,antigen_score_dict,bcr_score_dict=one_round_rank(s_target,f_antigens,output,antigen_score_dict,bcr_score_dict,
#                                                                                               subsample1=subsample1,subsample2=subsample2,seed=SEED,antigen_only=True)
#             print('Processing time for',str(subsample1),'BCRs is:',str(time.time()-start_time))
#             subsample1 = subsample1*10
#     elif subsample1 > BOTTOMLINE1 and (subsample2*1148) > BOTTOMLINE2:
#         break

# if len(s_target)>0:
#     print('Ranking in',str(BOTTOMLINE1),'background BCRs and/or',str(BOTTOMLINE2*1148),'background antigens...')
#     # antigen_score_dict,bcr_score_dict = generate_score_dict(background,antigen_ids,antigen_dict,CDR3h_dict,Vh_dict,NPY_DIR,model_mix,model_mix2,
#                                  # antigen_score_dict=antigen_score_dict,bcr_score_dict=bcr_score_dict,subsample1=subsample1,subsample2=subsample2,seed =SEED)
#     # res = calculate_rank(s_target,antigen_score_dict,bcr_score_dict,f_antigens,len_dict,model_mix,model_mix2)
#     output = pd.concat([output,f_res],axis =0)
#     percentage_completed = output.shape[0] / target.shape[0] * 100
#     print(f'Completed {percentage_completed:.2f}% of entries...')

# if ANTIGEN_ONLY and not BCR_ONLY:
#     with open(OUT_DIR+'/background_bcr_score_dict_'+str(BOTTOMLINE1)+'BCRs.pkl', 'wb') as outfile:
#         pickle.dump(bcr_score_dict, outfile)
# elif BCR_ONLY and not ANTIGEN_ONLY:
#     with open(OUT_DIR+'/background_antigen_score_dict_'+str(BOTTOMLINE2)+'Antigens.pkl', 'wb') as outfile:
#         pickle.dump(antigen_score_dict, outfile)
# else:
#     with open(OUT_DIR+'/background_antigen_score_dict_'+str(BOTTOMLINE1)+'BCRs_'+str(BOTTOMLINE2)+'Antigens.pkl', 'wb') as outfile:
#         pickle.dump(antigen_score_dict, outfile)
#     with open(OUT_DIR+'/background_bcr_score_dict_'+str(BOTTOMLINE1)+'BCRs_'+str(BOTTOMLINE2)+'Antigens.pkl', 'wb') as outfile:
#         pickle.dump(bcr_score_dict, outfile)


# In[ ]:




# bcr_score_dict


# In[ ]:





# In[ ]:




#     merged = target.merge(output, on=['record_id','BCR_id'], how='left')
#     merged.to_csv(OUT_DIR+'/merged_results.csv',index=False)


# In[ ]:
