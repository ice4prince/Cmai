#!/usr/bin/env python
# coding: utf-8

# In[69]:




import os
import sys
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
from sklearn.model_selection import train_test_split
print("Today's date:",date.today())
print(str(datetime.now()))


# In[ ]:




parser = argparse.ArgumentParser(description='Parameters for pair model.')

# Add a optional argument
parser.add_argument('--code', type=str, help='the Cmai directory',default = '/project/DPDS/Wang_lab/shared/BCR_antigen/code/Cmai')
parser.add_argument('--input',type = str, help = 'the input folder for the preprocessed input',default = 'data/intermediates')
parser.add_argument('--out',type = str, help = 'the directory for output files',default = '/project/DPDS/Wang_lab/shared/BCR_antigen/code/Cmai/data/example')
parser.add_argument('--species', action='store_true', help='match the species of background BCR to the target BCR. NOTE: the species MUST BE specified and unique in the target BCR input.')
parser.add_argument('--seed', type=int, help='the seed for the first 100 background BCRs. To use the prepared embeded 100 BCRs, keep the seed to default 1',default = 1)
parser.add_argument('--subsample', type=int, help='the initial sample size of background BCRs. The default is 100',default = 100)
parser.add_argument('--bottomline', type=int, help='the maximum size for subsample of background BCRs, which should no more than 1000000. The deafult is 10000',default = 10000)
parser.add_argument('--verbose', action='store_true', help='Enable verbose output, default is False.')

args = parser.parse_args()

CODE_DIR = args.code
INPUT_DIR = args.input
OUT_DIR = args.out
MATCHING_SPECIES = args.species
SEED = args.seed
SUBSAMPLE = args.subsample
BOTTOMLINE = args.bottomline
VERBOSE = args.verbose


# In[171]:




# CODE_DIR = '/project/DPDS/Wang_lab/shared/BCR_antigen/code/Cmai'
# #CUTOFF = 1800
# SEED = 1
# VERBOSE = False
# MATCHING_SPECIES = False
# SUBSAMPLE = 100
# BOTTOMLINE = 10000
# OUT_DIR = '/project/DPDS/Wang_lab/shared/BCR_antigen/code/Cmai/data/example'


# In[172]:


if BOTTOMLINE >1000000:
    print('The bottomline cannot be larger than 1,000,000.')
    exit()


# In[71]:


os.chdir(CODE_DIR)


# In[72]:



BACKGROUND = 'data/background/backgroundBCR.csv.gz'
NPY_DIR = INPUT_DIR+'/NPY' ###need to add a command to move the pair.npy under results/pred/ to the intermediates/
MODEL = 'models/model.pth'
INPUT = INPUT_DIR+'/processed_input.csv'

from wrapV.Vwrap import embedV ##input needs to be list of strings
from wrapCDR3.CDR3wrap import embedCDR3 ##input needs to be list of strings

CHANNEL_ANTIGEN = 600
CLIP = None
LR = 0.005
T = 0.005
LAMBDA = 0
w = 100


# In[73]:


print('system version:',sys.version)
device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")
print('device:',device)
torch.set_printoptions(precision=10)


# In[76]:


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

def build_BCR_dict(dataset,colname,precise = False):
    cols = dataset.filter(like = colname)
    uniq_keys = pd.unique(cols.values.ravel()).tolist()
    if colname == 'CDR3h':
        uniq_embedding,_,input_keys = embedCDR3(uniq_keys,precise = precise)
    elif colname == 'Vh':
        uniq_embedding,_,input_keys = embedV(uniq_keys,precise = precise)
    i = 0
    mydict = {}
    for key in input_keys:
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

def get_antigen_dict(df):
    antigen_pool = df['Antigen_id'].unique()
    antigen_dict = {}
    for antigen in antigen_pool:
        antigen_dict[antigen] = np.load(NPY_DIR+'/'+antigen+'.pair.npy')/w
    return antigen_dict


# In[81]:


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
        antigen_key = your_dict['Antigen_id']
#        print(antigen_key)
        aalens_key = self.lens_dict[antigen_key]
        cdr3_key = your_dict['BCR_CDR3h']
        v_key = your_dict['BCR_Vh']
        bcr_key = your_dict['BCR_id']
        index_key = your_dict['record_id']
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
            antigen_import = np.load(self.antigen_fpath_dict+'/'+antigen_name+'.pair.npy')/w
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
        return len(self.your_data_list)


# In[83]:


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
        if subsample_ratio < 1.0:
            if not subsample_ratio == 1/10000:
                self.seed = None
            return dataframe.sample(frac=subsample_ratio,random_state=self.seed)
        else:
            return dataframe

    def __getitem__(self, idx):
        bcr_dict = self.bcr_pool[idx]
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
        if v_seq not in self.v_dict:
            v_feat,*_ = embedV([v_seq],precise = precise)
            v_feat = v_feat[0]
            self.v_dict[v_seq]=v_feat
        else:
            v_feat = self.v_dict[v_seq]
        bcr_feat = np.concatenate((cdr3_feat,v_feat))
        return bcr_feat

    def pool_antigen(self,antigen_tensor,out_n_channel):
#        lengthen = antigen_input.shape[1]
        pooling_layer = nn.AdaptiveAvgPool2d((out_n_channel,out_n_channel))
        output = pooling_layer(antigen_tensor.permute(0,3,1,2)).permute(0,2,3,1)
        return output

    def __len__(self):
        return len(self.bcr_pool)


# In[86]:




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


# In[90]:




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


# In[111]:


def background_scores(dataloader,model):
    score_ls = []
    model.eval()
    for batch in dataloader:
        score = model(batch)
        score_ls.append(score.item())
    return score_ls


# In[125]:


def check_score(dataloader,model):
    res = pd.DataFrame(columns=['record_id', 'Antigen', 'BCR_id','Score'])
#    print(zip(index_idx,antigen_index))
#    print(b_pair.shape,w_pair.shape)
    model.eval()
    for batch in dataloader:
        pair, index_idx, antigen_key, bcr_key =batch ##Change the InLoader, when binary is False, w_pair = score
        out = model(pair)
    #    if binary:
        score = out.cpu().detach().tolist()
        df = pd.DataFrame({'record_id':index_idx,'Antigen':antigen_key,'BCR_id':bcr_key,'Score':score})
        res = pd.concat([res,df],axis=0,ignore_index=True)
    return(res)


# In[128]:


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


# In[141]:


def generate_score_dict(df,score_dict,antigen_dict,cdr3_dict,v_dict,model,subsample=1/10000,seed =SEED):
    for antigen_id, antigen in antigen_dict.items():
        background_loader = DataLoader(rankDataset(df, antigen, cdr3_dict, v_dict,subsample_ratio=subsample,seed=seed),1)
        score_background = background_scores(background_loader,model)
        score_dict[antigen_id]=score_background
    return(score_dict)


# In[143]:


def calculate_rank(df,score_dict,antigen_dict,len_dict,model):
    check_loader = DataLoader(checkDataset(df, antigen_dict, NPY_DIR,len_dict),1)
    res_check = check_score(check_loader,model)
    res_check['Rank'] = res_check.apply(lambda row: locate_rank(row['Score'], score_dict[row['Antigen']]), axis=1)
    return(res_check)


# In[74]:


target_file = pd.read_csv(INPUT) # required columns 'Vh','CDR3h', optional 'species'
background = pd.read_csv(BACKGROUND,compression='gzip') # with columns 'Vh','CDR3h','file','species'


# In[75]:


if 'BCR_species' in target_file.columns:
    unique_values = target_file['BCR_species'].dropna().unique()
    if len(unique_values) == 1 and unique_values[0] in {'human', 'mouse'}:
        species = unique_values[0]
        print('The species is:',species)
        if MATCHING_SPECIES:
            print('matching the background species...')
            background = background[background['species']==species]
    else:
        print('The target species are not unique or not in the background species.')
else:
    print('No species are specified!')


# In[178]:

target = target_file
# target = preprocess(target_file)
# target.to_csv(INPUT_DIR+'/filtered_input.csv')

# In[79]:


antigen_dict = get_antigen_dict(target)
for key,value in antigen_dict.items():
    print(key,value.shape)


# In[77]:


if SEED == 1:
    with open('data/background/default100_V_dict.pkl','rb') as f:
        Vh_dict = pickle.load(f)
    with open('data/background//default100_CDR3_dict.pkl','rb') as f:
        CDR3h_dict = pickle.load(f)
else:
    back100 = background.sample(frac=1/10000, random_state=SEED)
    Vh_dict = build_BCR_dict(back100,'Vh',precise = True)
    CDR3h_dict = build_BCR_dict(back100,'CDR3h',precise = True)
# with open('data/background/default100_V_dict.pkl','wb') as f:
#     pickle.dump(Vh_dict, f)
# with open('data/background//default100_CDR3_dict.pkl','wb') as f:
#     pickle.dump(CDR3h_dict, f)


# In[80]:


len_dict = target.drop_duplicates(subset='Antigen_id').set_index('Antigen_id')['Antigen_seq'].map(len).to_dict()


# In[92]:


model_mix = mix_model()
checkpoint = torch.load(MODEL)
model_mix.load_state_dict(checkpoint)
model_mix.to(device)
optimizer = torch.optim.Adam(model_mix.parameters(),lr=LR)
scheduler = ReduceLROnPlateau(optimizer,mode = 'min',factor=0.1,patience =10,verbose=False,threshold_mode='rel',threshold=1e-5)


# In[144]:


# score_dict = generate_score_dict(background,score_dict,antigen_dict,CDR3h_dict,Vh_dict,model_mix)
# res = calculate_rank(target,score_dict,antigen_dict,len_dict,model_mix)
# for antigen_id, antigen in antigen_dict.items():
#     background_loader = DataLoader(rankDataset(background, antigen, CDR3h_dict, Vh_dict,subsample_ratio=1/10000,seed=SEED),1)
#     score_background = background_scores(background_loader,model_mix)
#     score_dict[antigen_id]=score_background


# In[126]:


# target_dataset = checkDataset(target, antigen_dict, NPY_DIR,len_dict)
# check_loader = DataLoader(target_dataset,1)
# res_check = check_score(check_loader,model_mix)
# res_check['Rank'] = res_check.apply(lambda row: locate_rank(row['Score'], score_dict[row['Antigen']]), axis=1)


# In[134]:


# df_cutoffs = pd.DataFrame({'backsample': [100, 1000, 10000, 100000],
#                            'cutoffs': [0.2, 0.1, 0.05, 0.01]})
# df_cutoffs.to_csv('data/intermediates/able.txt', sep='\t', index=False)


# In[183]:


df = pd.read_csv('paras/cutoff_table.txt', sep='\t')
cutoffs_dict = df.set_index('backsample')['cutoffs'].to_dict()
print('threshold to enter the next level of ranking:',cutoffs_dict)


# In[164]:


# subratio=1/10000
# filtered_df = res[res['Rank']<cutoffs_dict[subratio*1000000]]
# filtered_antigen_dict = {k: v for k, v in antigen_dict.items() if k in filtered_df['Antigen'].values}
# filtered_target = target[target['record_id'].isin(filtered_df['record_id'])]


# In[184]:


score_dict = {}
subsample = SUBSAMPLE
s_target =  target
f_antigens = antigen_dict
output = pd.DataFrame()
while len(s_target)>0 and subsample<=BOTTOMLINE:
    print('Ranking for ',s_target['record_id'].to_list(),'...')
    score_dict = generate_score_dict(background,score_dict,f_antigens,CDR3h_dict,Vh_dict,model_mix,subsample=subsample/1000000,seed=SEED)
    res = calculate_rank(s_target,score_dict,f_antigens,len_dict,model_mix)
    output = pd.concat([output,res[res['Rank']>=cutoffs_dict[subsample]]],axis =0)
    f_res = res[res['Rank']<cutoffs_dict[subsample]]
    f_antigens = {k: v for k, v in f_antigens.items() if k in f_res['Antigen'].values}
    s_target = s_target[s_target['record_id'].isin(f_res['record_id'])]
    subsample = subsample*10


# In[187]:


output.to_csv(OUT_DIR+'/binding_results.csv')
