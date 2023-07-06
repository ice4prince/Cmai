#!/usr/bin/env python
# coding: utf-8

# In[11]:


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


# In[127]:


# parser = argparse.ArgumentParser(description='Parameters for pair model.')

# # Add a optional argument
# parser.add_argument('--code', type=str, help='the CLAnO directory',default = '/project/DPDS/Wang_lab/shared/BCR_antigen/code/CLAnO')
# parser.add_argument('--bcr',type = str, help = 'the input files for BCRs',default = '/project/DPDS/Wang_lab/shared/BCR_antigen/code/CLAnO/data/example/Jun_bcr.csv')
# parser.add_argument('--antigen',type = str, help = 'the fasta file of antigens to input',default = '/project/DPDS/Wang_lab/shared/BCR_antigen/code/CLAnO/data/example/Jun_bcr.csv')
# # parser.add_argument('--mode', type=str, default='binary', choices=['binary', 'continuous'], help='Your choice: binary or continuous. Default is binary.')
# parser.add_argument('--continuous', action='store_true', help='swtich the mode from binary to continuous, default mode is binary.')
# parser.add_argument('--species', action='store_true', help='match the species of background BCR to the target BCR. NOTE: the species MUST BE specified and unique in the target BCR input.')
# # parser.add_argument('--cut_off', type=int, help='the maximum of length of sequencing, over which will be deleted. default is 1800',default = 1800)
# parser.add_argument('--seed', type=int, help='the seed for the first 100 background BCRs. To use the prepared embeded 100 BCRs, keep the seed to default 1',default = 1)
# parser.add_argument('--Limit', type = float, help = 'the size limit for the antigens in gpu.',default = 5)
# parser.add_argument('--max_num', type = float, help = 'the size limit for the antigens in gpu.',default = 10)
# parser.add_argument('--group_size', type = int, help = 'the leap to change the antigen unless there are not enough entries.',default = 50)
# parser.add_argument('--verbose', action='store_true', help='Enable verbose output, default is False.')

# args = parser.parse_args()

# CLANO_DIR = args.code
# BCR_INPUT = args.bcr
# ANTIGEN_INPUT = args.antigen
# CUTOFF = args.cut_off
# SEED = args.seed
# LIMIT = args.Limit
# Max_num = args.max_num
# GROUP_SIZE = args.group_size
# VERBOSE = args.verbose
# CONT = args.continuous
# MATCHING_SPECIES = args.species



# In[126]:


CLANO_DIR = '/project/DPDS/Wang_lab/shared/BCR_antigen/code/CLAnO'
ANTIGEN_INPUT = "/project/DPDS/Wang_lab/shared/BCR_antigen/code/CLAnO/data/example/cleaned.antigen.Jun.fasta"
BCR_INPUT = "/project/DPDS/Wang_lab/shared/BCR_antigen/code/CLAnO/data/example/Jun_bcr.csv"
#CUTOFF = 1800
SEED = 1
GROUP_SIZE = 50
LIMIT = 5
Max_num = 10
VERBOSE = False
CONT = False
MATCHING_SPECIES = False


# In[114]:


os.chdir(CLANO_DIR)
BACKGROUND = 'data/background/backgroundBCR.csv'
NPY_DIR = 'data/intermediates/NPY' ###need to add a command to move the pair.npy under results/pred/ to the intermediates/
MODEL = 'models/binary_model.pth'

from wrapV.Vwrap import embedV ##input needs to be list of strings
from wrapCDR3.CDR3wrap import embedCDR3 ##input needs to be list of strings

CHANNEL_ANTIGEN = 600
CLIP = None
LR = 0.005
T = 0.005
LAMBDA = 0
w = 100


# In[115]:


print('system version:',sys.version)
device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")
print('device:',device)
torch.set_printoptions(precision=10)


# In[116]:


if CONT:
    MODE = 'continuous'
    print('mode switching ...')
else:
    MODE = 'binary'
print('Entering',MODE,'mode now...')


# In[130]:


target_bcr_file = pd.read_csv(BCR_INPUT) # required columns 'Vh','CDR3h', optional 'species'
background = pd.read_csv(BACKGROUND) # with columns 'Vh','CDR3h','file','species'


# In[118]:


if 'BCR_species' in target_bcr_file.columns:
    unique_values = target_bcr_file['BCR_species'].dropna().unique()
    if len(unique_values) == 1 and unique_values[0] in {'human', 'mouse'}:
        species = unique_values[0]
        print('The species is:',species)
        if MATCHING_SPECIES:
            print('matching the background species...')
            backgroud_bcr_file = backgroud_bcr_file[backgroud_bcr_file['species']==species]
    else:
        print('The target species are not unique or not in the background species.')
else:
    print('No species are specified!')


# In[120]:



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

def check_bad_bcr(seq):
    allowed_letters = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
    uppercase_string = seq.upper()
    return not any(char not in allowed_letters for char in uppercase_string)

def filter_bad_bcr(df):
    substrings = ['Vh', 'CDR3h', 'Antigen']
    some_columns = [col for col in df.columns if any(sub in col for sub in substrings)]
    for col in some_columns:
        df[col] = df[col].apply(lambda x: x.replace(' ', ''))
    mask = df[some_columns].applymap(check_bad_bcr)
    filtered_df = df[mask.all(axis=1)]
    return filtered_df

def preprocess(df):
    df = filter_bad_bcr(df)
#     if 'Antigen' in df.columns:
#         df = filter_big_antigens(df,cutoff)
#         df = df.assign(antigen_index = df['Project'] + '/' + df['id'].astype(str))
#     df = df.sort_values('antigen_index',)
    df = df.assign(record_id = ['record_' + str(s) for s in range(df.shape[0])])
    return df

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

# def predict_size(len,datatype = 'float32'):
#     # Define the shape of the tensor
#     shape = [1, len, len, 318]
#     element_size = np.dtype(datatype).itemsize
#     num_elements = math.prod(shape)
#     tensor_size = num_elements * element_size
# #    size_gb = tensor_size/(1024**3)
# #    print(f"The tensor takes up {size_gb:.2f} GB of memory.")
# #    size_gb = tensor_size/(1024*1024*1024)
# #    print(f"Tensor size: {size_gb:.2f} GB")
#     return tensor_size


# In[108]:


#backgroud = preprocess(backgroud_bcr_file,CUTOFF)


# In[112]:


#backgroud.to_csv(BACKGROUND,index=False)


# In[131]:


# 100 background BCR.
# Sample 1% of the rows


# In[134]:


if SEED == 1:
    with open('data/background/default100_V_dict.pkl','rb') as f:
        Vh_dict = pickle.load(f)
    with open('data/background//default100_CDR3_dict.pkl','rb') as f:
        CDR3h_dict = pickle.load(f)
else:
    back100 = background.sample(frac=1/10000, random_state=SEED)
    Vh_dict = build_BCR_dict(subsample,'Vh',precise = True)
    CDR3h_dict = build_BCR_dict(subsample,'CDR3h',precise = True)
# with open('data/background/default100_V_dict.pkl','wb') as f:
#     pickle.dump(Vh_dict, f)
# with open('data/background//default100_CDR3_dict.pkl','wb') as f:
#     pickle.dump(CDR3h_dict, f)


# In[135]:


antigen_dict = {}
with open(ANTIGEN_INPUT, "r") as handle:
    for record in SeqIO.parse(handle, "fasta"):
        # Access the ID and sequence of each record
        antigen_dict[record.id] = str(record.seq)
print(antigen_dict.keys())


# In[ ]:


##MARK HERE###


# In[209]:


#import torch
#from torch.utils.data import Dataset
#import numpy as np
class rankDataset(Dataset):
    def __init__(self, antigen, dataframe, antigen_dict, cdr3_dict, v_dict, antigen_fpath_dict,subsample_ratio=1/10000,seed=SEED):
        self.seed = seed
        self.antigen_fpath_dict = antigen_fpath_dict
        self.antigen_dict = antigen_dict
        self.cdr3_dict = cdr3_dict
        self.v_dict = v_dict
        self.background = self.subsample_data(dataframe, subsample_ratio=subsample_ratio)
        self.bcr_pool = self.background[['Vh','CDR3h']].to_dict(orient='records')
#        print(len(self.bcr_pool))
        self.antigen_feat = self.extract_antigen(antigen)[0].to(device)
        self.lengthen = len(self.antigen_dict[antigen])
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
        if MODE == 'binary':
            lengthen = CHANNEL_ANTIGEN
        else:
            lengthen = self.lengthen
        single_antigen_g = self.antigen_feat  
        single_BCR_g = torch.from_numpy(bcr_feat).to(device)
        BCR_t = torch.tile(single_BCR_g,(lengthen,lengthen,1))
        pair_feat_g = torch.cat((single_antigen_g,BCR_t),dim=2)
        del single_BCR_g,BCR_t
        torch.cuda.empty_cache()
        return pair_feat_g

    def extract_antigen(self,antigen_name):
        try:
            antigen_import = np.load(str(self.antigen_fpath_dict+'/'+antigen_name+'.pair.npy'))/w
#             if not antigen_import.shape[1] == self.lengthen:
#                 print(Fore.RED + 'antigen ' +str(antigen_name)+' embedding'+str(antigen_import.shape[1])+' is NOT in the correct shape '+str(self.lens_dict[antigen_name])+'!'+ Style.RESET_ALL)
#                 exit()            
            single_antigen = torch.from_numpy(antigen_import).float()
            if MODE == 'binary':
                single_antigen = self.pool_antigen(single_antigen,CHANNEL_ANTIGEN) ###ON CPU
        except ValueError:
            print('The embedding of antigen %s cannot be found!' % antigen_name)        
        return single_antigen

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


# In[208]:


class checkDataset(Dataset):
    def __init__(self, antigen, dataframe, antigen_dict, antigen_fpath_dict):
        self.antigen_fpath_dict = antigen_fpath_dict
        self.antigen_dict = antigen_dict
        self.data = dataframe
        self.bcr_pool = self.data[['ID','Vh','CDR3h']].to_dict(orient='records')
        self.antigen_feat = self.extract_antigen(antigen)[0].to(device)
        self.lengthen = len(self.antigen_dict[antigen])

    def __getitem__(self, idx):
        bcr_dict = self.bcr_pool[idx]
        index_key = bcr_dict['ID']
        v_key = bcr_dict['Vh']
        cdr3_key = bcr_dict['CDR3h']        
        bcr_feat = self.__embedding_BCR(cdr3_key,v_key,precise = True)
        pair_feat = self.__comb_embed_gpu(bcr_feat)
        return pair_feat,index_key

    def __comb_embed_gpu(self,bcr_feat):
        if MODE == 'binary':
            lengthen = CHANNEL_ANTIGEN
        else:
            lengthen = self.lengthen
        single_antigen_g = self.antigen_feat  
        single_BCR_g = torch.from_numpy(bcr_feat).to(device)
        BCR_t = torch.tile(single_BCR_g,(lengthen,lengthen,1))
        pair_feat_g = torch.cat((single_antigen_g,BCR_t),dim=2)
        del single_BCR_g,BCR_t
        torch.cuda.empty_cache()
        return pair_feat_g

    def extract_antigen(self,antigen_name):
        try:
            antigen_import = np.load(str(self.antigen_fpath_dict+'/'+antigen_name+'.pair.npy'))/w
#             if not antigen_import.shape[1] == self.lengthen:
#                 print(Fore.RED + 'antigen ' +str(antigen_name)+' embedding'+str(antigen_import.shape[1])+' is NOT in the correct shape '+str(self.lens_dict[antigen_name])+'!'+ Style.RESET_ALL)
#                 exit()            
            single_antigen = torch.from_numpy(antigen_import).float()
            if MODE == 'binary':
                single_antigen = self.pool_antigen(single_antigen,CHANNEL_ANTIGEN) ###ON CPU
        except ValueError:
            print('The embedding of antigen %s cannot be found!' % antigen_name)        
        return single_antigen

    def __embedding_BCR(self,cdr3_seq,v_seq,precise = True):
#         if cdr3_seq not in self.cdr3_dict:
        cdr3_feat,*_ = embedCDR3([cdr3_seq],precise = precise)
        cdr3_feat = cdr3_feat[0]
#             self.cdr3_dict[cdr3_seq]=cdr3_feat
#         else:
#             cdr3_feat = self.cdr3_dict[cdr3_seq]
#         if v_seq not in self.v_dict:
        v_feat,*_ = embedV([v_seq],precise = precise)
        v_feat = v_feat[0]
#             self.v_dict[v_seq]=v_feat
#         else:
#             v_feat = self.v_dict[v_seq]
        bcr_feat = np.concatenate((cdr3_feat,v_feat))
        return bcr_feat
    
    def pool_antigen(self,antigen_tensor,out_n_channel):
#        lengthen = antigen_input.shape[1]
        pooling_layer = nn.AdaptiveAvgPool2d((out_n_channel,out_n_channel))
        output = pooling_layer(antigen_tensor.permute(0,3,1,2)).permute(0,2,3,1)
        return output

    def __len__(self):
        return len(self.bcr_pool)


# In[175]:


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


# In[196]:


class mix_model(nn.Module):
    def __init__(self,mode='binary'):
        super(mix_model,self).__init__()
        if mode =='binary':
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
        else:
            self.model10 = nn.Sequential(
                nn.Linear(318,40),#.to(torch.float64),
                # in (1,len,len,318)
                # out (1,len,len.50)
                nn.LeakyReLU(0.1),
            )

            self.model11 = nn.Sequential(
                nn.AdaptiveAvgPool2d((CHANNEL_ANTIGEN,CHANNEL_ANTIGEN))
            )

            self.model12 = nn.Sequential(
                nn.Linear(40,30),#.to(torch.float64),
                nn.LeakyReLU(0.1),
                nn.Linear(30,20),#.to(torch.float64),
                # out (1,len,len,20)
                nn.LeakyReLU(0.1)
            )
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
    def forward(self,x,mode='binary'): ###because in getitem, return is .cuda(), now input is on gpu
#         x = torch.empty(0)
#         x = x.to(device)
#        x = x.permute(0,2,1,3)
#         print('after permute',x.shape)
        x0 = torch.empty(0)
        x0 = x0.to(device)
        if mode=='binary':
            x = self.model1(x)
            for i in range(len(x)):
                k = x[i]
                k = self.model2(k).unsqueeze(0)
                k = self.model2_1(k)
                x0 = torch.cat((x0, k), dim=0)
        else:
            x = self.model10(x)
            for i in range(len(x)):
                k = x[i]
                k=k.permute(2,0,1)
                k=self.model11(k)
                k=k.permute(1,2,0)
                k=self.model12(k)
                k = self.model2(k).unsqueeze(0)
                k = self.model2_1(k)
                x0 = torch.cat((x0, k), dim=0)
        x0 = F.normalize(x0)
        x0 = self.model3(x0).squeeze()
#        if binary:
        out  = x0
        return(out)


# In[142]:


def background_scores(dataloader,model):
    score_ls = []
    model.eval()
    for batch in dataloader:
        score = model(batch)
        score_ls.append(score.item())
    return score_ls

def check_scores(dataloader,model):
    score_dict = {}
    model.eval()
    for batch, index in dataloader:
        score_dict[index] = model(batch).item()
    return score_dict 


# In[146]:


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


# In[197]:


model_mix = mix_model(mode=MODE)
checkpoint = torch.load(MODEL)
model_mix.load_state_dict(checkpoint)
model_mix.to(device)
optimizer = torch.optim.Adam(model_mix.parameters(),lr=LR)
scheduler = ReduceLROnPlateau(optimizer,mode = 'min',factor=0.1,patience =10,verbose=False,threshold_mode='rel',threshold=1e-5)


# In[257]:


# MARK HERE: Is it Possible to remove the alpha1 beta1 alpha2 beta2 from checkpoint?
# REMEMBER to write if 1%, expand the subsample...


# In[195]:


# # Remove the parameters you don't want
# del checkpoint['alpha1']
# del checkpoint['alpha2']
# del checkpoint['beta1']
# del checkpoint['beta2']

# # Save the modified checkpoint
# torch.save(checkpoint,MODEL)


# In[213]:


for antigen in antigen_dict.keys():
    print('For antigen:',antigen)
#     background_loader =rankDataset(antigen, background, antigen_dict, CDR3h_dict, Vh_dict,NPY_DIR,subsample_ratio=1/10000,seed=1)
    background_dataset =rankDataset(antigen, background, antigen_dict, CDR3h_dict, Vh_dict,NPY_DIR,subsample_ratio=1/10000,seed=1)
    background_loader = DataLoader(background_dataset, 1)
    back_scores = background_scores(background_loader,model_mix)
    target_dataset =checkDataset(antigen, target_bcr_file, antigen_dict, NPY_DIR)
    target_loader =DataLoader(target_dataset, 1)
    target_score = check_scores(target_loader,model_mix)
    #print('the binding score of the tested bcr is:',target_score)
    rank = {}
    for bcr, score in target_score.items():
        rank[bcr[0]]= locate_rank(score,back_scores)
    print('the binding rank of the tested bcr is:',rank)


# In[ ]:


#mark here to export

