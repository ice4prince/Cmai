#!/usr/bin/env python
# coding: utf-8

# In[75]:


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
module_path = os.path.abspath(os.path.join('..'))
if module_path not in sys.path:
    sys.path.append(module_path)
from datetime import date
from datetime import datetime
import torch.optim as optim
from torch.optim.lr_scheduler import ReduceLROnPlateau
import argparse
print("Today's date:",date.today())
print(str(datetime.now()))


# In[61]:

parser = argparse.ArgumentParser(description='Parameters for pair model.')

# Add a optional argument
parser.add_argument('--batch_size', type=int, help='the batch size, default is 100',default = 100)
#parser.add_argument('--number_of_batch', type=int, help='the number of batch in each epoch, default is 300',default = 300)
parser.add_argument('--epoch', type=int, help='the maximum of epoch, default is 50',default = 100)
parser.add_argument('--tag', type=int, help='the tag to specify the running',default = 1)
parser.add_argument('--initial_learning_rate', type=float, help='the starting leanring rate, default is 0.05', default=0.01)
parser.add_argument('--delta', type=float, help='the minimum difference between better and worse, default is 0.01', default=0.01)
parser.add_argument('--Lambda', type=float, help='the parameter adjustable for loss to avoid INF and NAN', default=0)
parser.add_argument('--weight', type=float, help='the parameter adjustable to shrink embedding of antigens to input', default=0.03)
parser.add_argument('--cut_off', type=int, help='the maximum of length of sequencing, over which will be deleted. default is 1800',default = 1800)
parser.add_argument('--group_size', type = int, help = 'the leap to change the antigen unless there are not enough entries.',default = 100)
parser.add_argument('--Limit', type = float, help = 'the size limit for the antigens in gpu.',default = 3)
parser.add_argument('--verbose', action='store_true', help='Enable verbose output')

# Add a switch (flag) argument
#parser.add_argument('--verbose', action='store_true', help='Print verbose output')
# Parse the command-line arguments
args = parser.parse_args()

BATCH_SIZE = args.batch_size
#NUMBER_BATCH = args.number_of_batch
LR = args.initial_learning_rate
EPOCH = args.epoch
CUTOFF = args.cut_off
TAG = args.tag
T = args.delta
LAMBDA = args.Lambda
w = args.weight
NPY_DIR ='/project/DPDS/Wang_lab/shared/BCR_antigen/data/rfscript/Antigen_embed/CLEANED'
LIMIT = args.Limit
#NPY_DIR ='/project/DPDS/Wang_lab/shared/BCR_antigen/data/rfscript/Antigen_embed/CLEANED'
VERBOSE = args.verbose
GROUP_SIZE = args.group_size

# BATCH_SIZE = 1
# NUMBER_BATCH = 20
# LR = 0.01
# EPOCH = 50
# CUTOFF = 1800
# TAG = 1
# T = 0.01
# LAMBDA = 0
# w = 1
# LIMIT = 3
# NPY_DIR ='/project/DPDS/Wang_lab/shared/BCR_antigen/data/rfscript/Antigen_embed/CLEANED'
# VERBOSE = True
# GROUP_SIZE = 100


# In[3]:


device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")
print(device)


# In[4]:


#os.getcwd()
os.chdir('/project/DPDS/Wang_lab/shared/BCR_antigen/code/DeLAnO/')
#mp.set_sharing_strategy('file_system')
#mp.set_start_method('spawn')
from wrapV.Vwrap import embedV
from wrapCDR3.CDR3wrap import embedCDR3
torch.set_printoptions(precision=10)


# In[5]:


def get_now():
    now = str(datetime.now())
    tosec = now.split('.')[0]
    out_now = tosec.split(' ')[0]+'_'+tosec.split(' ')[1].replace(':',"-")
    return out_now

def filter_big_antigens(dataset,cutoff):
    dataset['aalens']=list(map(len,dataset['Antigen'].tolist()))
    data_filtered = dataset.loc[dataset['aalens']< cutoff]
    print('After removing antigens larger than '+str(cutoff)+', '+str(100*data_filtered.shape[0]/dataset.shape[0])+'% antigens remained.')
    return data_filtered

def preprocess(df,cutoff):
    for col in df[['BetterBCR_Vh','BetterBCR_CDR3h','WorseBCR_Vh','WorseBCR_CDR3h','Antigen']].columns:
        df[col] = df[col].apply(lambda x: x.replace(' ', ''))
    df = filter_big_antigens(df,cutoff)
    df = df.assign(antigen_index = df['Project'] + '/' + df['id'].astype(str))
    df = df.sort_values('antigen_index',)
    df = df.assign(record_id = ['record_' + str(s) for s in range(df.shape[0])])
    return df


def has_space(string):
    return " " in string
def check_space(df):
    df = df[['BetterBCR_Vh','BetterBCR_CDR3h','WorseBCR_Vh','WorseBCR_CDR3h','Antigen']]
    for col in df.columns:
        if df[col].str.contains(" ").any():
            print(Fore.RED +str(col)+' contains space!')
        else:
            print(Fore.GREEN +str(col)+' PASS!')
    print(Style.RESET_ALL)


def aalen(aa_embedding):
    return(aa_embedding.shape[1])

def build_BCR_dict(dataset,colname,precise = False):
    cols = dataset.filter(like = colname)
    uniq_keys = pd.unique(cols.values.ravel()).tolist()
    if colname == 'CDR3h':
        uniq_embedding,*_ = embedCDR3(uniq_keys,precise = precise)
    elif colname == 'Vh':
        uniq_embedding,*_ = embedV(uniq_keys,precise = precise)
    i = 0
    mydict = {}
    for key in uniq_keys:
        mydict[key] = uniq_embedding[i]
        i += 1
    return(mydict)


# In[7]:


# hard_file_train = pd.read_csv('/project/DPDS/Wang_lab/shared/BCR_antigen/data/hard_negative_train_downsampled.csv')
# hard_file_val = pd.read_csv('/project/DPDS/Wang_lab/shared/BCR_antigen/data/hard_negative_validation_downsampled.csv')


# In[8]:


# hard_neg_train = preprocess(hard_file_train,CUTOFF)
# hard_neg_val = preprocess(hard_file_val, CUTOFF)
# print('Train Data: ')
# check_space(hard_neg_train)
# print('Validation Data: ')
# check_space(hard_neg_val)


# In[11]:


# hard_neg_val.columns[[2,3,4,5,15,16,17]]


# In[12]:


# hard_neg_val.loc[:5,['BetterBCR_Vh', 'BetterBCR_CDR3h', 'WorseBCR_Vh', 'WorseBCR_CDR3h',
#        'aalens', 'antigen_index', 'record_id']]


# In[27]:


# hard_neg_val['antigen_index'].value_counts()


# In[102]:


# group_id_train = generate_group_index(hard_neg_train,GROUP_SIZE)
# group_id_val = generate_group_index(hard_neg_val,GROUP_SIZE)


# In[89]:


def get_my_data_list(dataset):
    selected_cols = ['BetterBCR_Vh', 'BetterBCR_CDR3h', 'WorseBCR_Vh', 'WorseBCR_CDR3h',
       'aalens', 'antigen_index', 'record_id']
    ds_to_dict = dataset[selected_cols]
#ds_to_dict = dataset[selected_cols].set_index('record_id')
    my_data_list = ds_to_dict.to_dict(orient='records')
    return my_data_list

def generate_group_index(df,group_size):
    groups = df.groupby('antigen_index').cumcount() // group_size + 1
    # add the group IDs to the original dataframe as a new column 'group'
    df['group'] = groups
    # print the resulting dataframe
    df['group_id'] = df.apply(lambda row: row['antigen_index'] + '_' + str(row['group']), axis=1)
    group_id = df['group_id'].tolist()
    return group_id


# In[33]:


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


# In[24]:


# def extract_antigen(antigen_dict,antigen_index,npy_dir = NPY_DIR,verbose = False):
#     if antigen_index in antigen_dict:
#         embedding = antigen_dict[antigen_index]
#     else:

#         try:
#             embedding = np.load(str(npy_dir+'/'+antigen_index+'.pair.npy'))/w
# #            print(npy.shape)
#             antigen_dict[antigen_index] = embedding
#         except ValueError:
#             print('The embedding of antigen %s cannot be found!' % embeded_Antigen_Files[i])
#         if verbose:
#             print(Fore.RED + 'New antigen added to dictionary:'+ antigen_index+Style.RESET_ALL)
#     return embedding


# In[28]:


# antigen_dict = {}
# antigen_in = {}


# In[29]:


# Shan_Beta = extract_antigen(antigen_dict,'Shan/Beta')
# Chappert_B5 = extract_antigen(antigen_dict,'Chappert/B5')


# In[40]:


# aalens_df = hard_neg_val.loc[:,['antigen_index','aalens']].groupby('antigen_index')['aalens'].mean()
# len_dict = aalens_df.to_dict()
# len_dict


# In[ ]:


#get_antigen_in


# In[63]:


# def check_size(antigen_dict,len_dict,new_antigen,limit_size,datatype = 'float32'):
#     dict_size = 0
#     for key in antigen_dict:
# #        print(key,":")
#         aalen = len_dict[key]
# #        print(aalen)
#         dict_size += predict_size(aalen,datatype = datatype)
# #    print('size in dict:',dict_size)
#     check_size = dict_size + predict_size(len_dict[new_antigen])
#     check_limit = limit_size*(1024**3)
# #    print('size all to check:',check_size)
#     return check_size < check_limit


# In[67]:


# def rotate_dict_in(antigen_dict,len_dict,antigen_index,limit_size=LIMIT):
#     while check_size(antigen_dict,len_dict,new_aantigen_index,limit_size = limit_size) and len(antigen_dict)>0:
#         key = random.choice(list(antigen_dict.keys()))
#         print(antigen_dict.keys())
#         del antigen_dict[key]
#     print('antigens still in:',antigen_dict.keys())


# In[69]:


# def get_antigen_in(antigen_in,antigen_dict,len_dict,antigen_index):
#     if not antigen_index in antigen_in:
#         if not antigen_index in antigen_dict:
#             antigen_to_in = extract_antigen(antigen_dict,antigen_index)
#         else:
#             antigen_to_in = antigen_dict[antigen_index]
#         ##check and import
#         rotate_dict_in(antigen_in,len_dict,antigen_index,limit_size = LIMIT)
#         antigen_in[antigen_index] = torch.from_numpy(antigen_to_in).to(device)
#     return(antigen_in[antigen_index])


# In[114]:


#DATATYPE = 'float16'


# In[122]:


#import torch
#from torch.utils.data import Dataset
#import numpy as np
class SampleDataset(Dataset):
    def __init__(self, antigen_fpath_dict, cdr3_dict, v_dict, your_data_list):
        self.antigen_fpath_dict = antigen_fpath_dict
#        self.antigen_fpath_dict = self.__read_files()
        self.your_data_list = your_data_list
        self.antigen_dict = {}
        self.cdr3_dict = cdr3_dict
        self.v_dict = v_dict
        self.antigen_in = {}
        self.lens_dict = {}

    def __getitem__(self, idx):
        your_dict = self.your_data_list[idx]
        antigen_key = your_dict['antigen_index']
        aalens_key = your_dict['aalens']
        betterCDR_key = your_dict['BetterBCR_CDR3h']
        worseCDR_key = your_dict['WorseBCR_CDR3h']
        betterV_key = your_dict['BetterBCR_Vh']
        worseV_key = your_dict['WorseBCR_Vh']
        index_key = your_dict['record_id']
        self.lens_dict[antigen_key] = aalens_key

#        antigen_feat = self.__extract_antigen(antigen_key)
        ##check whether antigen_index in antigen_in; if not, import it to

        better_feat = self.__embedding_BCR(betterCDR_key,betterV_key,precise = True)
        worse_feat = self.__embedding_BCR(worseCDR_key,worseV_key,precise = True)
#        better_pair = self.__comb_embed(antigen_feat,better_feat)
#        worse_pair = self.__comb_embed(antigen_feat,worse_feat)
        self.__get_antigen_in(antigen_key)
        better_pair = self.__comb_embed_gpu(antigen_key,better_feat)
        worse_pair = self.__comb_embed_gpu(antigen_key,worse_feat)
        better_out = better_pair.squeeze(dim=0)
        worse_out = worse_pair.squeeze(dim=0)
#        out_dict = {}
#        out_dict[index_key] = (better_out, worse_out)
#        return better_out.to(device), worse_out.to(device)#, aalen_key #IF RUN COMBINING_EMBED in cpu
    #IF RUN COMBINING_EMBED in GPU:
        return better_out, worse_out, index_key, antigen_key
    def __comb_embed_gpu(self,antigen_index,BCR_feat):
        len = self.lens_dict[antigen_index]
        single_antigen_g = self.antigen_in[antigen_index]
#        single_antigen_g = torch.from_numpy(single_antigen).to(device)
        single_BCR_g = torch.from_numpy(BCR_feat).to(device)
#        single_BCR_g = torch.from_numpy(BCR_feat).half().to(device)
        BCR_t = torch.tile(single_BCR_g,(1,len,len,1))
        pair_feat_g = torch.cat((single_antigen_g,BCR_t),dim=3)
        del single_antigen_g,single_BCR_g,BCR_t
        return pair_feat_g
#     def __comb_embed(self,single_antigen,BCR_feat):
#         '''process your data'''
#     #    singe_antigen = antigen_dict[antigen_name]
#         lengthen = aalen(single_antigen)
# #         if verbose ==True:
# #             if lengthen > 500:
# #                 print(Fore.RED + 'length of antigen over 500: '+str(lengthen))
# #                 print(Style.RESET_ALL)
#         BCR_tile = np.tile(BCR_feat,(1,lengthen,lengthen,1))
#         pair_np = np.concatenate((single_antigen,BCR_tile),axis=3)
#         pair_feat = torch.from_numpy(pair_np)
#         return(pair_feat)
    def __get_antigen_in(self,antigen_index):
        print('Next antigen:',antigen_index)
        if not antigen_index in self.antigen_in:
            if not antigen_index in self.antigen_dict:
                antigen_to_in = self.__extract_antigen(antigen_index)
            else:
                antigen_to_in = self.antigen_dict[antigen_index]
            ##check and import
            self.__rotate_dict_in(antigen_index,limit_size = LIMIT)
            #single_antigen_in =
            self.antigen_in[antigen_index] = torch.from_numpy(antigen_to_in).to(device)
#            self.antigen_in[antigen_index] = torch.from_numpy(antigen_to_in).half().to(device)

#        return single_antigen_in
    def __rotate_dict_in(self,antigen_index,limit_size=LIMIT):
        while self.__check_size(antigen_index,limit_size = LIMIT) and len(self.antigen_in)>0:
#        while self.__check_size(antigen_index,limit_size = LIMIT,datatype = 'float16') and len(self.antigen_in)>0:

            key = random.choice(list(self.antigen_in.keys()))
            print(self.antigen_in.keys())
            del self.antigen_in[key]
        print('antigens still in:',self.antigen_in.keys())

    def __check_size(self,antigen_index,limit_size= LIMIT,datatype = 'float32'):
        dict_size = 0
        for key in self.antigen_in:
    #        print(key,":")
            aalen = self.lens_dict[key]
    #        print(aalen)
            dict_size += predict_size(aalen,datatype = datatype)
    #    print('size in dict:',dict_size)
        check = dict_size + predict_size(self.lens_dict[antigen_index])
        check_limit = limit_size*(1024**3)
    #    print('size all to check:',check_size)
        return check < check_limit


    def __extract_antigen(self,antigen_name):
        if antigen_name in self.antigen_dict:
            single_antigen = self.antigen_dict[antigen_name]
        else:
            try:
                antigen_import = np.load(str(self.antigen_fpath_dict+'/'+antigen_name+'.pair.npy'))/w
                single_antigen = antigen_import ###ON CPU
#                single_antigen = torch.from_numpy(antigen_import).to(device) ###ON GPU
#            print(npy.shape)
                self.antigen_dict[antigen_name] = single_antigen
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
            df1 = pd.DataFrame([cdr3_seq])
            cdr3_feat,*_ = embedCDR3(df1[0],precise = precise)
            self.cdr3_dict[cdr3_seq]=cdr3_feat
        else:
            cdr3_feat = self.cdr3_dict[cdr3_seq]
        if v_seq not in self.v_dict:
            df2 = pd.DataFrame([v_seq])
            v_feat,*_ = embedV(df2[0],precise = precise)
            self.v_dict[v_seq]=v_feat
        else:
            v_feat = self.v_dict[v_seq]
        bcr_feat = np.concatenate((cdr3_feat,v_feat))
        return bcr_feat

    def __len__(self):
        return len(self.your_data_list)


# In[78]:


class GroupShuffleSampler(Sampler):
    def __init__(self, data_source, group_ids, shuffle=True):
        self.data_source = data_source
        self.group_ids = group_ids
        self.group_indices = self._group_indices()
        self.shuffle = shuffle

    def _group_indices(self):
        """
        Returns a dictionary where the keys are the unique group_ids and
        the values are the indices corresponding to each group_id in the data_source.
        """
        group_indices = {}
        for idx, group_id in enumerate(self.group_ids):
            if group_id not in group_indices:
                group_indices[group_id] = []
            group_indices[group_id].append(idx)
        return group_indices

    def __iter__(self):
        """
        Generates an iterator that yields indices of the data_source
        in shuffled order within groups and shuffled order of groups.
        """
        group_ids = list(self.group_indices.keys())
        if self.shuffle:
            np.random.shuffle(group_ids)
        indices = []
        for group_id in group_ids:
            if self.shuffle:
                np.random.shuffle(self.group_indices[group_id])
            indices += self.group_indices[group_id]
        return iter(indices)

    def __len__(self):
        return len(self.data_source)


# In[79]:


class SelfAttentionPooling(nn.Module):
    """
    Implementation of SelfAttentionPooling
    Original Paper: Self-Attention Encoding and Pooling for Speaker Recognition
    https://arxiv.org/pdf/2008.01077v1.pdf
    """
    def __init__(self, input_dim):
        super(SelfAttentionPooling, self).__init__()
        hidden_dim=10
        self.W1 = nn.Linear(input_dim, hidden_dim)
        self.relu = nn.ReLU()
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
        utter_rep = torch.sum(batch_rep * att_w, dim=1)

        return utter_rep

class mix_model(nn.Module):
    def __init__(self):
        super(mix_model,self).__init__()
        self.model1 = nn.Sequential(
            nn.Linear(318,50),#.to(torch.float64),
            # in (1,len,len,318)
            # out (1,len,len.50)
            nn.ReLU(),
            nn.Linear(50,20),#.to(torch.float64),
            nn.ReLU(),
            nn.Linear(20,20),#.to(torch.float64),
            # out (1,len,len,20)
            nn.ReLU()
        )
        self.model2 = SelfAttentionPooling(input_dim=20)
        self.model2_1 = SelfAttentionPooling(input_dim=20)
        # input_dim = hidden size (number of channels)
#        self.model2 = nn.MultiheadAttention(embed_dim=20,num_heads=N_HEADS,batch_first=True)
        self.model3 = nn.Sequential(
            nn.Linear(20,10),
            nn.ReLU(),
            nn.Linear(10,1)
        )

#        self.sigmoid = nn.Sigmoid()
    def forward(self,input): ###because in getitem, return is .cuda(), now input is on gpu
        x = torch.empty(0)
        x = x.to(device)
        for i in range(len(input)):
            k = input[i]
#            k = comb_embedding(antigen_index[i],BCR_dict[i],antigen_dict_in)
            #.to(device)
#            k = F.normalize(k)
#            print(k)
#            print('k shape:',k.shape,'k device:',k.device,'k datatype:',k.dtype)
            k = self.model1(k)
            ##in (1,len,len,318)
            #out (1,len,len,20)
            #print('After model1:',k.shape)
#        x = torch.permute(x,(0,3,1,2))
            k = self.model2(k.squeeze()).unsqueeze(0)
            #k = torch.mean(k,dim=1)
            #out: (1,len,20)
            #print('After model2:',k.shape)
            k = self.model2_1(k)
            #out: (1,20)
            #print('After attention:',x)
#            k = torch.mean(attn_output,dim=1)
            #out: (1,20)
#            print('x device',x.device,'k device',k.device)
            x = torch.cat((x, k), dim=0)
        #print('After mean:',x)
        #out (batch_size, 20)
        x = F.normalize(x)
        out = self.model3(x)
#        x = x/T
#        print(x.shape)
#        x = self.sigmoid(x)
#        out = torch.exp(x) ###I wanna try this because the output of model is too close to each other.
        return(out)

def loss_function(out_b,out_w):
    loss = torch.sum(torch.relu(out_b-out_w+ T)+LAMBDA*(out_b*out_b + out_w*out_w))
    return loss/out_b.shape[0]


# In[81]:


def train_epoch(dataloader,model,loss_fn,optimizer,verbose = False):
    train_loss = 0.0
    train_success = 0
    i = 0
    model.train()
    for batch in dataloader:
        idx_ls = []
#        print(' Batch:',i)
        if verbose:
            print(' Batch:',i)
            before_enter_model = time.time()
        better, worse, idx, antigen_name = batch
        idx_ls.append(idx[0])
        out_b = model(better)
        out_w = model(worse)
        loss = loss_fn(out_b,out_w)
        if BATCH_SIZE == 1:
            success = int(torch.gt(out_w,out_b))
        else:
            success = torch.sum(torch.gt(out_w,out_b)).item()/out_b.shape[0]
        train_success += success
        if math.isinf(loss) or math.isnan(loss):
            prob_bw = [out_b,out_w]
            print('ERROR: The loss is INF or NaN! '+str(prob_bw))
            break
        optimizer.zero_grad()
        loss.backward()
        optimizer.step()
        train_loss += loss.item()
        if verbose:
            print('record_id',idx[0],':')
            print(Fore.BLUE + 'out_b: '+str(out_b))
            print(Fore.BLUE + 'out_w: '+str(out_w)+Style.RESET_ALL)
            modeling_compelte = time.time()
            print('modeling time:'+str(modeling_compelte - before_enter_model))
            print('loss of the batch:',loss.item())
            print('Success or not of the current batch:',success)
            if i > 0:
                print('number of success / number of records included:',train_success/i)
        i += 1
        # if i > 500:
        #     break
    return train_loss/i,train_success/i,idx_ls


def val_epoch(dataloader,model,loss_fn,verbose = False):
    val_loss = 0.0
    val_success = 0
    i = 0
    model.eval()
    for batch in dataloader:
        idx_ls = []
        if verbose:
            print(' Batch:',i)
            before_enter_model = time.time()
        better, worse, idx, antigen_name = batch
        idx_ls.append(idx[0])
        out_b = model(better)
        out_w = model(worse)
        loss = loss_fn(out_b,out_w)
        val_loss += loss.item()
        if BATCH_SIZE == 1:
            success = int(torch.gt(out_w,out_b))
        else:
            success = torch.sum(torch.gt(out_w,out_b)).item()/out_b.shape[0]
        val_success += success
#        idx_ls.append[idx]
        if verbose:
            print('record_id',idx[0],':')
            print(Fore.BLUE + 'out_b: '+str(out_b))
            print(Fore.BLUE + 'out_w: '+str(out_w)+Style.RESET_ALL)
            modeling_compelte = time.time()
            print('modeling time:'+str(modeling_compelte - before_enter_model))
            print('loss of the batch:',loss.item())
            print('Success or not of the current batch:',success)
            if i > 0:
                print('number of success / number of records included:',val_success/i)
        i += 1
        # if i > 100:
        #     break
    return val_loss/i,val_success/i,idx_ls


# In[82]:


parent_dir = '/project/DPDS/Wang_lab/shared/BCR_antigen/data/output'
out_dir = os.path.join(parent_dir,get_now()+'_tag'+str(TAG))
try:
    os.mkdir(out_dir)
except FileExistsError:
    # directory already exists
    pass
print(out_dir,'is created!')
model_dir = out_dir +'/model_comb'
try:
    os.mkdir(model_dir)
except FileExistsError:
    pass
print(model_dir,'is created!')


# In[83]:


hard_file_train = pd.read_csv('/project/DPDS/Wang_lab/shared/BCR_antigen/data/hard_negative_train_downsampled.csv')
hard_file_val = pd.read_csv('/project/DPDS/Wang_lab/shared/BCR_antigen/data/hard_negative_validation_downsampled.csv')


# In[84]:


hard_neg_train = preprocess(hard_file_train,CUTOFF)
hard_neg_val = preprocess(hard_file_val, CUTOFF)
print('Train Data: ')
check_space(hard_neg_train)
print('Validation Data: ')
check_space(hard_neg_val)


# In[85]:


group_id_train = generate_group_index(hard_neg_train,GROUP_SIZE)
group_id_val = generate_group_index(hard_neg_val,GROUP_SIZE)


# In[87]:


hard_neg_df = pd.concat([hard_neg_train,hard_neg_val],axis = 0)
Vh_dict = build_BCR_dict(hard_neg_df,'Vh',precise = True)
CDR3h_dict = build_BCR_dict(hard_neg_df,'CDR3h',precise = True)


# In[90]:


my_train_list = get_my_data_list(hard_neg_train)
my_val_list = get_my_data_list(hard_neg_val)


# In[123]:


dataset_train = SampleDataset(NPY_DIR,CDR3h_dict,Vh_dict,my_train_list)
dataset_val = SampleDataset(NPY_DIR,CDR3h_dict,Vh_dict,my_val_list)


# In[124]:


sampler_tr = GroupShuffleSampler(dataset_train,group_id_train,shuffle=True)
dl_train = DataLoader(dataset_train,1,sampler=sampler_tr)


# In[105]:


# i = 0
# for batch in dl_train:
#     print('Batch:',i,'antigen:',batch[3][0],'shape:',batch[0].shape,'record_id',batch[2][0])
#     i += 1
#     if i > 3:
#         break


# In[126]:


sampler_vl = GroupShuffleSampler(dataset_val,group_id_val,shuffle=True)
dl_val = DataLoader(dataset_val,1,sampler=sampler_vl)


# In[111]:


# i = 0
# for batch in dl_val:
#     print('Batch:',i,'antigen:',batch[3][0],'shape:',batch[0].shape,'record_id',batch[2][0])
#     i += 1
#     if i > 3:
#         break


# In[120]:


model_mix = mix_model()#.half()
model_mix
#model_mix.to(device)
checkpoint = torch.load('/project/DPDS/Wang_lab/shared/BCR_antigen/data/output/2023-02-18_12-33-25_tag2/model_comb/Batch50BatchNumber_600Lr0.01Epoch54tag2_easy_neg.pth')
model_mix.load_state_dict(checkpoint)
#model_mix.half()
model_mix.to(device)
#for name, param in model_mix.named_parameters():
#    print(name, param.dtype)
optimizer = torch.optim.Adam(model_mix.parameters(),lr=LR)
scheduler = ReduceLROnPlateau(optimizer,mode = 'min',factor=0.1,patience =10,verbose=False,threshold_mode='rel',threshold=1e-5)


# In[121]:


#def main():
    # Set the start method to 'spawn'
# #mp.set_start_method('spawn')
# i = 0
# for record in dl_train:
#     start = time.time()
# #    print(len(record))
# #    print(record[0].shape,record[1].shape)
#     better, worse, idx, antigen_name = record
#     out_b = model_mix(better)
#     out_w = model_mix(worse)
# #    print('b shape and w shaple:',out_b.shape,out_w.shape)
#     loss_ex = loss_function(out_b,out_w)
#     print('loss:',loss_ex.item())
# #    print('gt:',int(torch.gt(out_w,out_b)))
#     success = torch.sum(torch.gt(out_w,out_b)).item()/out_b.shape[0]
#     print('success:',success)
#     print('the',i,'-th batch time:',time.time()-start)
#     i += 1
#     if i > 2:
#         break


# In[37]:


# temp = torch.randn((1,1300,1300,318))
# #temp.dtype
# tensor_size = temp.numel()*temp.element_size()
# print(f"The tensor takes up {tensor_size / (1024**3):.2f} GB of memory.")


# In[113]:


# len = 100
# shape = [1, len, len, 318]
# element_size = np.dtype('float32').itemsize
# num_elements = math.prod(shape)
# # Create a tensor with the specified shape
# #tensor = torch.randn((1,len,len,318))
# # Predict the size of the tensor
# #num_elements = tensor.numel()
# print(num_elements)
# #element_size = tensor.element_size()
# print(element_size)
# tensor_size = num_elements * element_size
# print(tensor_size/(1024**3))


# In[111]:





# In[35]:


# tensor_shape = (1, 1300, 1300, 318)
# #dtype = np.float32
# np.dtype(dtype).itemsize
# #tensor_size = np.prod(tensor_shape) * np.dtype(dtype).itemsize


# In[62]:


#NUMBER_BATCH = 20
print('batch size:',BATCH_SIZE,'EPOCH:', EPOCH,'Learning rate:',LR)
start_time = time.time()
print('Initialling ...')
#antigen_dict = {}
initial_val_loss,initial_val_accuracy,ini_val_index = val_epoch(dl_val,model_mix,loss_function,verbose = VERBOSE)
print('Initial validation:',time.time()-start_time)
validation_complete = time.time()
initial_train_loss, initial_train_accuracy,ini_train_index = val_epoch(dl_train,model_mix,loss_function,verbose = VERBOSE)
print('Initial training:',time.time()-validation_complete)
trian_index_ini = model_dir+"/Batch"+str(BATCH_SIZE)+"Lr"+str(LR)+"tag"+str(TAG)+"_ini_train_index.csv"
val_index_ini = model_dir+"/Batch"+str(BATCH_SIZE)+"Lr"+str(LR)+"tag"+str(TAG)+"_ini_val_index.csv"
with open(trian_index_ini, 'w', newline='') as file:
    writer = csv.writer(file)
    # Write each item in the list as a separate row in the CSV file
    for item in ini_train_index:
        writer.writerow([item])
with open(val_index_ini, 'w', newline='') as file:
    writer = csv.writer(file)
    # Write each item in the list as a separate row in the CSV file
    for item in ini_val_index:
        writer.writerow([item])


# In[20]:


train_LOSS = [initial_train_loss]
train_ACC = [initial_train_accuracy]
val_LOSS = [initial_val_loss]
val_ACC = [initial_val_accuracy]


# In[ ]:
print('Start Training...')
for epoch in range(EPOCH):
    print('epoch:',str(epoch))
    start_epoch_time = time.time()
    train_loss,train_accuracy,train_idx = train_epoch(dl_train,model_mix,loss_function,optimizer,verbose = VERBOSE)
    val_loss,val_accuracy,val_idx = val_epoch(dl_val,model_mix,loss_function,verbose = VERBOSE)
    print('Epoch '+str(epoch)+' train loss: '+str(train_loss)+' val loss: '+str(val_loss)+' val accuracy: '+str(val_accuracy))
    scheduler.step(val_loss) ##should I try scheduler.step(-val_accuracy)?
    train_LOSS.append(train_loss)
    train_ACC.append(train_accuracy)
    val_LOSS.append(val_loss)
    val_ACC.append(val_accuracy)
    state_dict=model_mix.state_dict()
    torch.save(state_dict,model_dir+"/Batch"+str(BATCH_SIZE)+"Lr"+str(LR)+"Epoch"+str(epoch)+"tag"+str(TAG)+"_hard_neg.pth")
    trian_index_file = model_dir+"/Batch"+str(BATCH_SIZE)+"Lr"+str(LR)+"Epoch"+str(epoch)+"tag"+str(TAG)+"_train_index.csv"
    val_index_file = model_dir+"/Batch"+str(BATCH_SIZE)+"Lr"+str(LR)+"Epoch"+str(epoch)+"tag"+str(TAG)+"_val_index.csv"
    with open(trian_index_file, 'w', newline='') as file:
        writer = csv.writer(file)
        # Write each item in the list as a separate row in the CSV file
        for item in train_idx:
            writer.writerow([item])
    with open(val_index_file, 'w', newline='') as file:
        writer = csv.writer(file)
        # Write each item in the list as a separate row in the CSV file
        for item in val_idx:
            writer.writerow([item])
    if math.isnan(train_loss) or math.isnan(val_loss):
        state_dict = model_mix.state_dict()
        torch.save(state_dict,model_dir+"/Error_model_epoch"+str(epoch)+"/Batch"+str(BATCH_SIZE)+"_Lr"+str(LR)+"_tag"+str(TAG)+".pth")
        break
        print('Epoch:',epoch,'train_loss:',train_LOSS[-1])
    if val_accuracy > 0.977:
        print('Validation Accuracy reached 0.977!')
        break
    print('All_through time: '+str(time.time()-start_epoch_time))

# In[ ]:


lossTable = pd.DataFrame({'Epoch':[i for i in range(len(train_LOSS))],'Train_loss':train_LOSS,'Val_loss':val_LOSS,'Val_accuracy':val_ACC})
print(lossTable.head())


# In[ ]:


fig, axs = plt.subplots(2, 2)
for ax in axs.flat:
    ax.label_outer()
lossTable.plot(x='Epoch',y='Train_loss',kind='line',ax=axs[0,0])
lossTable.plot(x='Epoch',y='Val_loss',kind='line',ax=axs[0,1])
lossTable.plot(x='Epoch',y='Val_accuracy',kind='line',ax=axs[1,0])
plt.savefig(out_dir+'/lossPlot_tag'+str(TAG)+'_batchSize'+str(BATCH_SIZE)+'_Epoch'+str(EPOCH)+'_Lr'+str(LR)+'.png',dpi=200)
lossTable.to_csv(out_dir+'/lossTable_tag'+str(TAG)+'_batchSize'+str(BATCH_SIZE)+'_Epoch'+str(EPOCH)+'_Lr'+str(LR)+'.csv')


# In[8]:


# # read in the text files and store them as dataframes
# df_combined = pd.DataFrame()
# colname_ls =['batch','idx','outb','outw','loss']

# for var in colname_ls:
#     print(var)
#     #df_name = 'df_'+ var
#     file_name ='/project/DPDS/Wang_lab/shared/BCR_antigen/code/initial_train_'+var+'.txt'
#     print(file_name)
#     #globals()[df_name]
#     df = pd.read_csv(file_name, sep='\t', header=None)
#     #print(df.head())
#     df_combined = pd.concat([df_combined,df],axis = 1)
#     #print(df_combined.head())
#     # combine the dataframes into a single dataframe with each file as a column
# #df_combined = pd.concat([df1, df2, df3], axis=1)
# print(df_combined.head())
# # save the combined dataframe to a CSV file
# #df_combined.to_csv('combined_file.csv', index=False)

# df_combined.columns = colname_ls
# df_combined['idx'] = df_combined['idx'].str.split(' ', expand=True)[1]
# df_combined['loss']= df_combined['loss'].str.split(':',expand=True)[1]
# df_combined['outb']= df_combined['outb'].str.split(':',expand=True)[1].str.split(',',expand=True)[0]
# df_combined['outw']= df_combined['outw'].str.split(':',expand=True)[1].str.split(',',expand=True)[0]
# df_combined['outb'] =df_combined['outb'].str.split('\[\[',expand=True)[1].str.split('\]\]',expand=True)[0]
# df_combined['outw'] =df_combined['outw'].str.split('\[\[',expand=True)[1].str.split('\]\]',expand=True)[0]
# df_combined.head()

# df_combined.to_csv('/project/DPDS/Wang_lab/shared/BCR_antigen/code/initial_train_summary.csv', index=False)


# In[9]:


# df_combined['success'] = df_combined['outb']<df_combined['outw']


# In[122]:


# df_combined['success'].value_counts()


# In[48]:


# import random
# import torch
# from torch.utils.data import Sampler

# class GroupSampler(Sampler):
#     def __init__(self, groups, group_size=1):
#         self.groups = groups
#         self.group_size = group_size
#         self.group_indices = {}
#         for i, group in enumerate(self.groups):
#             for j in range(0, len(group), self.group_size):
#                 self.group_indices.setdefault(i, []).extend(group[j:j+self.group_size])

#     def __iter__(self):
#         # shuffle the group indices
#         group_indices = list(self.group_indices.keys())
#         random.shuffle(group_indices)

#         # iterate over the shuffled groups and yield the shuffled indices within each group
#         for i in group_indices:
#             indices = self.group_indices[i]
#             random.shuffle(indices)
#             yield from indices

#     def __len__(self):
#         return sum(len(group) for group in self.groups)

# # create a dataset with groups
# groups = [[0, 1, 2], [3, 4, 5, 6], [7, 8]]
# dataset = torch.utils.data.TensorDataset(torch.randn(9, 3))

# # create a GroupSampler for the dataset
# sampler = GroupSampler(groups, group_size=2)

# # create a DataLoader with the GroupSampler
# dataloader = DataLoader(dataset, batch_size=2, sampler=sampler)

# # iterate over the DataLoader and print the indices of each batch
# for indices in dataloader:
#     print(indices)


# In[52]:


# import torch
# from torchvision.datasets import CIFAR10
# from torchvision import transforms

# # Download the CIFAR10 dataset and apply some basic transforms
# transform = transforms.Compose([
#     transforms.ToTensor(),
#     transforms.Normalize((0.5, 0.5, 0.5), (0.5, 0.5, 0.5))
# ])
# dataset = CIFAR10(root='./data', train=True, download=True, transform=transform)

# # Get the labels for each image in the dataset
# labels = [label for _, label in dataset]

# # Assign a group ID to each label
# group_ids = []
# for label in labels:
#     if label < 3:
#         group_ids.append(0)
#     elif label < 6:
#         group_ids.append(1)
#     else:
#         group_ids.append(2)


# In[69]:


#Create a GroupShuffleSampler and a dataloader
# dataloader = torch.utils.data.DataLoader(df, batch_size=1, sampler=sampler)


# In[ ]:


# from torch.utils.data import DataLoader, SubsetRandomSampler, TensorDataset
# # create a sample dataset with 10 samples and 3 groups
# data = torch.randn(10, 3)
# groups = [0, 0, 0, 1, 1, 2, 2, 2, 2, 2]

# # create a dataset object from the data and groups
# dataset = TensorDataset(data)

# # create a list of indices for each group
# group_indices = [torch.where(torch.tensor(groups) == i)[0] for i in set(groups)]

# # shuffle the order of the groups
# group_order = torch.randperm(len(group_indices))

# # loop over the groups in the shuffled order and create a SubsetRandomSampler for each group
# sampler = []
# for i in group_order:
#     group_sampler = SubsetRandomSampler(group_indices[i])
#     sampler.append(group_sampler)

# # create a dataloader object that generates batches using the samplers
# dataloader = DataLoader(dataset, batch_size=3, sampler=sampler)

# # loop over the batches and print their contents
# for batch in dataloader:
#     print(batch)


# In[60]:


# hard_neg_train.to_csv('/project/DPDS/Wang_lab/shared/BCR_antigen/data/hard_neg_train_with_recordid.csv',index=False)


# In[61]:


# hard_neg_val.to_csv('/project/DPDS/Wang_lab/shared/BCR_antigen/data/hard_neg_val_with_recordid.csv',index=False)


# In[52]:


# class CustomSampler(Sampler):
#     def __init__(self, data_source):
#         self.data_source = data_source
#         self.indices = []

#         # Sort the dataset by type and index
#         sorted_indices = np.lexsort((data_source.targets, np.arange(len(data_source))))

#         # Divide the dataset into bins with binwidth 100 or smaller if a type has less than 100 entries
#         bin_indices = []
#         bin_targets = []
#         bin_size = 0
#         for idx in sorted_indices:
#             target = data_source.targets[idx]
#             if bin_size == 0 or (target == bin_targets[0] and bin_size < 100):
#                 bin_indices.append(idx)
#                 bin_targets.append(target)
#                 bin_size += 1
#             else:
#                 np.random.shuffle(bin_indices)
#                 self.indices.extend(bin_indices)
#                 bin_indices = [idx]
#                 bin_targets = [target]
#                 bin_size = 1

#         # Shuffle the indices within the last bin
#         if len(bin_indices) > 0:
#             np.random.shuffle(bin_indices)
#             self.indices.extend(bin_indices)

#     def __iter__(self):
#         return iter(self.indices)

#     def __len__(self):
#         return len(self.indices)


# In[8]:


# def my_collate(batch):
#     # Get the maximum length of the tensors in the batch
#     max_len = max([x[2] for x in batch])

#     # Pad each tensor with zeros to the maximum length
#     #for x in batch:
#     #    print(x[0].shape,x[1].shape,x[2])
#     t1_pad = [F.pad(x[0], [0,0,0,max_len - x[2],0,max_len - x[2]]) for x in batch]
#     t2_pad = [F.pad(x[1], [0,0,0, max_len - x[2], 0, max_len - x[2]]) for x in batch]

#     # Stack the padded tensors and lengths
#     t1_padded = torch.stack(t1_pad, dim=0)
#     t2_padded = torch.stack(t2_pad, dim=0)
#     lengths = torch.tensor([x[2] for x in batch])

#     return t1_padded, t2_padded, lengths


# In[43]:


# import pandas as pd

# # create a sample dataframe with a column 'A' containing values 'a', 'b', 'c', and 'd'
# df = pd.DataFrame({'A': ['a']*150 + ['b']*200 + ['c']*230 + ['d']*90})

# # group the dataframe by column 'A' and assign a unique group ID to each group
# groups = df.groupby('A').cumcount() // 100 + 1

# # add the group IDs to the original dataframe as a new column 'group'
# df['group'] = groups

# # print the resulting dataframe
# df['A_group'] = df.apply(lambda row: row['A'] + '_' + str(row['group']), axis=1)
# # df['A_group'].value_counts()

# # create a dictionary of indices for each group
# group_indices = {}
# for group in df['A_group'].unique():
#     indices = df[df['A_group']==group].index.tolist()
#     group_indices[group] = indices

# # new_dict = {i: group_indices[key] for i, key in enumerate(group_indices.keys())}
# # print(new_dict.keys())# # create a custom sampler that shuffles the order of groups and the order of samples within each group


# In[46]:


# class GroupShuffleSampler(Sampler):
#     def __init__(self, group_indices):
#         self.group_indices = group_indices
#         self.group_order = list(group_indices.keys())
#         random.shuffle(self.group_order)
#         print(self.group_order)
#     def __iter__(self):
#         for i in self.group_order:
#             group_indices = self.group_indices[i]
#             group_sampler = SubsetRandomSampler(group_indices)
#             for index in group_sampler:
#                 yield index

#     def __len__(self):
#         return len(self.group_indices)
# # create a dataloader object that generates batches using the custom sampler
# sampler = GroupShuffleSampler(group_indices)
# dataloader = DataLoader(df, batch_size=5, sampler=sampler)

# # loop over the batches and print their contents
# for batch in dataloader:
#     print(batch)


# In[ ]:


# def train_epoch(model,device,dataloader,loss_fn,optimizer):
#     train_loss, train_euclidean = 0.0, 0.0
#     model.train()
#     for sequence in dataloader:
#         input_seq = sequence[0].to(device)
#         encoded, decoded = model(input_seq)
# #        print('encoded:',encoded)
# #        print('decoded:',decoded)
#         loss = loss_fn(decoded,input_seq)
# #        print('loss:',loss)
#         optimizer.zero_grad()
#         loss.backward()
#         optimizer.step()
#         train_loss += loss.item() * sequence[0].size()[0]
#         train_euclidean += torch.sum(torch.diagonal(torch.cdist(decoded.view(-1,PADDING*5),input_seq.view(-1,PADDING*5)))).item()
#     return train_loss/len(dataloader.sampler), train_euclidean/len(dataloader.sampler)

# def val_epoch(model,device,dataloader,loss_fn):
#     val_loss, val_euclidean = 0.0, 0.0
#     model.eval()
#     for sequence in dataloader:
#         input_seq = sequence[0].to(device)
#         encoded, decoded = model(input_seq)
#         loss = loss_fn(decoded,input_seq)
#         val_loss+=loss.item() * sequence[0].size()[0]
#         val_euclidean += torch.sum(torch.diagonal(torch.cdist(decoded.view(-1,PADDING*5),input_seq.view(-1,PADDING*5)))).item()
#     return val_loss/len(dataloader.sampler), val_euclidean/len(dataloader.sampler)
