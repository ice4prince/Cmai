#!/usr/bin/env python
# coding: utf-8

# In[ ]:





# In[2]:




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
from tqdm import tqdm
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
from sklearn.model_selection import train_test_split
print("Today's date:",date.today())
print(str(datetime.now()))


# In[ ]:






# In[ ]:






parser = argparse.ArgumentParser(description='Parameters for pair model.')

# Add a optional argument
parser.add_argument('--input', type=str, help='the input file tag',default = '3diff_easy_negative')
parser.add_argument('--batch_size', type=int, help='the batch size, default is 100',default = 1)
parser.add_argument('--seed', type=int, help='the seed for dataloader, default is None',default = None)
parser.add_argument('--epoch', type=int, help='the maximum of epoch, default is 50',default = 50)
parser.add_argument('--tag', type=int, help='the tag to specify the running',default = 1)
parser.add_argument('--initial_learning_rate', type=float, help='the starting leanring rate, default is 0.05', default=0.01)
parser.add_argument('--delta', type=float, help='the minimum difference between better and worse, default is 0.01', default=0.01)
parser.add_argument('--Lambda', type=float, help='the parameter adjustable for loss to avoid INF and NAN', default=0)
parser.add_argument('--weight', type=float, help='the parameter adjustable to shrink embedding of antigens to input', default=100)
parser.add_argument('--cut_off', type=int, help='the maximum of length of sequencing, over which will be deleted. default is 1800',default = 1800)
parser.add_argument('--group_size', type = int, help = 'the leap to change the antigen unless there are not enough entries.',default = 100)
parser.add_argument('--Limit', type = float, help = 'the size limit for the antigens in gpu.',default = 5)
parser.add_argument('--max_num', type = float, help = 'the size limit for the antigens in gpu.',default = 10)
parser.add_argument('--model',type = str, help = 'the model file to import',default = None)
parser.add_argument('--resize',type = int, help = 'the length after resizing',default = 500)
parser.add_argument('--small_sample',type = int, help = 'whether downsample for experiment.',default = None)
parser.add_argument('--inter_sub',type = float, help = 'the ratio of subsampling for training and internal validation.',default = 0.01)
parser.add_argument('--outer_sub',type = float, help = 'the ratio of subsampling for the external validation',default = 1)
# parser.add_argument('--model',type = str, help = 'the model file to import',default = '/project/DPDS/Wang_lab/shared/BCR_antigen/data/output/2023-02-18_12-33-25_tag2/model_comb/Batch50BatchNumber_600Lr0.01Epoch54tag2_easy_neg.pth')
parser.add_argument('--verbose', action='store_true', help='Enable verbose output')

# Add a switch (flag) argument
#parser.add_argument('--verbose', action='store_true', help='Print verbose output')
# Parse the command-line arguments
args = parser.parse_args()

BATCH_SIZE = args.batch_size
INPUT = args.input
#NUMBER_BATCH = args.number_of_batch
LR = args.initial_learning_rate
EPOCH = args.epoch
CUTOFF = args.cut_off
TAG = args.tag
SEED = args.seed
T = args.delta
LAMBDA = args.Lambda
w = args.weight
NPY_DIR ='/project/DPDS/Wang_lab/shared/BCR_antigen/data/rfscript/Antigen_embed/CLEANED'
LIMIT = args.Limit
#NPY_DIR ='/project/DPDS/Wang_lab/shared/BCR_antigen/data/rfscript/Antigen_embed/CLEANED'
VERBOSE = args.verbose
GROUP_SIZE = args.group_size
MODEL = args.model
Max_num = args.max_num
Small_sample = args.small_sample
EX_SUBSAMPLE = args.outer_sub
IN_SUBSAMPLE = args.inter_sub
CHANNEL_ANTIGEN = args.resize

# In[247]:


# INPUT = '3diff_easy_negative'
# BATCH_SIZE = 1
# LR = 0.005
# EPOCH = 20
# CUTOFF = 1800
# TAG = 1
# SEED = 44
# T = 0.05
# LAMBDA = 0
# w = 100
# LIMIT = 5
# NPY_DIR ='/project/DPDS/Wang_lab/shared/BCR_antigen/data/rfscript/Antigen_embed/CLEANED'
# VERBOSE = False
# GROUP_SIZE = 100
# Max_num = 2
# Small_sample = 100000
# EX_SUBSAMPLE = 0.1
# IN_SUBSAMPLE = 0.01
# CHANNEL_ANTIGEN = 500
# MODEL = 'full_validation_2023-04-11/model_comb/selected_3diff_epoch14_0.71accu.pth'


# In[248]:




print('batch size:',BATCH_SIZE,
      '\ninput label:',INPUT,
      '\nEpoch max:',EPOCH,
      '\nLearning rate:',LR,
      '\nCut off:',CUTOFF,
      '\ntag:',TAG,
      '\ndelta:',T,
      '\nLambda:',LAMBDA,
      '\nweight:',w,
      '\nlimit:',LIMIT,
      '\nverbose:',VERBOSE,
      '\ngroup_size:',GROUP_SIZE,
      '\nmodel:',MODEL,
      '\nmax number of antigen:',Max_num,
      '\number of channels:',CHANNEL_ANTIGEN,
     '\nsmall_sample:',Small_sample)


# In[249]:


print(sys.version)

device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")
print(device)


# In[250]:


#os.getcwd()
os.chdir('/project/DPDS/Wang_lab/shared/BCR_antigen/code/DeLAnO/')
#mp.set_sharing_strategy('file_system')
#mp.set_start_method('spawn')
from wrapV.Vwrap import embedV
from wrapCDR3.CDR3wrap import embedCDR3
torch.set_printoptions(precision=10)


# In[251]:




input_ori = pd.read_csv('/project/DPDS/Wang_lab/shared/BCR_antigen/data/'+INPUT+'.csv',index_col = 0)
total_N_row = input_ori.shape[0]

#exVal_ori = pd.read_csv('/project/DPDS/Wang_lab/shared/BCR_antigen/data/results_hard_negative/results_hard_negative.csv')
exVal_ori = pd.read_csv('/project/DPDS/Wang_lab/shared/BCR_antigen/data/results_hard_negative/results_hard_negative_filtered.csv',index_col=0)
total_N_ex = exVal_ori.shape[0]


# In[252]:




# keep_antigen = ['46472/COV_HKU1','46472/COV_MERS','46472/COV_OC43','46472/COV_SARS','46472/COV_SARS2','54042-4/CoV_229E','54042-4/CoV_MERS','54042-4/CoV_NL63','54042-4/CoV_OC43','54042-4/CoV_SARS1','54042-4/HA_NC99','Chappert/B5','Dugan/NP','Dugan/NPSARS2','Dugan/ORF8','Dugan/ORF8SARS2',
# 'Dugan/Spike','Dugan/SPKSARS2','Libra-seq/HA_Anhui','Libra-seq/HA_indo','Libra-seq/HA_NC99','Makowski/HGFR','Makowski/ova','Mason/HER2','Shan/Alpha','Shan/Beta','Shan/Delta','Shan/DeltaPlus','Shan/Epsilon','Shan/Eta','Shan/Gamma','Shan/N439K','Shan/WT614G','Shiakolas/CoV_ACE2',]


# In[253]:




# exVal_ori['antigen_index'] = exVal_ori['Project']+'/'+exVal_ori['id']
# mask= exVal_ori['antigen_index'].isin(bad_antigen)
# exVal_filtered = exVal_ori[mask]


# In[254]:




# exVal_filtered.reset_index().iloc[:,1:-1].to_csv('/project/DPDS/Wang_lab/shared/BCR_antigen/data/results_hard_negative/results_hard_negative_filtered.csv')


# In[255]:


with open('/project/DPDS/Wang_lab/shared/BCR_antigen/data/toInput/easy_V_dict.pkl','rb') as f:
    Vh_dict = pickle.load(f)
with open('/project/DPDS/Wang_lab/shared/BCR_antigen/data/toInput/easy_CDR3_dict.pkl','rb') as f:
    CDR3h_dict = pickle.load(f)
with open('/project/DPDS/Wang_lab/shared/BCR_antigen/data/toInput/hard_V_dict.pkl', 'rb') as f:
    Vh_hard_dict = pickle.load(f)
with open('/project/DPDS/Wang_lab/shared/BCR_antigen/data/toInput/hard_CDR3_dict.pkl', 'rb') as f:
    CDR3h_hard_dict = pickle.load(f)


# In[256]:




# import gc
# del model_mix
# torch.cuda.empty_cache()
# gc.collect()


# In[257]:


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

def check_bad_bcr(seq):
    allowed_letters = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
    uppercase_string = seq.upper()
    return not any(char not in allowed_letters for char in uppercase_string)

def filter_bad_bcr(df):
    mask = df[['BetterBCR_Vh', 'BetterBCR_CDR3h', 'WorseBCR_Vh', 'WorseBCR_CDR3h']].applymap(check_bad_bcr)
    filtered_df = df[mask.all(axis=1)]
    return filtered_df

def preprocess(df,cutoff):
    for col in df[['BetterBCR_Vh','BetterBCR_CDR3h','WorseBCR_Vh','WorseBCR_CDR3h','Antigen']].columns:
        df[col] = df[col].apply(lambda x: x.replace(' ', ''))
    df =  filter_bad_bcr(df)
    df = filter_big_antigens(df,cutoff)
    df = df.assign(antigen_index = df['Project'] + '/' + df['id'].astype(str))
    df = df.sort_values('antigen_index',)
    df = df.assign(record_id = ['record_' + str(s) for s in range(df.shape[0])])
    return df

# def has_space(string):
#     return " " in string
def check_bad_letters(df):
    df = df[['BetterBCR_Vh','BetterBCR_CDR3h','WorseBCR_Vh','WorseBCR_CDR3h','Antigen']]
    for col in df.columns:
        if not all(df[col].apply(check_bad_bcr)):
            print(Fore.RED +str(col)+' contains uncommon aa or symbols!')
        else:
            print(Fore.GREEN +str(col)+' PASS!')
    print(Style.RESET_ALL)


# def aalen(aa_embedding):
#     return(aa_embedding.shape[1])

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


# In[258]:




if Small_sample is not None:
    input_df = input_ori.sample(n = Small_sample)
    exVal = exVal_ori.sample(n = Small_sample)
else:
    input_df = input_ori
    exVal = exVal_ori


# In[259]:




input_df = preprocess(input_df,CUTOFF)
print(INPUT,'Data:')
check_bad_letters(input_df)

exVal = preprocess(exVal,CUTOFF)
print('External Validation Data:')
check_bad_letters(exVal)


# In[260]:


# def block_accuracy(cat_ls,num_ls):
#     category_data = {}
#     for cat, num in zip(cat_ls, num_ls):
#         if cat not in category_data:
#             category_data[cat] = {'sum': 0, 'count': 0}
#         category_data[cat]['sum'] += num
#         category_data[cat]['count'] += 1
#         category_data[cat]['mean'] = category_data[cat]['sum']/category_data[cat]['count']
#     mean_of_values = sum(entry['mean'] for entry in category_data.values()) / len(category_data)
#     return category_data,mean_of_values
# def accu_add_new_entry(category_data,add):
#     if add[0] not in category_data:
#         category_data[add[0]] = {'sum': 0, 'count': 0, 'mean':0}
#     category_data[add[0]]['sum'] += add[1]
#     category_data[add[0]]['count'] += 1
#     category_data[add[0]]['mean'] = category_data[add[0]]['sum']/category_data[add[0]]['count']
#     mean_of_values = sum(entry['mean'] for entry in category_data.values()) / len(category_data)
#     return category_data,mean_of_values


# In[261]:


def batch_result(batch,model,loss_fn):
    b_pair,w_pair,index_idx,antigen_index =batch
#    print(zip(index_idx,antigen_index))
#    print(b_pair.shape,w_pair.shape)
    out_b = model_mix(b_pair)
    out_w = model_mix(w_pair)
    loss = loss_fn(out_b,out_w)
    success = torch.gt(out_w,out_b).tolist()
    outb = out_b.tolist()
    outw = out_w.tolist()
    df = pd.DataFrame({'record_id':index_idx,'antigen':antigen_index,'out_b':outb,'out_w':outw,'success':success})
    return(df,loss)


# In[262]:


# def extract_antigen(self,antigen_name,verbose=False):
#     if antigen_name in self.antigen_dict:
#         single_antigen = self.antigen_dict[antigen_name]
#     else:
#         try:
#             antigen_import = np.load(str(self.antigen_fpath_dict+'/'+antigen_name+'.pair.npy'))/w
#             if not antigen_import.shape[1] == self.lens_dict[antigen_name]:
#                 print(Fore.RED + 'antigen ' +str(antigen_name)+' embedding'+str(antigen_import.shape[1])+' is NOT in the correct shape '+str(self.lens_dict[antigen_name])+'!'+ Style.RESET_ALL)
#                 exit()
#             single_antigen = antigen_import ###ON CPU
# #                single_antigen = torch.from_numpy(antigen_import).to(device) ###ON GPU
# #            print(npy.shape)
#             self.antigen_dict[antigen_name] = single_antigen
#             if verbose:
#                 print(Fore.RED + 'New antigen added to dictionary:'+antigen_name+Fore.RESET)
#         except ValueError:
#             print('The embedding of antigen %s cannot be found!' % antigen_name)
# #             if verbose:
# #                 print(Fore.RED + 'New antigen added to dictionary:',antigen_name)
# #                 print(Fore.RED + 'Number of antigens included:',len(self.antigen_dict))
# #                 print(Style.RESET_ALL)
#     return single_antigen


# In[263]:


#import torch
#from torch.utils.data import Dataset
#import numpy as np
#import torch
#from torch.utils.data import Dataset
#import numpy as np
class SampleDataset(Dataset):
    def __init__(self, dataframe, antigen_fpath_dict, cdr3_dict, v_dict, lens_dict, subsample_ratio = 1.0):
        self.dataframe = self.subsample_data(dataframe, subsample_ratio)
        self.group_ids = self.generate_group_index(group_size=GROUP_SIZE)  # Adjust the group_size as needed
        self.your_data_list = self.get_my_data_list()
        self.antigen_fpath_dict = antigen_fpath_dict
#        self.antigen_fpath_dict = self.__read_files()
#        self.your_data_list = your_data_list
        self.antigen_dict = {}
        self.cdr3_dict = cdr3_dict
        self.v_dict = v_dict
        self.antigen_in = {}
        self.lens_dict = lens_dict

    def subsample_data(self, dataframe, subsample_ratio):
        if subsample_ratio < 1.0:
            return dataframe.sample(frac=subsample_ratio)
        else:
            return dataframe

    def __getitem__(self, idx):
        your_dict = self.your_data_list[idx]
        antigen_key = your_dict['antigen_index']
        aalens_key = your_dict['aalens']
        betterCDR_key = your_dict['BetterBCR_CDR3h']
        worseCDR_key = your_dict['WorseBCR_CDR3h']
        betterV_key = your_dict['BetterBCR_Vh']
        worseV_key = your_dict['WorseBCR_Vh']
        index_key = your_dict['record_id']
#        self.lens_dict[antigen_key] = aalens_key
#        antigen_feat = self.__extract_antigen(antigen_key)
        ##check whether antigen_index in antigen_in; if not, import it to

        better_bcr = self.__embedding_BCR(betterCDR_key,betterV_key,precise = True)
        worse_bcr = self.__embedding_BCR(worseCDR_key,worseV_key,precise = True)
#        print('better shape',better_feat.shape)
#         better_pair = self.__comb_embed(antigen_key,better_feat)
#         worse_pair = self.__comb_embed(antigen_key,worse_feat)
        self.__get_antigen_in(antigen_key)
        antigen_feat = self.antigen_in[antigen_key]
#        print('antigen shape',self.antigen_in[antigen_key].shape)
        better_pair = self.__comb_embed_gpu(antigen_key,better_bcr)
        worse_pair = self.__comb_embed_gpu(antigen_key,worse_bcr)
        return better_pair, worse_pair, index_key, antigen_key#better_feat, worse_feat,

#        better_out = better_pair.squeeze(dim=0)
#        worse_out = worse_pair.squeeze(dim=0)
        # print('better out dtype:',better_out.dtype)
        # print('worse out dtype:',worse_out.dtype)
#        out_dict = {}
#        out_dict[index_key] = (better_out, worse_out)
#        return better_out.to(device), worse_out.to(device)#, aalen_key #IF RUN COMBINING_EMBED in cpu
    #IF RUN COMBINING_EMBED in GPU:

    def get_my_data_list(self,selected_cols = ['BetterBCR_Vh', 'BetterBCR_CDR3h', 'WorseBCR_Vh', 'WorseBCR_CDR3h','Antigen','aalens', 'antigen_index', 'record_id']):
        ds_to_dict = self.dataframe[selected_cols]
    #ds_to_dict = dataset[selected_cols].set_index('record_id')
        my_data_list = ds_to_dict.to_dict(orient='records')
        return my_data_list

    def generate_group_index(self, group_size):
        df = self.dataframe.copy()
        groups = df.groupby('antigen_index').cumcount() // group_size + 1
        df['group'] = groups
        df['group_id'] = df.apply(lambda row: row['antigen_index'] + '_' + str(row['group']), axis=1)
        group_id = df['group_id'].tolist()
        return group_id

    def __comb_embed_gpu(self,antigen_name,BCR_feat):
#        print('antigen to comb:',antigen_name)
        lengthen = CHANNEL_ANTIGEN
        #rint('The current antigen is: ',antigen_name)
#        print('length get from dict_len:',lengthen)
        single_antigen_g = self.antigen_in[antigen_name][0]
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

    def __get_antigen_in(self,antigen_name):
        #print('Next antigen:',antigen_name)
        if not antigen_name in self.antigen_in:
            if not antigen_name in self.antigen_dict:
                antigen_to_in = self.extract_antigen(antigen_name)
            else:
                antigen_to_in = self.antigen_dict[antigen_name]
            ##check and import
            self.__rotate_dict_in(antigen_name,limit_size = LIMIT)
            #single_antigen_in =
            try:

                antigen_tensor = torch.from_numpy(antigen_to_in).to(device)
                self.antigen_in[antigen_name] = self.pool_antigen(antigen_tensor,CHANNEL_ANTIGEN) ###ON CPU

            except RuntimeError as e:
                if "CUDA out of memory" in str(e):
                    print("CUDA out of memory. Clearing cache and trying again...")
                    for key in self.antigen_in:
                        self.antigen_in.pop(key, None)
                    self.antigen_in.clear()
                    torch.cuda.empty_cache()
                    # Try the operation again
                    try:
                        self.antigen_in[antigen_name] = torch.from_numpy(antigen_to_in).to(device)
                    except RuntimeError as e:
                        if "CUDA out of memory" in str(e):
                            print("Still out of memory after clearing cache.")
                        else:
                            raise
                else:
                    raise
#            self.antigen_in[antigen_name] = torch.from_numpy(antigen_to_in).half().to(device)

#        return single_antigen_in
    def __rotate_dict_in(self,antigen_name,limit_size=LIMIT,verbose=False):
#        while self.__check_size(antigen_name,limit_size = LIMIT) and len(self.antigen_in)>0:
        while self.__check_size(antigen_name,limit_size = LIMIT,datatype = 'float32') and len(self.antigen_in)>0 or len(self.antigen_in)> Max_num:

            key = random.choice(list(self.antigen_in.keys()))
#            print(self.antigen_in.keys())
            self.antigen_in.pop(key, None)
            if verbose:
                print(Fore.RED +'the antigen:'+str(key)+'is deleted!'+ Fore.RESET)
            torch.cuda.empty_cache()
#        if len(self.antigen_in) >3:
            if verbose:
                print(Fore.RED +'antigens still in:'+str(self.antigen_in.keys())+ Fore.RESET)

    def pool_antigen(self,antigen_input,out_n_channel):
#        lengthen = antigen_input.shape[1]
        pooling_layer = nn.AdaptiveAvgPool2d((out_n_channel,out_n_channel))
        output = pooling_layer(antigen_input.permute(0,3,1,2)).permute(0,2,3,1)
        return output
    def __check_size(self,antigen_name,limit_size= LIMIT,datatype = 'float32'):
        dict_size = 0
        for key in self.antigen_in:
    #        print(key,":")
            aalen = self.lens_dict[key]
    #        print(aalen)
            dict_size += predict_size(aalen,datatype = datatype)
    #    print('size in dict:',dict_size)
        check = dict_size + predict_size(self.lens_dict[antigen_name])
        check_limit = limit_size*(1024**3)
    #    print('size all to check:',check_size)
        return check > check_limit


    def extract_antigen(self,antigen_name,verbose=False):
        if antigen_name in self.antigen_dict:
            single_antigen = self.antigen_dict[antigen_name]
        else:
            try:
                antigen_import = np.load(str(self.antigen_fpath_dict+'/'+antigen_name+'.pair.npy'))/w
                if not antigen_import.shape[1] == self.lens_dict[antigen_name]:
                    print(Fore.RED + 'antigen ' +str(antigen_name)+' embedding'+str(antigen_import.shape[1])+' is NOT in the correct shape '+str(self.lens_dict[antigen_name])+'!'+ Style.RESET_ALL)
                    exit()
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
        return len(self.your_data_list)


# In[264]:


class GroupShuffleSampler(Sampler):
    def __init__(self, data_source, group_ids, shuffle=True):
        self.data_source = data_source
        self.group_ids = group_ids
        self.group_indices = self._group_indices()
        self.shuffle = shuffle

    def _group_indices(self):
        """b
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


# In[265]:


class InLoader:
    def __init__(self, dataframe, antigen_fpath_dict, cdr3_dict, v_dict, lens_dict,train_ratio=0.8, subsample_ratio=.01, seed=SEED):
        self.dataframe = dataframe
        self.train_ratio = train_ratio
        self.subsample_ratio = subsample_ratio
        self.seed = seed
        self.antigen_fpath_dict = antigen_fpath_dict
        self.cdr3_dict = cdr3_dict
        self.v_dict = v_dict
        self.lens_dict = lens_dict
        self.train_df, self.val_df = self.split_data()

    def split_data(self):
        train_df, val_df = train_test_split(
            self.dataframe, train_size=self.train_ratio, random_state=self.seed)
        return train_df, val_df

    def get_train_dataset(self):
        return SampleDataset(self.train_df, self.antigen_fpath_dict, self.cdr3_dict, self.v_dict, self.lens_dict, subsample_ratio=self.subsample_ratio)

    def get_val_dataset(self):
        return SampleDataset(self.val_df, self.antigen_fpath_dict, self.cdr3_dict, self.v_dict, self.lens_dict, subsample_ratio=self.subsample_ratio)

    def new_epoch(self):
        print("Starting a new epoch...")
        self.train_df, self.val_df = self.split_data()


# In[266]:


# ##MARK HERE###

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

        utter_rep = torch.sum(batch_rep* att_w, dim=1)

        return utter_rep

class mix_model(nn.Module):
    def __init__(self):
        super(mix_model,self).__init__()
        self.model1 = nn.Sequential(
            nn.Linear(318,50),#.to(torch.float64),
            # in (1,len,len,318)
            # out (1,len,len.50)
            nn.SELU(),
            nn.Linear(50,20),#.to(torch.float64),
            nn.SELU(),
            nn.Linear(20,20),#.to(torch.float64),
            # out (1,len,len,20)
            nn.SELU(),
#            nn.Dropout(p=0.2),
        )
        self.model2 = SelfAttentionPooling(input_dim=20,hidden_dim=10)
        self.model2_1 = SelfAttentionPooling(input_dim=20,hidden_dim=10)
        # input_dim = hidden size (number of channels)
#        self.model2 = nn.MultiheadAttention(embed_dim=20,num_heads=N_HEADS,batch_first=True)
        self.model3 = nn.Sequential(
#            nn.Linear(50,20),
#            nn.ReLU(),
            nn.Linear(20,10),
            nn.SELU(),
            nn.Linear(10,1)
        )
#         self.model4 = nn.Sequential(
#             nn.Linear(50,20),
#             nn.ReLU()
#         )

#        self.sigmoid = nn.Sigmoid()
    def forward(self,x): ###because in getitem, return is .cuda(), now input is on gpu
#         x = torch.empty(0)
#         x = x.to(device)
#        x = x.permute(0,2,1,3)
#         print('after permute',x.shape)
        x = self.model1(x)
#         print('after model1',x.shape)
        x0 = torch.empty(0)
        x0 = x0.to(device)
        for i in range(len(x)):
            k = x[i]
            k = self.model2(k).unsqueeze(0)
#             print('after model2',k.shape)
            k = self.model2_1(k)
#             print('after model2_1',k.shape)
            x0 = torch.cat((x0, k), dim=0)
#         print('after loop:',x0.shape)
        x0 = F.normalize(x0)
        out = self.model3(x0)
#         print('after model3',out.shape)


#         for i in range(len(pair)):
#             k = pair[i]
# #            print('start in model:',k.shape,k)
# #            k = comb_embedding(antigen_index[i],BCR_dict[i],antigen_dict_in)
#             #.to(device)
# #            k = F.normalize(k)
# #            print(k)
# #            print('k shape:',k.shape,'k device:',k.device,'k datatype:',k.dtype)
#             k = k.permute(1,0,2)
#             k = self.model1(k)
# #            print('after model1:',k.shape, k)
#             ##in (1,len,len,318)
#             #out (1,len,len,20)
#             #print('After model1:',k.shape)
# #        x = torch.permute(x,(0,3,1,2))
#             k = self.model2(k.squeeze()).unsqueeze(0)
# #            print('after model2:',k.shape,k)
#             #k = torch.mean(k,dim=1)
#             #out: (1,len,20)
#             #print('After model2:',k.shape)
#             k = self.model2_1(k)
# #            print('after model2_1:',k.shape, k)
# #            print('type before cat',k.dtype)
#             #out: (1,20)
#             #print('After attention:',x)
# #            k = torch.mean(attn_output,dim=1)
#             #out: (1,20)
# #            print('x device',x.device,'k device',k.device)
#             x = torch.cat((x, k), dim=0)
#            print('after cat:',x.shape, x)
#        print('type before cat',x.dtype)
        #print('After mean:',x)
        #out (batch_size, 20)
#        x = F.normalize(x)
#         out = self.model3(x)
#        bcr = bcr.to(device)
#        x1 = torch.cat((bcr,x),dim =1)
#        print('before model3:',x.shape,x)
#        out = self.model3(x1)
#        x1 = self.model4(x1)
#        x2 = torch.cat((bcr,x1),dim =1)
#        out = self.model3(x2)
#        print('after model3:',out.shape,out)
#        x = x/T
#        print(x.shape)
#        x = self.sigmoid(x)
#        out = torch.exp(x) ###I wanna try this because the output of model is too close to each other.
        return(out.squeeze())

def loss_function(out_b,out_w):
    loss = torch.sum(torch.selu(out_b-out_w+ T)+LAMBDA*(out_b*out_b + out_w*out_w))
    return loss/out_b.numel()


# In[267]:


# model = mix_model().to(device)
# test = torch.randn((10,100,100,318)).to(device)
# out = model(test)
# out.shape


# In[268]:


# out_antigen.shape
# output_tensor = F.interpolate(out_antigen.squeeze(0), size=50, mode='linear', align_corners=True).unsqueeze(0)


# In[269]:


# model2 = SelfAttentionPooling(input_dim = 100,hidden_dim = 50)
# self1_antigen = model2(out_antigen)
# print(self1_antigen.shape)
# output_tensor = F.interpolate(out_antigen, size=50, mode='linear', align_corners=True)
# print(output_tensor.shape)
# model3 = SelfAttentionPooling(input_dim = 50,hidden_dim = 20)
# out=model3(output_tensor)
# print(out.shape)


# In[270]:


# class mix_model(nn.Module):
#     def __init__(self):
#         super(mix_model,self).__init__()
#         self.model1 = nn.Sequential(
#             nn.Conv2d(in_channels=318, out_channels=128, kernel_size=3, stride=2, padding=1),
#             # in (batch size, 318,100,100)
#             # out (batch size, 128,50,50)
#             nn.ReLU(),
#             nn.MaxPool2d(kernel_size=2, stride=2, padding=0),
#             # in (batch size, 128,100,100)
#             # out (batch size, 128,25,25)
#         )
#         self.model2 = nn.Sequential(
#             nn.Conv2d(in_channels=128, out_channels=64, kernel_size=3, stride=2, padding=0),
#             # out (batchsize,64,12,12)
#             nn.ReLU(),
#             nn.MaxPool2d(kernel_size=2, stride=2, padding=0),
#             # out  (batchsize,64,6,6)
#         )
# #             nn.Conv2d(in_channels=318, out_channels=144, kernel_size=3 stride=2 padding=1)
# #             # in (batch_size,318,100,100)
# #             # out (1,144,50,50)
# #             nn.ReLU(),
# #             nn.MaxPool2d(kernel_size=2, stride = 2,padding =0)
#             # in (batch,144,50,50)
#             # out ()
# #             nn.Linear(174,100),#.to(torch.float64),
# #             # in (batch_size,100,100,318)
# #             # out (1,len,len.100)
# #             nn.ReLU(),
# #             nn.Linear(100,50),#.to(torch.float64),
# #             # out: (1,len,len,50)
# #             nn.ReLU(),
# #             nn.Linear(50,20),#.to(torch.float64),
# #             # out (1,len,len,50)
# #             nn.ReLU(),
# #             nn.Dropout(p=0.2),
# #        )
# #         self.model2 = SelfAttentionPooling(input_dim=20,hidden_dim = 10)
# #         self.model2_1 = SelfAttentionPooling(input_dim=20,hidden_dim = 10)
#         # input_dim = hidden size (number of channels)
# #        self.model2 = nn.MultiheadAttention(embed_dim=20,num_heads=N_HEADS,batch_first=True)
#         self.model3 = nn.Sequential(
#             nn.Linear(64 * 6 * 6, 512),
#             nn.ReLU(),
# #            nn.Dropout(p=0.2),
#             nn.Linear(512, 20),
# #            nn.Linear(50,20),
# #            nn.ReLU(),
# #             nn.Linear(20,10),
# #             nn.ReLU(),
# #            nn.Dropout(p=0.2),
# # #            nn.SELU(),
# #             nn.Linear(10,1),
#              nn.Dropout(p=0.2),
#         )
#         self.model4 = nn.Sequential(
#             nn.Linear(20,10),
#             nn.ReLU(),
# #            nn.Dropout(p=0.2),
#             nn.ReLU(),
#             nn.Linear(10,1),
#             nn.Dropout(p=0.2),
#         )

# #        self.sigmoid = nn.Sigmoid()
#     def forward(self,pair): ###because in getitem, return is .cuda(), now input is on gpu
# #         x = torch.empty(0)
# #         x = x.to(device)
#         x = pair.permute(0,3,1,2)
#         x = self.model1(x)
#         x = self.model2(x)
# #        print(x.shape)
#         x = x.reshape(x.shape[0], -1)
# # #         for i in range(len(pair)):
# # #             k = pair[i]
# # #             k = k.permute(0,3,1,2)
# # #             k = self.model1(k)
# # #             k = self.model2(k.squeeze()).unsqueeze(0)
# # #             k = self.model2_1(k)
# # #             x = torch.cat((x, k), dim=0)
#         x = self.model3(x)
#         out = self.model4(x)
#         return(out.squeeze())

# def loss_function(out_b,out_w):
#     loss = torch.sum(torch.selu(out_b-out_w+ T)+LAMBDA*(out_b*out_b + out_w*out_w))
#     return loss/out_b.numel()


# In[271]:


# model_mix = mix_model()
# for batch in train_loader:
#     b_pair,w_pair,idx,antigen = batch
#     print(b_pair.shape)
#     b_out=model_mix(b_pair)
#     break


# In[272]:


# test = torch.randn((10,100,100,318))
# x = test.permute(0,3,1,2)
# model1 = nn.Sequential(
#     nn.Conv2d(in_channels=318, out_channels=128, kernel_size=3, stride=2, padding=1),
#     # in (batch size, 318,100,100)
#     # out (batch size, 128,50,50)
#     nn.ReLU(),
#     nn.MaxPool2d(kernel_size=2, stride=2, padding=0),
#     # in (batch size, 128,100,100)
#     # out (batch size, 128,25,25)
# )
# model2 = nn.Sequential(
#     nn.Conv2d(in_channels=128, out_channels=64, kernel_size=3, stride=2, padding=0),
#     # out (batchsize,64,12,12)
#     nn.ReLU(),
#     nn.MaxPool2d(kernel_size=2, stride=2, padding=0),
#     # out  (batchsize,64,6,6)
# )
# x = model1(x)
# print(x.shape)
# t = model2(x)
# print(t.shape)
# t1=t.reshape(t.shape[0], -1)
# print(t1.shape)

# model3 = nn.Sequential(
#     nn.Linear(64 * 6 * 6, 512),
#     nn.ReLU(),
# #            nn.Dropout(p=0.2),
#     nn.Linear(512, 20),
# #            nn.Linear(50,20),
# #            nn.ReLU(),
# #             nn.Linear(20,10),
# #             nn.ReLU(),
# #            nn.Dropout(p=0.2),
# # #            nn.SELU(),
# #             nn.Linear(10,1),
#      nn.Dropout(p=0.2),
# )
# model4 = nn.Sequential(
#     nn.Linear(20,10),
#     nn.ReLU(),
# #            nn.Dropout(p=0.2),
#     nn.ReLU(),
#     nn.Linear(10,1),
#     nn.Dropout(p=0.2),
# )
# t = model3(t1)
# print(t.shape)
# t = model4(t)
# print(t.shape)


# In[273]:


# class SelfAttentionPooling(nn.Module):
#     """
#     Implementation of SelfAttentionPooling
#     Original Paper: Self-Attention Encoding and Pooling for Speaker Recognition
#     https://arxiv.org/pdf/2008.01077v1.pdf
#     """
#     def __init__(self, input_dim,hidden_dim):
#         super(SelfAttentionPooling, self).__init__()
# #        hidden_dim=10
#         self.W1 = nn.Linear(input_dim, hidden_dim)
#         self.relu = nn.ReLU()
#         self.W2 = nn.Linear(hidden_dim, 1)

#     def forward(self, batch_rep):
#         """
#         input:
#             batch_rep : size (N, T, H), N: batch size, T: sequence length, H: Hidden dimension

#         attention_weight:
#             att_w : size (N, T, 1)

#         return:
#             utter_rep: size (N, H)
#         """
#         att_w = F.softmax(self.W2(self.relu(self.W1(batch_rep))).squeeze(-1),dim=1).unsqueeze(-1)

#         utter_rep = torch.sum(batch_rep* att_w, dim=1)

#         return utter_rep

# class mix_model(nn.Module):
#     def __init__(self):
#         super(mix_model,self).__init__()
#         self.model1 = nn.Sequential(
#             nn.Linear(174,100),#.to(torch.float64),
#             # in (1,len,len,144)
#             # out (1,len,len.100)
#             nn.ReLU(),
#             nn.Linear(100,50),#.to(torch.float64),
#             # out: (1,len,len,50)
#             nn.ReLU(),
#             nn.Linear(50,20),#.to(torch.float64),
#             # out (1,len,len,50)
#             nn.ReLU(),
#             nn.Dropout(p=0.2),
#         )
#         self.model2 = SelfAttentionPooling(input_dim=20,hidden_dim = 10)
#         self.model2_1 = SelfAttentionPooling(input_dim=20,hidden_dim = 10)
#         # input_dim = hidden size (number of channels)
# #        self.model2 = nn.MultiheadAttention(embed_dim=20,num_heads=N_HEADS,batch_first=True)
#         self.model3 = nn.Sequential(
# #            nn.Linear(50,20),
# #            nn.ReLU(),
#             nn.Linear(20,10),
#             nn.ReLU(),
#             nn.Dropout(p=0.2),
# #            nn.SELU(),
#             nn.Linear(10,1),
#             nn.Dropout(p=0.2),
#         )
# #         self.model4 = nn.Sequential(
# #             nn.Linear(50,20),
# #             nn.ReLU()
# #         )

# #        self.sigmoid = nn.Sigmoid()
#     def forward(self,pair): ###because in getitem, return is .cuda(), now input is on gpu
#         x = torch.empty(0)
#         x = x.to(device)
#         for i in range(len(pair)):
#             k = pair[i]
#             k = k.permute(0,2,1,3)
#             k = self.model1(k)
#             k = self.model2(k.squeeze()).unsqueeze(0)
#             k = self.model2_1(k)
#             x = torch.cat((x, k), dim=0)
#         out = self.model3(x)
#         return(out)

# def loss_function(out_b,out_w):
#     loss = torch.sum(torch.selu(out_b-out_w+ T)+LAMBDA*(out_b*out_b + out_w*out_w))
#     return loss/out_b.shape[0]


# In[274]:


def train_epoch(dataloader,model,loss_fn,optimizer,verbose = False):
    train_loss = 0.0
    i = 0
    res = pd.DataFrame(columns=['record_id', 'antigen', 'out_b','out_w','success'])
    model.train()
    for batch in tqdm(dataloader):
        optimizer.zero_grad()
        df,loss=batch_result(batch,model,loss_fn)
        res = pd.concat([res,df],axis=0,ignore_index=True)
        train_loss += loss.item()
        accu_biased = sum(res['success'])/res.shape[0]
        grouped = res.groupby('antigen')['success'].mean()
        accu_unbiased = grouped.mean()
        i += 1
        if verbose:
            print('Batch',i,'batch_loss',loss.item(),'loss_until_now',train_loss/i,'accu_biased',accu_biased,'accu_unbiased',accu_unbiased)
        if math.isinf(loss) or math.isnan(loss):
            prob_bw = [out_b.item(),out_w.item()]
            print('ERROR: The loss is INF or NaN! '+str(prob_bw))
            break
        loss.backward()
        optimizer.step()
    return train_loss/i,accu_biased,accu_unbiased,res,grouped

def val_epoch(dataloader,model,loss_fn,verbose = False):
    val_loss = 0.0
    i = 0
    res = pd.DataFrame(columns=['record_id', 'antigen', 'out_b','out_w','success'])
    model.eval()
    for batch in tqdm(dataloader):
        df,loss=batch_result(batch,model,loss_fn)
        res = pd.concat([res,df],axis=0,ignore_index=True)
        val_loss += loss.item()
        accu_biased = sum(res['success'])/res.shape[0]
        grouped = res.groupby('antigen')['success'].mean()
        accu_unbiased = grouped.mean()
        i += 1
        if verbose:
            print('Batch',i,'batch_loss',loss.item(),'loss_until_now',train_loss/i,'accu_biased',accu_biased,'accu_unbiased',accu_unbiased)
        if math.isinf(loss) or math.isnan(loss):
            prob_bw = [out_b.tolist(),out_w.tolist()]
            print('ERROR: The loss is INF or NaN! '+str(prob_bw))
            break
    return val_loss/i,accu_biased,accu_unbiased,res,grouped

# def val_epoch_export(dataloader,model,loss_fn,verbose = False):
#     loss_ls = []
#     success_ls = []
#     idx_ls = []
#     batch_ls = []
#     antigen_ls = []
#     i = 0
#     model.eval()
#     for batch in dataloader:
#         if verbose:
#             print(' Batch:',i)
#             before_enter_model = time.time()
#         batch_ls.append(i)
#         better, worse, better_bcr, worse_bcr, idx, antigen_name = batch
#         idx_ls.append(idx[0])
#         antigen_ls.append(antigen_name[0])
#         out_b = model(better,better_bcr)
#         out_w = model(worse,worse_bcr)
#         loss = loss_fn(out_b,out_w)
#         loss_ls.append(loss.item())
#         if BATCH_SIZE == 1:
#             success = int(torch.gt(out_w,out_b))
#         else:
#             success = torch.sum(torch.gt(out_w,out_b)).item()/out_b.shape[0]
#         success_ls.append(success)
# #        idx_ls.append[idx]
#         if verbose:
#             print('record_id',idx[0],':')
#             print(Fore.BLUE + 'out_b: '+str(out_b.item()))
#             print(Fore.BLUE + 'out_w: '+str(out_w.item())+Style.RESET_ALL)
#             modeling_compelte = time.time()
#             #print('modeling time:'+str(modeling_compelte - before_enter_model))
#             print('loss of the batch:',loss.item())
#             print('Success or not of the current batch:',success)
#         i += 1
#         print('number of success / number of records included:',sum(success_ls)/i)
# #         if i > 30000:
# #             break
#         df = pd.DataFrame({'batch': batch_ls, 'idx': idx_ls, 'antigen': antigen_ls, 'loss': loss_ls, 'success': success_ls})
#     return df

def summary_cohort(*dfs):
    if len(dfs) == 1:
        df = dfs[0]
    elif len(dfs) > 1:
        df = pd.concat(dfs)
    else:
        print('No Dataframe provided!')
        exit()
    df[['cohort','antigen']]= df['antigen'].str.split('/',1,expand = True)
    return df.groupby('cohort')['success'].mean()


# In[275]:





#prob_ls = prob_ls = ['Q83WB0','P08069.1']
# prob_ls = ['P04275.3', 'P03303.3', 'P10636.4', 'Q5SHQ8.1', 'P12490.1', 'P05106.2', 'Q58HT7.1', 'O15409.2', 'P05981.1', 'P05556.2', 'P10646.1', 'P07946.2', 'P17763.2', 'P00749.2', 'Q2N0S5', 'P03300.3', 'P09616.2', 'P06213.4', 'Q83WB0', 'P10721.1', 'P00451.1', 'P18177.3', 'P13201.1', 'P16471.1', 'P06756.2', 'Q71F56.1', 'P10451.1', 'P00740.2', 'P13664.1', 'P01901.1','P30443.1', 'P00722.2', 'P97287.3', 'P16586.1', 'P07564.2', 'Q82446.1', 'O92972.2', 'P16154.2', 'Q81971', 'O14672.1', 'P08069.1','P03951.1','P06213.4','P08069.1']
#input_df = input_df[~input_df['id'].isin(prob_ls)]
#print(input_df.shape[0]/total_N_row*100,'% records remained after filtering bad antigens!')


# In[276]:






#print(input_df.shape[0]/3629600*100,'% remained.')


# In[277]:




# Vh_dict = build_BCR_dict(input_df[['BetterBCR_Vh','BetterBCR_CDR3h']],'Vh',precise = False)
# CDR3h_dict = build_BCR_dict(input_df[['BetterBCR_Vh','BetterBCR_CDR3h']],'CDR3h',precise = False)
# with open('/project/DPDS/Wang_lab/shared/BCR_antigen/data/toInput/easy_V_dict.pkl','wb') as f:
#     pickle.dump(Vh_dict, f)
# with open('/project/DPDS/Wang_lab/shared/BCR_antigen/data/toInput/easy_CDR3_dict.pkl','wb') as f:
#     pickle.dump(CDR3h_dict, f)


# In[278]:




# Vh_hard_dict = build_BCR_dict(exVal,'Vh',precise = True)
# CDR3h_hard_dict = build_BCR_dict(exVal,'CDR3h',precise = True)
# with open('/project/DPDS/Wang_lab/shared/BCR_antigen/data/toInput/hard_V_dict.pkl', 'wb') as f:
#     pickle.dump(Vh_hard_dict, f)
# with open('/project/DPDS/Wang_lab/shared/BCR_antigen/data/toInput/hard_CDR3_dict.pkl', 'wb') as f:
#     pickle.dump(CDR3h_hard_dict, f)


# In[279]:




exVal = exVal.reindex(columns=input_df.columns)
input_cb = pd.concat([input_df, exVal], axis=0)

aalens_df = input_cb.groupby('antigen_index')['aalens'].mean().astype(int)
len_dict = aalens_df.to_dict()
#len_dict


# In[280]:


#BATCH_SIZE = 50


# In[281]:


inputloader = InLoader(input_df, NPY_DIR,CDR3h_dict,Vh_dict,len_dict,train_ratio=0.7, subsample_ratio=IN_SUBSAMPLE, seed=None)

# Get the train and validation datasets
train_dataset = inputloader.get_train_dataset()
val_dataset = inputloader.get_val_dataset()
ex_val_dataset = SampleDataset(exVal, NPY_DIR, CDR3h_hard_dict, Vh_hard_dict, len_dict, subsample_ratio=EX_SUBSAMPLE)

# Access the group_ids for the train and validation datasets
train_group_ids = train_dataset.group_ids
val_group_ids = val_dataset.group_ids
ex_group_ids = ex_val_dataset.group_ids

# Create the PyTorch data loaders with the GroupShuffleSampler
train_loader = DataLoader(train_dataset, BATCH_SIZE, sampler=GroupShuffleSampler(train_dataset, train_group_ids))
val_loader = DataLoader(val_dataset, BATCH_SIZE, sampler=GroupShuffleSampler(val_dataset, val_group_ids))
ex_loader = DataLoader(ex_val_dataset, 1, sampler=GroupShuffleSampler(ex_val_dataset,ex_group_ids))


# In[282]:


# Shan_N439 = train_dataset.extract_antigen('Shan/N439K')
# Shan_N439.shape
# example_tensor = torch.from_numpy(Shan_N439)
# # input_tensor = example_tensor.squeeze(0)


# In[283]:


# example_permute = example_tensor.permute(0,1,2,3)
# print(example_permute[0,:,:,:].shape)
# pooling_layer = nn.AdaptiveAvgPool1d(50)
# pooling_layer(example_permute[0,:,:,:]).shape


# In[284]:




model_mix = mix_model()#.half()
model_mix
#model_mix.to(device)
#checkpoint = torch.load('/project/DPDS/Wang_lab/shared/BCR_antigen/data/output/2023-02-18_12-33-25_tag2/model_comb/Batch50BatchNumber_600Lr0.01Epoch54tag2_easy_neg.pth')
#checkpoint = torch.load('/project/DPDS/Wang_lab/shared/BCR_antigen/data/output/'+ MODEL)
#model_mix.load_state_dict(checkpoint)
# for name, param in model_mix.named_parameters():
#     print('model parameters dtype:')# print the name and data type of each parameter
#     print(f"{name}: {param.dtype}")
# for param in model_mix.parameters():
#     param.data = param.data.half()
#     if param.grad is not None:
#         param.grad.data = param.grad.data.half()
# set the model to use float16
#model = model.half()
# model_mix.half()
model_mix.to(device)
#for name, param in model_mix.named_parameters():
#    print(name, param.dtype)
optimizer = torch.optim.Adam(model_mix.parameters(),lr=LR)
scheduler = ReduceLROnPlateau(optimizer,mode = 'min',factor=0.1,patience =5,verbose=False,threshold_mode='rel',threshold=1e-4)


# In[285]:


# i = 0
# loop_loss = 0.0
# res = pd.DataFrame(columns=['record_id', 'antigen', 'out_b','out_w','success'])
# for batch in val_loader:
#     model_mix.train()
#     df,loss=batch_result(batch,model_mix,loss_function)
# #     b_pair,w_pair,index_idx,antigen_index =batch
# # #    print(zip(index_idx,antigen_index))
# # #    print(b_pair.shape,w_pair.shape)
# #     out_b = model_mix(b_pair)
# #     out_w = model_mix(w_pair)
# #     success = torch.gt(out_b,out_w).tolist()
# #     outb = out_b.tolist()
# #     outw = out_w.tolist()
# #     df = pd.DataFrame({'record_id':index_idx,'antigen':antigen_index,'out_b':outb,'out_w':outw,'success':success})
#     res = pd.concat([res,df],axis=0,ignore_index=True)
#     loop_loss += loss.item()
#     accu_biased = sum(res['success'])/res.shape[0]
#     grouped = res.groupby('antigen')['success'].mean()
#     accu_unbiased = grouped.mean()
#     i += 1
#     print('batch_loss',loss.item(),'loss_until_now',loop_loss/i,'accu_biased',accu_biased,'accu_unbiased',accu_unbiased)
#     #batch_results =  {key: (antigen, outb, outw, success) for key, antigen,outb,outw, success in zip(index_idx, antigen_index, outb, outw, torch.gt(out_b,out_w).tolist())}
#     #total_results = {**total_results,**batch_results}
#     i += 1
#     loss.backward()
#     optimizer.step()
#     if i > 10:
#         break
#print(res)


# In[286]:


# print(res)
# grouped = res.groupby('antigen')['success'].mean()
# print(grouped)
# grouped.mean()


# In[287]:


# def loss_function(out_b,out_w):
#     loss = torch.sum(torch.selu(out_b-out_w+ T)+LAMBDA*(out_b*out_b + out_w*out_w))
#     return loss/out_b.shape[0]


# In[288]:


# val_loss = 0.0
# i = 0
# res = pd.DataFrame(columns=['record_id', 'antigen', 'out_b','out_w','success'])
# model_mix.eval()

# for batch in ex_loader:
#     b_pair,w_pair,index_idx,antigen_index =batch
# #    print(zip(index_idx,antigen_index))
# #    print(b_pair.shape,w_pair.shape)
#     out_b = model_mix(b_pair)
#     out_w = model_mix(w_pair)
#     print(out_b)
# #     loss = loss_fn(out_b,out_w)
# #     success = torch.gt(out_b,out_w).tolist()
# #     outb = out_b.tolist()
# #     outw = out_w.tolist()
# #     df = pd.DataFrame({'record_id':index_idx,'antigen':antigen_index,'out_b':outb,'out_w':outw,'success':success})
#     break


# In[289]:


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


# In[290]:


with open(out_dir+"/para.txt", "w") as file:
    # Write the text to the file
    file.write(f"batch size: {BATCH_SIZE}\ninput label:{INPUT}\nEpoch max:{EPOCH}\nLearning rate:{LR}\nCut off:{CUTOFF}\ntag:{TAG}\ndelta:{T}\nLambda:{LAMBDA}\nweight:{w}\nlimit:{LIMIT}\nverbose:{VERBOSE}\ngroup_size:{GROUP_SIZE}\nmodel:{MODEL}\nmax number of antigen:{Max_num}\number of channels:{CHANNEL_ANTIGEN}\nsmall_sample:{Small_sample}\ninter_sub:{IN_SUBSAMPLE}\nouter_sub:{EX_SUBSAMPLE}")


# In[ ]:


print('batch size:',BATCH_SIZE,'EPOCH:', EPOCH,'Learning rate:',LR)
start_time = time.time()
print('Initialling ...')
init_train_loss,init_train_accu_biased,init_train_accu_unbiased,*_ = val_epoch(train_loader,model_mix,loss_function,verbose = VERBOSE)
init_val_loss,init_val_accu_biased, init_val_accu_unbiased,*_ = val_epoch(val_loader,model_mix,loss_function,verbose = VERBOSE)
init_ex_loss,init_ex_accu_biased,init_ex_accu_unbiased,res_init,grouped_ini = val_epoch(ex_loader,model_mix,loss_function,verbose = VERBOSE)
print('Initial: '+'Epoch: -1 train loss: '+str(init_train_loss)+' val loss: '+str(init_val_loss)+' val accu biased:'+str(init_val_accu_biased)+' val accuracy: '+str(init_val_accu_unbiased))
print('Initial external val loss: '+str(init_ex_loss)+' ex average accuracy biased: '+str(init_ex_accu_biased)+' ex accuracy unbiased:'+str(init_ex_accu_unbiased))
print('Initial training:',time.time()-start_time)


# In[ ]:


res_init.to_csv(out_dir+'/initial_accu_all.csv')
grouped_ini.to_csv(out_dir+'/initial_acc_per.csv')


# In[ ]:


train_LOSS = [init_train_loss]
train_ave_SUCCESS = [init_train_accu_biased]
train_ACC = [init_train_accu_unbiased]
val_LOSS = [init_val_loss]
val_ave_SUCCESS = [init_val_accu_biased]
val_ACC = [init_val_accu_unbiased]
ex_LOSS = [init_ex_loss]
ex_ave_SUCCESS = [init_ex_accu_biased]
ex_ACC = [init_ex_accu_unbiased]


# In[ ]:


print('Start Training...')
for epoch in range(EPOCH):
    inputloader.new_epoch()
    print('epoch:',str(epoch))
    start_epoch_time = time.time()
    train_loss,train_accu_biased,train_accuracy,*_ = train_epoch(train_loader,model_mix,loss_function,optimizer,verbose = VERBOSE)
    val_loss,val_accu_biased,val_accuracy,*_ = val_epoch(val_loader,model_mix,loss_function,verbose = VERBOSE)
    ex_loss,ex_accu_biased,ex_accuracy,res_ex,group_ex = val_epoch(ex_loader,model_mix,loss_function,verbose=VERBOSE)
#    accu_block_df = pd.DataFrame(accu_dict_ex)
    if epoch%5 == 0:
        res_ex.to_csv(out_dir+'/ex_accuracy_tag'+str(TAG)+'_Epoch'+str(epoch)+'_Lr'+str(LR)+'_all.csv')
        group_ex.to_csv(out_dir+'/ex_accuracy_tag'+str(TAG)+'_Epoch'+str(epoch)+'_Lr'+str(LR)+'_per_antigen.csv')
    print('Epoch '+str(epoch)+' train loss: '+str(train_loss)+' val loss: '+str(val_loss)+' val accu biased: '+str(val_accu_biased)+' val accuracy: '+str(val_accuracy)+'\n'
          +' external val loss: '+str(ex_loss)+' ex accu biased: '+str(ex_accu_biased)+' external val accuracy: '+str(ex_accuracy))
    scheduler.step(val_loss) ##should I try scheduler.step(-val_accuracy)?
    train_LOSS.append(train_loss)
    train_ave_SUCCESS.append(train_accu_biased)
    train_ACC.append(train_accuracy)
    val_LOSS.append(val_loss)
    val_ave_SUCCESS.append(val_accu_biased)
    val_ACC.append(val_accuracy)
    ex_LOSS.append(ex_loss)
    ex_ave_SUCCESS.append(ex_accu_biased)
    ex_ACC.append(ex_accuracy)
    state_dict=model_mix.state_dict()
    torch.save(state_dict,model_dir+"/Batch"+str(BATCH_SIZE)+"Lr"+str(LR)+"Epoch"+str(epoch)+"tag"+str(TAG)+"_3diff_v10.pth")
    # trian_index_file = model_dir+"/Batch"+str(BATCH_SIZE)+"Lr"+str(LR)+"Epoch"+str(epoch)+"tag"+str(TAG)+"_train_index.csv"
    # val_index_file = model_dir+"/Batch"+str(BATCH_SIZE)+"Lr"+str(LR)+"Epoch"+str(epoch)+"tag"+str(TAG)+"_val_index.csv"
    # ex_index_file = model_dir+"/Batch"+str(BATCH_SIZE)+"Lr"+str(LR)+"Epoch"+str(epoch)+"tag"+str(TAG)+"_ex_val_index.csv"
    # with open(trian_index_file, 'w', newline='') as file:
    #     writer = csv.writer(file)
    #     # Write each item in the list as a separate row in the CSV file
    #     for item in train_idx:
    #         writer.writerow([item])
    # with open(val_index_file, 'w', newline='') as file:
    #     writer = csv.writer(file)
    #     # Write each item in the list as a separate row in the CSV file
    #     for item in val_idx:
    #         writer.writerow([item])
    # with open(ex_index_file, 'w', newline='') as file:
    #     writer = csv.writer(file)
    #     # Write each item in the list as a separate row in the CSV file
    #     for item in ex_idx:
    #         writer.writerow([item])
    if math.isnan(train_loss) or math.isnan(val_loss):
        state_dict = model_mix.state_dict()
        torch.save(state_dict,model_dir+"/Error_model_epoch"+str(epoch)+"/Batch"+str(BATCH_SIZE)+"_Lr"+str(LR)+"_tag"+str(TAG)+".pth")
        break
        print('Epoch:',epoch,'train_loss:',train_LOSS[-1])
#     if val_accuracy > 0.98:
#         print('Validation Accuracy reached 0.98!')
#         break
    print('All_through time: '+str(time.time()-start_epoch_time))


# In[ ]:


# print("Model parameters:")
# for name, param in model_mix.named_parameters():
#     print(f"{name}: {param}")


# In[ ]:


# for batch in ex_loader:
#     better_pair,worse_pair,index_id,antigen =batch
#     out_b = model_mix(better_pair)
#     out_w = model_mix(worse_pair)
#     print(out_b,out_w)
#     break


# In[ ]:




#STOP HERE and DO NOT OVERWRITE lossTable!! ATTACH THE NEW RESUlTS to it!
epoch_ls = [i for i in range(len(train_LOSS)-1)]
epoch_ls.insert(0,-1)
lossTable =pd.DataFrame({'Epoch':epoch_ls,'Train_loss':train_LOSS,'Train_accuracy':train_ACC,
                        'Val_loss':val_LOSS,'Val_biased_accuracy':val_ave_SUCCESS,'Val_accuracy':val_ACC,
                        'Ex_val_loss':ex_LOSS,'Ex_biased_accuracy':val_ave_SUCCESS,'Ex_val_accuracy':ex_ACC})#Sequential


# In[ ]:


print(lossTable)


# In[ ]:




fig, axs = plt.subplots(3, 3)
for ax in axs.flat:
    ax.tick_params(axis='both', which='both', left=True, labelleft=True, bottom=True, labelbottom=True)
lossTable.plot(x='Epoch',y='Train_loss',kind='line',ax=axs[0,0])
lossTable.plot(x='Epoch',y='Train_accuracy',kind='line',ax=axs[0,1])
lossTable.plot(x='Epoch',y='Val_loss',kind='line',ax=axs[1,0])
lossTable.plot(x='Epoch',y='Val_biased_accuracy',kind='line',ax=axs[1,1])
lossTable.plot(x='Epoch',y='Val_accuracy',kind='line',ax=axs[1,2])
lossTable.plot(x='Epoch',y='Ex_val_loss',kind='line',ax=axs[2,0])
lossTable.plot(x='Epoch',y='Ex_biased_accuracy',kind='line',ax=axs[2,1])
lossTable.plot(x='Epoch',y='Ex_val_accuracy',kind='line',ax=axs[2,2])
plt.savefig(out_dir+'/lossPlot_tag'+str(TAG)+'_Epoch'+str(EPOCH)+'_Lr'+str(LR)+'.png',dpi=200)
lossTable.to_csv(out_dir+'/lossTable_tag'+str(TAG)+'_Epoch'+str(EPOCH)+'_Lr'+str(LR)+'.csv')


# In[ ]:




# from tqdm import tqdm
# import time

# # A list of items you want to iterate through
# #items = range(100)

# # Wrap the iterable with tqdm to display the progress bar
# i = 0
# for batch in tqdm(val_loader):
#     if VERBOSE:
#         print(' Batch:',i)
#         before_enter_model = time.time()
# #    batch_ls.append(i)
#     better, worse, better_bcr, worse_bcr, idx, antigen_name = batch
# #    idx_ls.append(idx[0])
# #    antigen_ls.append(antigen_name[0])
#     out_b = model_mix(better,better_bcr)
#     out_w = model_mix(worse,worse_bcr)
#     loss = loss_function(out_b,out_w)
# #    loss_ls.append(loss.item())
#     if BATCH_SIZE == 1:
#         success = int(torch.gt(out_w,out_b))
#     else:
#         success = torch.sum(torch.gt(out_w,out_b)).item()/out_b.shape[0]
# #    success_ls.append(success)
#     if math.isinf(loss) or math.isnan(loss):
#         prob_bw = [out_b,out_w]
#         print('ERROR: The loss is INF or NaN! '+str(prob_bw))
#         break
# #    loss.backward()
# #    optimizer.step()
# #    train_loss += loss.item()
#     i += 1
# for batch in dataloader:
#     if verbose:
#         print(' Batch:',i)
#         before_enter_model = time.time()
#     batch_ls.append(i)
#     better, worse, better_bcr, worse_bcr, idx, antigen_name = batch
#     idx_ls.append(idx[0])
#     antigen_ls.append(antigen_name[0])
#     out_b = model(better,better_bcr)
#     out_w = model(worse,worse_bcr)
#     loss = loss_fn(out_b,out_w)
#     loss_ls.append(loss.item())
#     if BATCH_SIZE == 1:
#         success = int(torch.gt(out_w,out_b))
#     else:
#         success = torch.sum(torch.gt(out_w,out_b)).item()/out_b.shape[0]
#     success_ls.append(success)
