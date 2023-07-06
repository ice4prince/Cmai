#!/usr/bin/env python
# coding: utf-8

# In[2]:


import os
import sys
import pandas as pd
import numpy as np
import random
from torch.utils.data import DataLoader, TensorDataset, SubsetRandomSampler
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


# In[3]:


torch.cuda.empty_cache()


# In[4]:


#os.getcwd()
os.chdir('/project/DPDS/Wang_lab/shared/BCR_antigen/code/DeLAnO/')

from wrapV.Vwrap import embedV
from wrapCDR3.CDR3wrap import embedCDR3


# In[3]:


# Create an ArgumentParser object
parser = argparse.ArgumentParser(description='Parameters for pair model.')

# Add a optional argument
parser.add_argument('--batch_size', type=int, help='the batch size, default is 100',default = 100)
parser.add_argument('--number_of_batch', type=int, help='the number of batch in each epoch, default is 300',default = 300)
parser.add_argument('--epoch', type=int, help='the maximum of epoch, default is 50',default = 50)
parser.add_argument('--tag', type=int, help='the tag to specify the running',default = 1)
parser.add_argument('--initial_learning_rate', type=float, help='the starting leanring rate, default is 0.05', default=0.05)
parser.add_argument('--delta', type=float, help='the minimum difference between better and worse, default is 0.01', default=0.01)
parser.add_argument('--Lambda', type=float, help='the parameter adjustable for loss to avoid INF and NAN', default=0.03)
parser.add_argument('--weight', type=float, help='the parameter adjustable to shrink embedding of antigens to input', default=0.03)
parser.add_argument('--cut_off', type=int, help='the maximum of length of sequencing, over which will be deleted. default is 1800',default = 1800)
parser.add_argument('--verbose', action='store_true', help='Enable verbose output')

# Add a switch (flag) argument
#parser.add_argument('--verbose', action='store_true', help='Print verbose output')
# Parse the command-line arguments
args = parser.parse_args()


# In[4]:


BATCH_SIZE = args.batch_size
NUMBER_BATCH = args.number_of_batch
LR = args.initial_learning_rate
EPOCH = args.epoch
CUTOFF = args.cut_off
TAG = args.tag
T = args.delta
LAMBDA = args.Lambda
w = args.weight
NPY_DIR ='BATCH_SIZEoject/DPDS/Wang_lab/shared/BCR_antigen/data/rfscript/Antigen_embed/CLEANED/'


# In[5]:


BATCH_SIZE = 50
NUMBER_BATCH = 20
LR = 0.01
EPOCH = 50
CUTOFF = 1800
TAG = 1
T = 0.01
LAMBDA = 0
w = 1
LIMIT = 10
NPY_DIR ='/project/DPDS/Wang_lab/shared/BCR_antigen/data/rfscript/Antigen_embed/CLEANED/'
VERBOSE = True


# In[6]:


device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")
print(device)


# In[7]:


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


# In[8]:


def extract_antigen(antigen_dict,antigen_index,npy_dir = NPY_DIR,verbose = False):
    if antigen_index in antigen_dict:
        embedding = antigen_dict[antigen_index]
    else:
        
        try:
            embedding = np.load(str(npy_dir+antigen_index+'.pair.npy'))/w
#            print(npy.shape)
            antigen_dict[antigen_index] = embedding
        except ValueError:
            print('The embedding of antigen %s cannot be found!' % embeded_Antigen_Files[i])
        if verbose:
            print(Fore.RED + 'New antigen added to dictionary:',antigen_index)
            print(Fore.RED + 'Number of antigens included:',len(antigen_dict))
            print(Style.RESET_ALL)
    return embedding


# In[9]:


def build_BCR_dict(dataset,colname):
    cols = dataset.filter(like = colname)
    uniq_keys = pd.unique(cols.values.ravel()).tolist()
    if colname == 'CDR3h':
        uniq_embedding,*_ = embedCDR3(uniq_keys)
    elif colname == 'Vh':
        uniq_embedding,*_ = embedV(uniq_keys)
    i = 0
    mydict = {}
    for key in uniq_keys:
        mydict[key] = uniq_embedding[i]
        i += 1
    return(mydict)


# In[10]:


def embedding_BCR(dataset,CDR3_dict,V_dict):
    BCR_b = np.empty((0,30),dtype=np.float32)
#    print('BCR_b intial type:',BCR_b.dtype)
    BCR_w = np.empty((0,30),dtype=np.float32)
    for i in dataset.index:
        V_better = V_dict[dataset.loc[i,'BetterBCR_Vh']]
        CDR3_better = CDR3_dict[dataset.loc[i,'BetterBCR_CDR3h']]
#        print('V_better:',V_better.dtype,'CDR3_better:',CDR3_better.dtype)
        i_b = np.concatenate((CDR3_better,V_better),axis=0)
#        print('i_b:',i_b.dtype)
#        print('i_b reshape:',i_b.reshape(1,30).dtype)
        BCR_b = np.concatenate((BCR_b,i_b.reshape(1,30)),axis=0)
#        print('BCR_b:',BCR_b.dtype)
        V_worse = V_dict[dataset.loc[i,'WorseBCR_Vh']]
        CDR3_worse = CDR3_dict[dataset.loc[i,'WorseBCR_CDR3h']]
        i_w = np.concatenate((CDR3_worse,V_worse),axis=0)
        BCR_w = np.concatenate((BCR_w,i_w.reshape(1,30)),axis=0)
#        break
#     CDR3_better,*_= embedCDR3(dataset['BetterBCR_CDR3h'].tolist())
#     CDR3_worse,*_= embedCDR3(dataset['WorseBCR_CDR3h'].tolist())
#     V_Better, *_ = embedV(dataset['BetterBCR_Vh'].tolist())
#     V_Worse, *_ = embedV(dataset['WorseBCR_Vh'].tolist())
#     CDR3_b,CDR3_w,V_b,V_w = list(map(np.vstack,(CDR3_better,CDR3_worse,V_Better,V_Worse)))
#     BCR_b = np.concatenate((CDR3_b,V_b),axis=1)
#     BCR_w = np.concatenate((CDR3_w,V_w),axis=1)
    return(BCR_b,BCR_w)


# In[11]:


def predict_size(len):
    # Define the shape of the tensor
    shape = (1, len, len, 318)
    # Create a tensor with the specified shape
    tensor = torch.from_numpy(np.ones(shape))
    # Predict the size of the tensor
    num_elements = tensor.numel()
    element_size = tensor.element_size()
    tensor_size = num_elements * element_size
#    size_gb = tensor_size/(1024*1024*1024)
#    print(f"Tensor size: {size_gb:.2f} GB")
    return(tensor_size)


# In[29]:


def check_batch_size(batch,limit):
    antigen_index = batch['antigen_index'].unique()
#    print(antigen_index)
    antigen_len_dict = batch.groupby('antigen_index')['aalens'].min().to_dict()
    total_size = 0
    for value in antigen_len_dict.values():
        total_size += predict_size(value)
#    print('batch size:',total_size/(1024*1024*1024))
    return total_size<limit*1024*1024*1024


# In[13]:


def generate_index(dataset,batch_size):
#batch_size = BATCH_SIZE
#dataset = hard_neg_train.loc[hard_neg_train.index[:100]]
    nrow = dataset.shape[0]
#    print(nrow)
    l = list(range(nrow))
#    print(l)
    #random.shuffle(l)
    #n = batch_size
    x = []
    for i in range(0,len(l),batch_size):
        each = l[i:i+batch_size]
        random.shuffle(each)
        x.append(each)
#print(x[0])
    return(x)


# In[28]:


#dataset = hard_neg_train
def batch_generator(dataset,batch_size,limit,CDR3_dict,V_dict):
#    while True:
    print('generating the batch index...')
    batch_index = generate_index(dataset,batch_size)
    for i in range(len(batch_index)):
        batch = dataset.iloc[batch_index[i]]
#            batch['antigen_index']=batch['Project']+'/'+batch['id']
#            if not check_batch_size(batch,limit):
#                break
        record_id = batch['record_id'].tolist()
        antigen_set = batch['antigen_index'].unique()
        antigen_index_dict = batch.set_index('record_id')['antigen_index'].to_dict()
        BCR_b, BCR_w = embedding_BCR(batch,CDR3_dict,V_dict)
        BCR_b_g = torch.from_numpy(BCR_b).to(device)
        BCR_w_g = torch.from_numpy(BCR_w).to(device)
        BCR_b_dict = {key: BCR_b_g[i,:] for i, key in enumerate(record_id)}
        BCR_w_dict = {key: BCR_w_g[i,:] for i, key in enumerate(record_id)}
        batch_out = {'record_id':record_id,'BCR_b':BCR_b_dict,'BCR_w':BCR_w_dict,'antigen_index':antigen_index_dict,'antigen_set':antigen_set}
        yield batch_out
        del batch_out
#        else:
#            break


# In[15]:


def pick_n_delete(antigen_dict,antigen_set):
    antigen_can_be_replaced = [key for key in antigen_dict.keys() if key not in antigen_set]
    antigen_to_delete = random.choice(antigen_can_be_replaced)
    del antigen_dict[antigen_to_delete]
    print('delete antigen:',antigen_to_delete)


# In[16]:


def antigen_to_device(antigen,antigen_set,antigen_dict,limit = LIMIT):
    while True:
        try:
            gpu_antigen = torch.from_numpy(antigen).to(device)
            allocated_size = torch.cuda.memory_allocated()/(1024**3)
            if allocated_size < limit:
                break
            else:
                pick_n_delete(antigen_dict,antigen_set)
        except RuntimeError as e:
    # If the allocation fails, free up some memory and try again
            if 'CUDA out of memory' in str(e):
                torch.cuda.empty_cache()
                pick_n_delete(antigen_dict,antigen_set)
            else:
                raise e
        
    return gpu_antigen


# In[17]:


def rotate_antigens(antigen_set,dict_in,dict_out):
    for single_antigen in antigen_set:
        if not single_antigen in dict_in:
            new_antigen = extract_antigen(dict_out,single_antigen)
            dict_in[single_antigen] = antigen_to_device(new_antigen,antigen_set,dict_in)
#            dict_in[single_antigen] =  antigen_in
            print('antigens in:',dict_in.keys())
#            print('antigen dictionary size:',sys.getsizeof(dict_in)/(1024**3))
#            print('memory allocated',torch.cuda.memory_allocated()/(1024**3))
    return dict_in          


# In[18]:


def comb_embedding(antigen_name,single_BCR,antigen_dict):
    single_antigen = antigen_dict[antigen_name]
    len = aalen(single_antigen)
#    single_antigen_g = torch.from_numpy(single_antigen).to(device)
    BCR_t = torch.tile(single_BCR,(1,len,len,1))
    Cb_BCR_antigen = torch.cat((single_antigen,BCR_t),dim=3).to(device)
    return Cb_BCR_antigen


# In[30]:


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
    def forward(self,antigen_dict_in:dict,antigen_index:dict,BCR_dict:dict,record_id:list):
        x = torch.empty(0)
        x = x.to(device)
        for i in record_id:
            k = comb_embedding(antigen_index[i],BCR_dict[i],antigen_dict_in)
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


# In[146]:


# dataset = hard_neg_train.loc[hard_neg_train.index[:100]]


# In[164]:


# batch_gen = batch_generator(dataset,BATCH_SIZE,LIMIT,CDR3h_dict,Vh_dict)
# i = 0
# for batch in batch_gen:
#     print('batch',i)
#     print('antigens in this batch:',batch['antigen_set'])
#     rotate_antigens(batch['antigen_set'],antigen_dict_in,antigen_dict_out)
# #    k = comb_embedding(batch['antigen_set'][0],batch['BCR_b'][batch['record_id'][0]],antigen_dict_in)
# #    print(k.dtype)
# #    t_model(k)
# #    print(batch['record_id'])
#     out_b = model_mix(antigen_dict_in,batch['antigen_index'],batch['BCR_b'],batch['record_id'])
#     out_w = model_mix(antigen_dict_in,batch['antigen_index'],batch['BCR_w'],batch['record_id'])
# #    break
        
# #    batch['BCR_w']
#     i += 1


# In[31]:


def train_epoch(dataset,antigen_dict_in,antigen_dict_out,model,loss_fn,optimizer,batch_size,verbose = False):
    train_loss = 0.0
    i = 0
    model.train()
    batch_gen = batch_generator(dataset,batch_size,LIMIT,CDR3h_dict,Vh_dict)
    for batch in batch_gen:
        if verbose:
            print(' Batch:',i)
            print('antigens in this batch:',batch['antigen_set'])
            before_enter_model = time.time()
        rotate_antigens(batch['antigen_set'],antigen_dict_in,antigen_dict_out)   
#        b_dict, w_dict = load_single_batch(batch,antigen_dict,verbose,npy_dir = NPY_DIR)
        optimizer.zero_grad()
        out_b = model(antigen_dict_in,batch['antigen_index'],batch['BCR_b'],batch['record_id'])
        out_w = model(antigen_dict_in,batch['antigen_index'],batch['BCR_w'],batch['record_id'])
        if verbose:
            print(Fore.BLUE + 'out_b: '+str(out_b[:3]))
            print(Fore.BLUE + 'out_w: '+str(out_w[:3]))
            modeling_compelte = time.time()
            print(Style.RESET_ALL)
            print('modeling time:'+str(modeling_compelte - before_enter_model))
        loss = loss_fn(out_b,out_w)
        if verbose:
            print('loss of the batch:',loss.item())
        if math.isinf(loss) or math.isnan(loss):
            prob_bw = [out_b,out_w]
            print('ERROR: The loss is INF or NaN! '+str(prob_bw))
            break
        loss.backward()
        optimizer.step()
        train_loss += loss.item()
        i += 1
#         if i > number_of_batches:
#             break
    return batch_gen,train_loss/i


# In[73]:


def val_epoch(dataset,antigen_dict_in,antigen_dict_out,model,loss_fn,batch_size,verbose = False):
    val_loss = 0.0
    val_success = 0
    i = 0
    model.eval()
    batch_gen = batch_generator(dataset,batch_size,LIMIT,CDR3h_dict,Vh_dict)
    for batch in batch_gen:
        if verbose:
            print(' Batch:',i)
            print('antigens in this batch:',batch['antigen_set'])
            before_enter_model = time.time()
        rotate_antigens(batch['antigen_set'],antigen_dict_in,antigen_dict_out)  
#        b_dict, w_dict = load_single_batch(batch,antigen_dict,verbose,npy_dir = NPY_DIR)
        optimizer.zero_grad()
        out_b = model(antigen_dict_in,batch['antigen_index'],batch['BCR_b'],batch['record_id'])
        out_w = model(antigen_dict_in,batch['antigen_index'],batch['BCR_w'],batch['record_id'])
        if verbose:
            print(Fore.BLUE + 'out_b: '+str(out_b[:3]))
            print(Fore.BLUE + 'out_w: '+str(out_w[:3]))
            modeling_compelte = time.time()
            print(Style.RESET_ALL)
            print('modeling time:'+str(modeling_compelte - before_enter_model))
        success = torch.sum(torch.gt(out_w,out_b)).item()/out_b.shape[0]
        val_success += success
        loss = loss_fn(out_b,out_w)
        val_loss += loss.item()
        if verbose:
            print('loss of the batch:',loss.item())
        i += 1
#         if i > number_of_batches:
#             break
    return batch_gen,val_loss/i,val_success/i


# In[21]:


hard_file_train = pd.read_csv('/project/DPDS/Wang_lab/shared/BCR_antigen/data/hard_negative_train_downsampled.csv')
hard_file_val = pd.read_csv('/project/DPDS/Wang_lab/shared/BCR_antigen/data/hard_negative_validation_downsampled.csv')


# In[22]:


hard_neg_train = preprocess(hard_file_train,CUTOFF)
hard_neg_val = preprocess(hard_file_val, CUTOFF)
#hard_neg_train = hard_neg_train.assign(record_id = ['record_' + str(s) for s in range(hard_neg_train.shape[0])])
#hard_neg_val = hard_neg_val.assign(record_id = ['record_' + str(s) for s in range(hard_neg_val.shape[0])])
print('Train Data: ')
check_space(hard_neg_train)
print('Validation Data: ')
check_space(hard_neg_val)


# In[23]:


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


# In[24]:


hard_neg_df = pd.concat([hard_neg_train,hard_neg_val],axis = 0)
Vh_dict = build_BCR_dict(hard_neg_df,'Vh')
CDR3h_dict = build_BCR_dict(hard_neg_df,'CDR3h')


# In[25]:


antigen_dict_in = {}
antigen_dict_out = {}


# In[26]:


model_mix = mix_model()
#model_mix.to(device)
checkpoint = torch.load('/project/DPDS/Wang_lab/shared/BCR_antigen/data/output/2023-02-18_12-33-25_tag2/model_comb/Batch50BatchNumber_600Lr0.01Epoch54tag2_easy_neg.pth')
model_mix.load_state_dict(checkpoint)
model_mix.to(device)
#for name, param in model_mix.named_parameters():
#    print(name, param.dtype)
optimizer = torch.optim.Adam(model_mix.parameters(),lr=LR)
scheduler = ReduceLROnPlateau(optimizer,mode = 'min',factor=0.1,patience =10,verbose=False,threshold_mode='rel',threshold=1e-5)


# In[19]:


# model_mix = mix_model()
# #model_mix.to(device)
# checkpoint = torch.load('/project/DPDS/Wang_lab/shared/BCR_antigen/data/output/2023-02-18_12-33-25_tag2/model_comb/Batch50BatchNumber_600Lr0.01Epoch54tag2_easy_neg.pth')
# model_mix.load_state_dict(checkpoint)
# optimizer = torch.optim.Adam(model_mix.parameters(),lr=LR)
# scheduler = ReduceLROnPlateau(optimizer,mode = 'min',factor=0.1,patience =10,verbose=False,threshold_mode='rel',threshold=1e-5)


# In[32]:


#NUMBER_BATCH = 20
print('number of batches:',NUMBER_BATCH,' batch size: ',BATCH_SIZE,' EPOCH: ', EPOCH,' Learning rate: ',LR)
start_time = time.time()
print('Initialling ...')
antigen_dict = {}
_,initial_train_loss,_ = val_epoch(hard_neg_train,antigen_dict_in,antigen_dict_out,model_mix,loss_function,BATCH_SIZE,verbose = VERBOSE)
_,initial_val_loss,initial_val_accuracy = val_epoch(hard_neg_val,antigen_dict_in,antigen_dict_out,model_mix,loss_function,BATCH_SIZE,verbose = VERBOSE)
print(time.time()-start_time)


# In[33]:


for key in ['46472/COV_MERS', '46472/COV_OC43', '46472/COV_HKU1', '46472/COV_SARS', '46472/COV_SARS2', '46472/CZA97', '54042-4/CoV_OC43', 'Dugan/HA-H3-Kansas', '54042-4/CoV_229E', '54042-4/CoV_NL63', '46472/ZM197']:
    del antigen_dict_in[key]
    print(torch.cuda.memory_allocated()/(1024**3))


# In[20]:


train_LOSS = [initial_train_loss]
val_LOSS = [initial_val_loss]
val_ACC = [initial_val_accuracy]


# In[ ]:


for epoch in range(EPOCH):
    print('epoch:'+str(epoch))
    start_epoch_time = time.time()
    batch_gen_train,train_loss = train_epoch(hard_neg_train,antigen_dict,model_mix,loss_function,optimizer,NUMBER_BATCH,verbose = VERBOSE)
    batch_gen_val,val_loss, val_accuracy = val_epoch(hard_neg_val,antigen_dict,model_mix,loss_function,NUMBER_BATCH,verbose = VERBOSE)
    print('Epoch '+str(epoch)+' train loss: '+str(train_loss)+' val loss: '+str(val_loss)+' val accuracy: '+str(val_accuracy))
    scheduler.step(val_loss) ##should I try scheduler.step(-val_accuracy)maybe?
    train_LOSS.append(train_loss)
    val_LOSS.append(val_loss)
    val_ACC.append(val_accuracy)
    state_dict=model_mix.state_dict()
    torch.save(state_dict,model_dir+"/Batch"+str(BATCH_SIZE)+'BatchNumber_'+str(NUMBER_BATCH)+"Lr"+str(LR)+"Epoch"+str(epoch)+"tag"+str(TAG)+"_easy_neg.pth")
    id_train = model_dir+'/Index_Batch'+str(BATCH_SIZE)+'BatchNumber_'+str(NUMBER_BATCH)+"Lr"+str(LR)+"Epoch"+str(epoch)+"tag"+str(TAG)+'.txt'
    id_val = model_dir+'/Index_Batch'+str(BATCH_SIZE)+'BatchNumber_'+str(NUMBER_BATCH)+"Lr"+str(LR)+"Epoch"+str(epoch)+"tag"+str(TAG)+'.txt'
    export_index(id_train,batch_gen_train,min(hard_neg_train.shape[0]//BATCH_SIZE,NUMBER_BATCH))
    export_index(id_val,batch_gen_val,min(hard_neg_val.shape[0]//BATCH_SIZE,NUMBER_BATCH))
    #with open(index_train, "w", newline="") as f1:
    #    writer = csv.writer(f1)
    #    writer.writerows(batch_index_train)
    #with open(index_val, "w", newline="") as f2:
    #    writer = csv.writer(f2)
    #    writer.writerows(batch_index_val)
    if math.isnan(train_loss):
        state_dict = model_mix.state_dict()
        torch.save(state_dict,model_dir+"/Error_model_epoch"+str(epoch)+"/Batch"+str(BATCH_SIZE)+'_BatchNumber'+str(NUMBER_BATCH)+"_Lr"+str(LR)+"_tag"+str(TAG)+".pth")
        break
        print('Epoch:',epoch,'train_loss:',train_LOSS[-1])
    if val_accuracy > 0.97:
        print('Validation Accuracy reached 0.97!')
        break
    print('All_through time: '+str(time.time()-start_epoch_time))


# In[42]:


# def export_index(outfile,batch_gen,number_of_batches):
#     with open(outfile, 'w', newline='') as f:
#         writer = csv.writer(f)
#         # Write the header row with column names
#         writer.writerow([''] + ['record_{}'.format(i) for i in range(10)])

#         # Iterate over the generated list of lists and write them to the CSV
#         i = 0
#         for batch in batch_gen:
#             row_name = 'Epoch {}'.format(i)
#             writer.writerow([row_name] + batch['record_id'])
#             i += 1
#             if i > number_of_batches:
#                 break


# In[31]:


# batch_gen = batch_generator(hard_neg_val, 100)
# with open('my_data1.csv', 'w', newline='') as f:
#     writer = csv.writer(f)
#     # Write the header row with column names
#     writer.writerow([''] + ['record_{}'.format(i) for i in range(10)])
    
#     # Iterate over the generated list of lists and write them to the CSV
#     i = 0
#     for batch in batch_gen:
#         row_name = 'Epoch {}'.format(i)
#         writer.writerow([row_name] + batch['record_id'])
#         i += 1
#         if i >3:
#             break


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
plt.savefig(out_dir+'/lossPlot_tag'+str(TAG)+'_batchNumber'+str(NUMBER_BATCH)+'_batchSize'+str(BATCH_SIZE)+'_Epoch'+str(EPOCH)+'_Lr'+str(LR)+'.png',dpi=200)
lossTable.to_csv(out_dir+'/lossTable_tag'+str(TAG)+'_batchNumber'+str(NUMBER_BATCH)+'_batchSize'+str(BATCH_SIZE)+'_Epoch'+str(EPOCH)+'_Lr'+str(LR)+'.csv')


# In[34]:


# #def generate_batch(dataset,batch_size):
# batch_size = BATCH_SIZE
# nrow = dataset.shape[0]
# l = list(range(nrow))
# #random.shuffle(l)
# n = batch_size
# x = []
# for i in range(0,len(l),n):
#     each = l[i:i+n]
#     random.shuffle(each)
#     x.append(each)
# print(x)
# #x = [l[i:i + n] for i in range(0, len(l), n)]

# #return(x)


# In[41]:


#batch_gen = batch_generator(hard_neg_val,100)
#len(list(batch_gen))
#export_index('my_gen.csv',batch_gen,3)
#hard_neg_val.shape[0]//100


# In[ ]:


# def loading_data(dataset,batch_index,which_batch,npy_dir = NPY_DIR):
# #data_loader = DataLoader(easy_neg,batch_size = BATCH_SIZE,shuffle = False)
# #    batch_index = generate_batch(dataset,batch_size)
# #print(batch_index[0])
# #    start_time = time.time()
#     data_load = dataset.iloc[batch_index[which_batch]]
# #    load_complete = time.time()
#     record_id = data_load['record_id'].tolist()
#     antigen_index = data_load['Project']+'/'+data_load['id']
#     antigen_index = antigen_index.tolist()
#     #print('loading time: '+str(load_complete - start_time))
# #    print(data_load.head())
#     BCR_b, BCR_w = embedding_BCR(data_load)
# #    bcr_embedding_complete = time.time()
#     #print('BCR embedding time: '+str(bcr_embedding_complete - load_complete))
# #    antigen_dict = [antigens[x] for x in antigen_index]
# #    antigen_embedding_complete = time.time()
#     #print('antigen embedding time: '+str(antigen_embedding_complete - bcr_embedding_complete))
#     return (record_id,BCR_b,BCR_w,antigen_index)


# In[ ]:


# import numpy as np

# def batch_generator(data, batch_size):
#     num_batches = len(data) // batch_size
#     for i in range(num_batches):
#         batch = data[i * batch_size : (i + 1) * batch_size]
#         yield batch
#         del batch

# # Load data
# data = np.load("data.npy")

# # Set batch size
# batch_size = 32

# # Create batch generator
# batch_gen = batch_generator(data, batch_size)

# # Loop over batches
# for batch in batch_gen:
#     # Process batch
#     ...


# In[13]:


# antigen_dict = {}
# antigen_embed_example = extract_antigen(antigen_dict,'Mason/HER2')
# print(antigen_dict.keys())
# antigen_embed_example = extract_antigen(antigen_dict,'Libra-seq/HA')
# print(antigen_dict.keys())
# antigen_embed_example = extract_antigen(antigen_dict,'Mason/HER2')
# print(antigen_dict.keys())


# In[57]:


# BATCH_SIZE = 10
# batch_gen_val = batch_generator(hard_neg_val,BATCH_SIZE)
# antigen_dict = {}
# i = 0
# for batch in batch_gen_val:
#     print('batch:',i)
#     b_dict, w_dict = load_single_batch(batch,antigen_dict)
#     print(antigen_dict.keys())
#     out_b = model_mix(b_dict)
#     out_w = model_mix(w_dict)
#     print('out_b:',out_b[:3],'out_w:',out_w[:3])
# #    val_success += success
#     loss = loss_function(out_b,out_w)
#     print('Loss: ',loss.item())
#     i += 1
#     if i > 3:
#         break


# In[56]:


# start_time = time.time()
# single_antigen = extract_antigen(antigen_dict,antigen_index0[0])
# extract_antigen_time = time.time()
# cb_b = comb_embedding(single_antigen,BCR_b0[0])
# combing_time = time.time()
# print(aalen(single_antigen))
# print('extracting antigen:',extract_antigen_time - start_time)
# print('comb_embedding:',combing_time - extract_antigen_time)
# print('All_through:',combing_time - start_time)


# In[57]:


# start_time = time.time()
# len = aalen(single_antigen)
# single_antigen = extract_antigen(antigen_dict,antigen_index0[0])
# single_antigen_g = torch.from_numpy(single_antigen).to(device)
# extract_antigen_time = time.time()
# BCR_0 = torch.from_numpy(BCR_b0[0]).to(device)
# BCR_t = torch.tile(BCR_0,(1,len,len,1))
# Cb_BCR_antigen = torch.cat((single_antigen_g,BCR_t),dim=3)
# combing_time = time.time()
# print(BCR_t.shape)
# print(Cb_BCR_antigen.shape)
# print('extracting antigen:',extract_antigen_time - start_time)
# print('comb_embedding:',combing_time - extract_antigen_time)
# print('All_through:',combing_time - start_time)


# In[66]:


# start_time = time.time()
# cb_ex = comb_embedding(single_antigen,BCR_b0[0])
# print(time.time() -start_time)
# print(cb_ex.device)


# In[47]:


# print(antigen_dict.keys())
# start_time= time.time()
# antigen1=extract_antigen(antigen_dict,'Shan/Alpha',verbose = True)
# print(antigen1.shape[1])
# antigen1_time = time.time()
# antigen2=extract_antigen(antigen_dict,'Shan/Gamma',verbose =True)
# print(antigen2.shape[1])
# antigen2_time = time.time()
# antigen3=extract_antigen(antigen_dict,'Libra-seq/HA',verbose = True)
# print(antigen3.shape[1])
# antigen3_time = time.time()
# print('Antigen1 time:',antigen1_time - start_time)
# print('Antigen2 time:',antigen2_time -  antigen1_time)
# print('Antigen3 time:',antigen3_time - antigen2_time)


# In[8]:


# def batch_generator(dataset,batch_size,npy_dir = NPY_DIR):
#     batch_index = generate_batch(dataset,batch_size)
#     num_batches = len(batch_index)
#     for i in range(num_batches):
#         data_load = dataset.iloc[batch_index[i]]
#         record_id = data_load['record_id'].tolist()
#         antigen_index = data_load['Project']+'/'+data_load['id']
#         antigen_index = antigen_index.tolist()
#         BCR_b, BCR_w = embedding_BCR(data_load)
#         batch = {'record_id':record_id,'BCR_b':BCR_b,'BCR_w':BCR_w,'antigen_index':antigen_index}
#         yield batch
#         del batch


# In[11]:


# def load_single_batch(batch,antigen_dict,verbose = False,npy_dir = NPY_DIR):
#     record_id = batch['record_id']
#     BCR_b = batch['BCR_b']
#     BCR_w = batch['BCR_w']
#     antigen_index = batch['antigen_index']
#     b_dict = {}
#     w_dict = {}
#     for i in range(len(record_id)):
# #        if verbose:
# #            print('record:',i)
#         start_time = time.time()
#         antigen = extract_antigen(antigen_dict,antigen_index[i],npy_dir,verbose)
#         print('antigen length:',antigen.shape[1])
#         complete_extract_antigen = time.time()
#         print('antigen time:',complete_extract_antigen - start_time)
#         cb_b = comb_embedding(antigen,BCR_b[i,])
#         cb_w = comb_embedding(antigen,BCR_w[i,])
#         print('combing time:',time.time()-complete_extract_antigen)
#         torch_b = torch.tensor(cb_b)
#         torch_w = torch.tensor(cb_w)
#         #print(cb_b.shape,cb_w.shape)
#         b_dict[record_id[i]] = torch_b
#         w_dict[record_id[i]] = torch_w
#     return(b_dict,w_dict)


# In[34]:


# batch_gen = batch_generator(hard_neg_val,25)
# for batch in batch_gen:
#     antigen_index0 = batch['antigen_index']
#     BCR_b0 = batch['BCR_b']
#     print(batch['antigen_index'])
#     break


# In[ ]:





# In[40]:


# print('before start:',torch.cuda.memory_allocated()/(1024**3))
# antigen_dict_in = {}
# antigen_dict_out = {}
# antigen_set = ['46472/COV_MERS', '46472/COV_HKU1', '46472/COV_OC43']
# rotate_antigens(antigen_set,antigen_dict_in,antigen_dict_out)
# print('IN antigens',antigen_dict_in.keys())
# print('OUT antigens:',antigen_dict_out.keys())
# print(torch.cuda.memory_allocated()/(1024**3))
# rotate_antigens(['46472/COV_SARS', '46472/COV_OC43'],antigen_dict_in,antigen_dict_out)
# print('IN antigens',antigen_dict_in.keys())
# print('OUT antigens:',antigen_dict_out.keys())
# print(torch.cuda.memory_allocated()/(1024**3))
# rotate_antigens(['Chappert/B5', 'Makowski/HGFR', 'Makowski/ova', 'Shan/Alpha'],antigen_dict_in,antigen_dict_out)
# print('IN antigens',antigen_dict_in.keys())
# print('OUT antigens:',antigen_dict_out.keys())
# print(torch.cuda.memory_allocated()/(1024**3))
# rotate_antigens(hard_neg_val['antigen_index'].unique(),antigen_dict_in,antigen_dict_out)
# print('IN antigens',antigen_dict_in.keys())
# print('OUT antigens:',antigen_dict_out.keys())
# print(torch.cuda.memory_allocated()/(1024**3))
# #rotate_antigens(hard_neg_train['antigen_index'].unique(),antigen_dict_in,antigen_dict_out)
# #print('IN antigens',antigen_dict_in.keys())
# #print('OUT antigens:',antigen_dict_out.keys())


# In[49]:


# torch.cuda.empty_cache()


# In[25]:


# #hard_neg_val['antigen_index']
# hard_neg_val['antigen_index']=hard_neg_val['Project']+'/'+hard_neg_val['id']
# antigen_index = list(set(hard_neg_val['antigen_index']))
# print(len(antigen_index))
# print(antigen_index)
# print(hard_neg_val['antigen_index'].nunique())


# In[28]:


# hard_neg_train.shape
# hard_neg_val.shape


# In[72]:


# def batch_generator(dataset,batch_size,npy_dir = NPY_DIR):
#     batch_index = generate_batch(dataset,batch_size)
#     num_batches = len(batch_index)
#     for i in range(num_batches):
#         data_load = dataset.iloc[batch_index[i]]
#         record_id = data_load['record_id'].tolist()
#         antigen_index = data_load['Project']+'/'+data_load['id']
#         antigen_index = antigen_index.tolist()
#         BCR_b, BCR_w = embedding_BCR(data_load)
#         batch = {'record_id':record_id,'BCR_b':BCR_b,'BCR_w':BCR_w,'antigen_index':antigen_index}
#         yield batch
#         del batch


# In[83]:


#import torch

# len = 1000
# # Define the shape of the tensor
# shape = (1, len, len, 318)

# # Create a tensor with the specified shape
# tensor = torch.randn(*shape)

# # Predict the size of the tensor
# num_elements = tensor.numel()
# element_size = tensor.element_size()
# tensor_size = num_elements * element_size

# # Print the predicted tensor size
# print(f"Tensor size: {tensor_size / (1024^3)} GB")


# In[153]:


# Get the total memory capacity of the first GPU
# total_memory = torch.cuda.get_device_properties(0).total_memory

# # Get the current GPU memory usage
# allocated_memory = torch.cuda.memory_allocated()
# reserved_memory = torch.cuda.memory_reserved()

# # Calculate the free memory
# free_memory = total_memory - allocated_memory - reserved_memory

# # Print the memory usage
# print(f"Total memory: {total_memory}")
# print(f"Allocated memory: {allocated_memory}")
# print(f"Reserved memory: {reserved_memory}")
# print(f"Free memory: {free_memory}")

# # Release any cached memory that is not in use by PyTorch tensors
# torch.cuda.empty_cache()

# # Get the current GPU memory usage after releasing cached memory
# allocated_memory = torch.cuda.memory_allocated()
# reserved_memory = torch.cuda.memory_reserved()

# # Calculate the free memory after releasing cached memory
# free_memory = total_memory - allocated_memory - reserved_memory

# # Print the memory usage after releasing cached memory
# print(f"Total memory: {total_memory}")
# print(f"Allocated memory: {allocated_memory}")
# print(f"Reserved memory: {reserved_memory}")
# print(f"Free memory: {free_memory}")


# In[34]:


### a test for memory use
# antigen_index_ex = dataset['antigen_index'].unique().tolist()
# antigen_dict = {}
# dict_in = {}
# for single_antigen in antigen_index_ex:
#     print(single_antigen)
#     new_antigen = extract_antigen(antigen_dict,single_antigen)
#     print('before to dict',torch.cuda.memory_allocated()/(1024**3))
#     dict_in[single_antigen]=torch.from_numpy(new_antigen).to(device)
#     print('after to dict',torch.cuda.memory_allocated()/(1024**3))
#     del dict_in[single_antigen]
#     print('after delete the item in dict',dict_in.keys())
#     print(torch.cuda.memory_allocated()/(1024**3))
#     del new_antigen
#     print('after delete new_antigen:',torch.cuda.memory_allocated()/(1024**3))


# In[51]:


# import torch
# import torch.nn as nn

# # Define the model
# class MyModel(nn.Module):
#     def __init__(self, input_size, hidden_size, output_size):
#         super(MyModel, self).__init__()
#         self.linear1 = nn.Linear(input_size, hidden_size)
#         self.linear2 = nn.Linear(hidden_size, output_size)

#     def forward(self, inputs, index_list):
#         x = inputs[index_list[0]]
#         x = self.linear1(x)
#         x = nn.functional.relu(x)
#         x = self.linear2(x)
#         return x

# # Define the input dictionary
# input_dict = {
#     "example_1": torch.randn(20).to(device),
#     "example_2": torch.randn(20).to(device),
#     "example_3": torch.randn(20).to(device)
# }

# # Define the list of indices
# index_list = ["example_2", "example_1", "example_3"]

# # Create an instance of the model
# model = MyModel(input_size=20, hidden_size=50, output_size=2).to(device)

# # Forward pass through the model
# output = model(input_dict, index_list)

# # Print the output tensor
# print(output)


# In[95]:


# trial for model
# class trial_model(nn.Module):
#     def __init__(self):
#         super(trial_model,self).__init__()
#         self.model1 = nn.Sequential(
#             nn.Linear(318,50).to(torch.float64),
#             # in (1,len,len,318)
#             # out (1,len,len.50)
#             nn.ReLU(),
#             nn.Linear(50,20).to(torch.float64),
#             nn.ReLU(),
#             nn.Linear(20,20).to(torch.float64),
#             # out (1,len,len,20)
#             nn.ReLU())
# #        self.sigmoid = nn.Sigmoid()
#     def forward(self, a:dict, index_list, b):
#         k = comb_embedding(index_list[0],b,a)
#         print(k.device,k.dtype)
#         out = self.model1(k)
#         return out
        


# In[96]:


# t_model = trial_model()
# t_model.to(device)
# for name, param in t_model.named_parameters():
#     print(name, param.dtype)


# In[129]:


# CDR3h_dict['CARDREDVYNAYPGSLDLW']


# In[156]:


# #hard_neg_train['antigen_index'][0]
# antigen = antigen_dict_in['46472/COV_MERS']
# print('antigen dtype:',antigen.dtype)
# BCR_b, BCR_w = embedding_BCR(hard_neg_train.loc[:10,:],CDR3h_dict,Vh_dict)
# print('BCR_b dtype:',BCR_b.dtype)
# BCR_b_g = torch.from_numpy(BCR_b).to(device)
# BCR_w_g = torch.from_numpy(BCR_w).to(device)
# print('bcr dtype:',BCR_b_g.dtype)

