#!/usr/bin/env python
# coding: utf-8

# In[1]:


from torch.utils.data import Dataset


# In[2]:


import os
import sys
import pandas as pd
import numpy as np
import random
from torch.utils.data import Dataset, DataLoader, TensorDataset, SubsetRandomSampler
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


# In[24]:
# Create an ArgumentParser object
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
NPY_DIR ='/project/DPDS/Wang_lab/shared/BCR_antigen/data/rfscript/Antigen_embed/CLEANED/'

# BATCH_SIZE = 1
# NUMBER_BATCH = 20
# LR = 0.01
# EPOCH = 50
# CUTOFF = 1800
# TAG = 1
# T = 0.01
# LAMBDA = 0
# w = 1
# LIMIT = 10
# NPY_DIR ='/project/DPDS/Wang_lab/shared/BCR_antigen/data/rfscript/Antigen_embed/CLEANED'
# VERBOSE = True


# In[36]:


device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")
print(device)


# In[4]:


#os.getcwd()
os.chdir('/project/DPDS/Wang_lab/shared/BCR_antigen/code/DeLAnO/')
#mp.set_sharing_strategy('file_system')
#mp.set_start_method('spawn')
from wrapV.Vwrap import embedV
from wrapCDR3.CDR3wrap import embedCDR3


# In[37]:


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
#    df = df.sort_values('antigen_index',)
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


# In[38]:


def get_my_data_list(dataset):
    selected_cols = dataset.columns[[2,3,4,5,15,16]].tolist()
    ds_to_dict = dataset[selected_cols]
#ds_to_dict = dataset[selected_cols].set_index('record_id')
    my_data_list = ds_to_dict.to_dict(orient='records')
    return my_data_list


# In[47]:


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

    def __getitem__(self, idx):
        your_dict = self.your_data_list[idx]
        antigen_key = your_dict['antigen_index']
        betterCDR_key = your_dict['BetterBCR_CDR3h']
        worseCDR_key = your_dict['WorseBCR_CDR3h']
        betterV_key = your_dict['BetterBCR_Vh']
        worseV_key = your_dict['WorseBCR_Vh']
#        aalen_key = your_dict['aalens']

        antigen_feat = self.__extract_antigen(antigen_key)
        better_feat = self.__embedding_BCR(betterCDR_key,betterV_key)
        worse_feat = self.__embedding_BCR(worseCDR_key,worseV_key)
        better_pair = self.__comb_embed(antigen_feat,better_feat)
        worse_pair = self.__comb_embed(antigen_feat,worse_feat)
        better_out = better_pair.squeeze(dim=0)
        worse_out = worse_pair.squeeze(dim=0)

        return better_out.to(device), worse_out.to(device)#, aalen_key #IF RUN COMBINING_EMBED in cpu
    #IF RUN COMBINING_EMBED in GPU:
#        return better_pair, worse_pair
#     def __comb_embed_gpu(single_antigen,BCR_feat):
#         len = aalen(single_antigen)
#         single_antigen_g = torch.from_numpy(single_antigen).cuda()
#         single_BCR_g = torch.from_numpy(BCR_feat).cuda()
#         BCR_t = torch.tile(single_BCR_g,(1,len,len,1))
#         pair_feat_g = torch.cat((single_antigen_g,BCR_t),dim=3)
#         return pair_feat_g
    def __comb_embed(self,single_antigen,BCR_feat):
        '''process your data'''
    #    singe_antigen = antigen_dict[antigen_name]
        lengthen = aalen(single_antigen)
#         if verbose ==True:
#             if lengthen > 500:
#                 print(Fore.RED + 'length of antigen over 500: '+str(lengthen))
#                 print(Style.RESET_ALL)
        BCR_tile = np.tile(BCR_feat,(1,lengthen,lengthen,1))
        pair_np = np.concatenate((single_antigen,BCR_tile),axis=3)
        pair_feat = torch.from_numpy(pair_np)
        return(pair_feat)

    def __extract_antigen(self,antigen_name):
        if antigen_name in self.antigen_dict:
            single_antigen = self.antigen_dict[antigen_name]
        else:
            try:
                single_antigen = np.load(str(self.antigen_fpath_dict+'/'+antigen_name+'.pair.npy'))/w
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
    def __embedding_BCR(self,cdr3_seq,v_seq):
        if cdr3_seq not in self.cdr3_dict:
            df1 = pd.DataFrame([cdr3_seq])
            cdr3_feat,*_ = embedCDR3(df1[0])
            self.cdr3_dict[cdr3_seq]=cdr3_feat
        else:
            cdr3_feat = self.cdr3_dict[cdr3_seq]
        if v_seq not in self.v_dict:
            df2 = pd.DataFrame([v_seq])
            v_feat,*_ = embedV(df2[0])
            self.v_dict[v_seq]=v_feat
        else:
            v_feat = self.v_dict[v_seq]
        bcr_feat = np.concatenate((cdr3_feat,v_feat))
        return bcr_feat

    def __len__(self):
        return len(self.your_data_list)


# In[8]:


def my_collate(batch):
    # Get the maximum length of the tensors in the batch
    max_len = max([x[2] for x in batch])

    # Pad each tensor with zeros to the maximum length
    #for x in batch:
    #    print(x[0].shape,x[1].shape,x[2])
    t1_pad = [F.pad(x[0], [0,0,0,max_len - x[2],0,max_len - x[2]]) for x in batch]
    t2_pad = [F.pad(x[1], [0,0,0, max_len - x[2], 0, max_len - x[2]]) for x in batch]

    # Stack the padded tensors and lengths
    t1_padded = torch.stack(t1_pad, dim=0)
    t2_padded = torch.stack(t2_pad, dim=0)
    lengths = torch.tensor([x[2] for x in batch])

    return t1_padded, t2_padded, lengths


# In[40]:


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


# In[56]:


def train_epoch(dataloader,model,loss_fn,optimizer,device,verbose = False):
    train_loss = 0.0
    i = 0
    model.train()
    for batch in dataloader:
        if verbose:
            print(' Batch:',i)
            before_enter_model = time.time()
        better, worse = batch
        out_b = model(better)
        out_w = model(worse)
        loss = loss_fn(out_b,out_w)
        if math.isinf(loss) or math.isnan(loss):
            prob_bw = [out_b,out_w]
            print('ERROR: The loss is INF or NaN! '+str(prob_bw))
            break
        optimizer.zero_grad()
        loss.backward()
        optimizer.step()
        train_loss += loss.item()
        if verbose:
            print(Fore.BLUE + 'out_b: '+str(out_b[:3]))
            print(Fore.BLUE + 'out_w: '+str(out_w[:3])+Style.RESET_ALL)
            modeling_compelte = time.time()
            print('modeling time:'+str(modeling_compelte - before_enter_model))
            print('loss of the batch:',loss.item())
        i += 1
#        if i > LIMIT:
#            break
    return train_loss/i


def val_epoch(dataloader,model,loss_fn,verbose = False):
    val_loss = 0.0
    val_success = 0
    i = 0
    model.eval()
    for batch in dataloader:
        if verbose:
            print(' Batch:',i)
            before_enter_model = time.time()
        better, worse = batch
        out_b = model(better)
        out_w = model(worse)
        loss = loss_fn(out_b,out_w)
        val_loss += loss.item()
        if BATCH_SIZE == 1:
            success = int(torch.gt(out_w,out_b))
        else:
            success = torch.sum(torch.gt(out_w,out_b)).item()/out_b.shape[0]
        val_success += success

        if verbose:
            print(Fore.BLUE + 'out_b: '+str(out_b[:3]))
            print(Fore.BLUE + 'out_w: '+str(out_w[:3])+Style.RESET_ALL)
            modeling_compelte = time.time()
            print('modeling time:'+str(modeling_compelte - before_enter_model))
            print('loss of the batch:',loss.item())
        i += 1
#         if i > number_of_batches:
#             break
    return val_loss/i,val_success/i


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


# In[42]:


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


# In[10]:


hard_file_train = pd.read_csv('/project/DPDS/Wang_lab/shared/BCR_antigen/data/hard_negative_train_downsampled.csv')
hard_file_val = pd.read_csv('/project/DPDS/Wang_lab/shared/BCR_antigen/data/hard_negative_validation_downsampled.csv')


# In[11]:


hard_neg_train = preprocess(hard_file_train,CUTOFF)
hard_neg_val = preprocess(hard_file_val, CUTOFF)
print('Train Data: ')
check_space(hard_neg_train)
print('Validation Data: ')
check_space(hard_neg_val)


# In[12]:


hard_neg_df = pd.concat([hard_neg_train,hard_neg_val],axis = 0)
Vh_dict = build_BCR_dict(hard_neg_df,'Vh')
CDR3h_dict = build_BCR_dict(hard_neg_df,'CDR3h')


# In[13]:


my_train_list = get_my_data_list(hard_neg_train)
my_val_list = get_my_data_list(hard_neg_val)


# In[48]:


dataset_train = SampleDataset(NPY_DIR,CDR3h_dict,Vh_dict,my_train_list)
dataset_val = SampleDataset(NPY_DIR,CDR3h_dict,Vh_dict,my_val_list)


# In[49]:


dl_train = DataLoader(dataset_train, 1, shuffle=True, collate_fn=None)#num_workers=BATCH_SIZE,
dl_val = DataLoader(dataset_val,1,shuffle=True, collate_fn=my_collate)


# In[43]:


model_mix = mix_model()
#model_mix.to(device)
checkpoint = torch.load('/project/DPDS/Wang_lab/shared/BCR_antigen/data/output/2023-02-18_12-33-25_tag2/model_comb/Batch50BatchNumber_600Lr0.01Epoch54tag2_easy_neg.pth')
model_mix.load_state_dict(checkpoint)
model_mix.to(device)
#for name, param in model_mix.named_parameters():
#    print(name, param.dtype)
optimizer = torch.optim.Adam(model_mix.parameters(),lr=LR)
scheduler = ReduceLROnPlateau(optimizer,mode = 'min',factor=0.1,patience =10,verbose=False,threshold_mode='rel',threshold=1e-5)


# In[55]:


#def main():
    # Set the start method to 'spawn'
#mp.set_start_method('spawn')
# i = 0
# for record in dl_train:
#     start = time.time()
# #    print(len(record))
# #    print(record[0].shape,record[1].shape)
#     better, worse = record
#     out_b = model_mix(better)
#     out_w = model_mix(worse)
#     print('b shape and w shaple:',out_b.shape,out_w.shape)
#     loss_ex = loss_function(out_b,out_w)
#     print('loss:',loss_ex.item())
#     print('gt:',int(torch.gt(out_w,out_b)))
#     success = torch.sum(torch.gt(out_w,out_b)).item()/out_b.shape[0]
#     print('success:',success)
#     print('the',i,'-th batch time:',time.time()-start)
#     i += 1
#     if i > 2:
#         break


# In[ ]:


#NUMBER_BATCH = 20
print('batch size:',BATCH_SIZE,'EPOCH:', EPOCH,'Learning rate:',LR)
start_time = time.time()
print('Initialling ...')
#antigen_dict = {}
initial_train_loss = val_epoch(dl_train,model_mix,loss_function,verbose = False)
initial_val_loss,initial_val_accuracy = val_epoch(dl_val,model_mix,loss_function,verbose = False)
print(time.time()-start_time)


# In[20]:


train_LOSS = [initial_train_loss]
val_LOSS = [initial_val_loss]
val_ACC = [initial_val_accuracy]


# In[ ]:


for epoch in range(EPOCH):
    print('epoch:',str(epoch))
    start_epoch_time = time.time()
    train_loss = train_epoch(dl_train,model_mix,loss_function,verbose = False)
    val_loss,val_accuracy = val_epoch(dl_val,model_mix,loss_function,verbose = False)
    print('Epoch '+str(epoch)+' train loss: '+str(train_loss)+' val loss: '+str(val_loss)+' val accuracy: '+str(val_accuracy))
    scheduler.step(val_loss) ##should I try scheduler.step(-val_accuracy)?
    train_LOSS.append(train_loss)
    val_LOSS.append(val_loss)
    val_ACC.append(val_accuracy)
    state_dict=model_mix.state_dict()
    torch.save(state_dict,model_dir+"/Batch"+str(BATCH_SIZE)+"Lr"+str(LR)+"Epoch"+str(epoch)+"tag"+str(TAG)+"_hard_neg.pth")
    if math.isnan(train_loss):
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
