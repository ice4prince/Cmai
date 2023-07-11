import torch
import torchvision
from torch.utils.data import random_split
import numpy as np
import pandas as pd
import csv
import os
from collections import Counter
import torch.nn as nn
import math
import random
import matplotlib.pyplot as plt
from scipy import stats
from sklearn.model_selection import RepeatedKFold
from torch.utils.data import DataLoader, TensorDataset, SubsetRandomSampler
import seaborn as sns
import matplotlib
import matplotlib.colors as mc
import argparse
#torch.cuda.empty_cache()

torch.cuda.is_available()

parser = argparse.ArgumentParser(add_help=True)
parser.add_argument('--input','-i',type=str,required=True,help='The input file.')
parser.add_argument('--output','-o',type=str,required=True,help='The folder to save the results.')
#parser.add_argument('--type','-t',type=str,default='sequence',help='which type of input it is, sequence or gene id. The defualt is sequence')
#parser.add_argument('--species','-s',type=str,default = 'human', help = 'when input file is gene ids, species must be defined: human or mouse, in which human is defualt.')
args=parser.parse_args()
#parser.add_argument('--batch',type=int,default=1000,help ='The batch size input')
#parser.add_argument('--lr',type=float,default=0.001,help = 'the learning rate' )
#parser.add_argument('--epoch',type=int,required=500,help = 'the maximum of epochs')
#parser.add_argument('--mode','-m',type = 'str',required = True,help = 'somatic or germline')
#parser.add_argument('--geneID','-g',type=str,default=NA, help = 'A file lising the gene ids of the V genes to input.')

InFile = args.input
Output = args.output
#type = args.type
#species = args.species
#BATCH_SIZE = args.batch
#LR = args.lr
#EPOCH = args.epoch
#MODE = 'somatic'
#SEED = 100

if not os.path.exists(Output):
    os.makedirs(Output)
#    os.makedirs(os.path.join(Output,'pths'))
    os.makedirs(os.path.join(Output,'heatmaps'))

dir = os.path.realpath(os.path.dirname(__file__))
aa_dict_dir=os.path.join(dir,'Atchley_factors.csv')
checkpoint = torch.load(dir+'/model.pth')
# Hyper Parameters
#EPOCH = 500

#SAMPLE_SIZE = 100000
#BATCH_SIZE = 2000
#VAL_SIZE = 5000
#LR = 0.001
CHANNEL0 =3*5
CHANNEL1 = 5*5
CHANNEL2 = 50*5
CHANNEL3 = 35*5
CHANNEL4 = 25*5
CHANNEL5 = 15
PADDING = 45
N_TRANS = 2

def preprocess(df):
    prcentNA=df['aminoAcid'].tolist().count('NA')/len(df['aminoAcid'].tolist())
    if prcentNA != 0:
        print('NA sequence is deleted:',prcentNA)
    df=df.dropna(subset=['aminoAcid'])
    dfseq=df['aminoAcid'].tolist()
    for i in range(len(dfseq)):
      #print(dfseq[i])
      if dfseq[i].count('X')>5 or dfseq[i].count('*')>5:
        dfseq[i] = 'NA'
    if dfseq.count('NA') !=0:
        print('Bad sequence to delete:',dfseq.count('NA')/len(dfseq))
    df['aminoAcid']=dfseq
    df=df.dropna(subset=['aminoAcid'])
    return(df)
def aamapping(peptideSeq,aa_dict):#Transform aa seqs to Atchley's factors.
    peptideList = []
    if len(peptideSeq)>PADDING:
        print('Length: '+str(len(peptideSeq))+' over bound!')
        peptideSeq=peptideSeq[0:PADDING]
    for aa_single in peptideSeq:
        try:
            peptideList.append(aa_dict[aa_single])
        except KeyError:
            print('Not proper aaSeqs: '+peptideSeq)
            peptideList.append([0]*5)
    for i in range(0, PADDING - len(peptideSeq)):
        peptideList.append([0]*5)
    return peptideList
def seqMap(dataset,aa_dict):#Wrapper of aamapping
    all_seq = []
    for seq in dataset:
        newLine = aamapping(seq,aa_dict)
        all_seq.append(newLine)
    print('aaMap done!')
    return all_seq


# In[84]:


class AutoEncoder(nn.Module):
    def __init__(self):
        super(AutoEncoder, self).__init__()
        self.layer1 =nn.Sequential(
            nn.Conv2d(in_channels=1,out_channels = CHANNEL0,kernel_size=(3,1),stride=(2,1)),
        # input: BATCH_SIZE*1*45*5
        # output: BATCH_SIZE*CHANNEL0*(45-3+1)/2*1
            nn.ReLU(),
            )
            # nn.BatchNorm2d(30),
        self.layer2 =nn.Sequential(
            nn.MaxPool2d((2,1)),
                    # in: BATCH_SIZE*CHANNEL0*22*5
                    # out: BATCH_SIZE*CHANNEL0*11*5
            nn.Conv2d(in_channels=CHANNEL0,out_channels=CHANNEL1,kernel_size=(2,1),stride=(2,1)),
            nn.ReLU(),
                    # in: BATCH_SIZE*CHANNEL0*11*5
                    # out: BATCH_SIZE*CHANNEL1*5*5
            nn.Conv2d(in_channels= CHANNEL1, out_channels=CHANNEL2,kernel_size=(3,1)),
                    # input dim: batch_size * CHANNEL1 * 5 * 5,
                    # output dim:batch_size * CHANNEL2 * 3 * 5
            nn.ReLU(),
                    # nn.BatchNorm2d(30),
        #    nn.MaxPool2d((2,1)),
                    # output dim: batch_size * CHANNEL * 8 * 5
        #    nn.Conv2d(CHANNEL2, CHANNEL3,(3,1)),
             # input dim: batch_size * CHANNEL2 * 8 * 5 ,
            # output dim:batch_size * CHANNEL3 * 6 * 5
        #    nn.ReLU(),
                    # nn.BatchNorm2d(30),
        #    nn.MaxPool2d((2,1)),
                    # output dim: batch_size* CHANNEL3 * 3 *5
            )
        self.layer3 = nn.Sequential(
            nn.Flatten(),
            # out: batch_size * CHANNEL2*3*5
            nn.Linear(CHANNEL2*3*5, CHANNEL4),
            nn.ReLU(),
            nn.Linear(CHANNEL4, CHANNEL5)
            )

        self.vae1 = nn.Linear(CHANNEL5, CHANNEL5)

        self.vae2 = nn.Linear(CHANNEL5, CHANNEL5)



        self.N = torch.distributions.Normal(0, 1)



        #self.flatten_trans = nn.Flatten(start_dim=1, end_dim=2), #input dim: PADDING * 1 * batch_size * 5, output dim: PADDING * batch_size * 5

        #transformer required dim: seq * batch * feature

        self.transformer_encoder = nn.TransformerEncoder(nn.TransformerEncoderLayer(d_model=CHANNEL0*5, nhead=5), num_layers=N_TRANS) # output dim: PADDING * batch_size *5



        #self.unflatten_trans = nn.Unflatten(0, (PADDING,1))# output dim: PADDING * 1 * batch_size * 5


        # self.linear1 = nn.Linear(in_features=EMBEDDING*N_CONV, out_features=BOTTLENECK)

        self.decoder = nn.Sequential(
            nn.Linear(in_features=CHANNEL5, out_features=CHANNEL4),
            nn.ReLU(),
            nn.Linear(in_features=CHANNEL4, out_features=CHANNEL2*3*5),
            nn.ReLU(),
            nn.Unflatten(1,(CHANNEL2,3,5)), # out dim: CHANNEL2 * 3 * 5
            nn.ConvTranspose2d(in_channels=CHANNEL2,out_channels=CHANNEL1,kernel_size=(6,1),stride =(3,1)
                               ,padding=(1,0)
                               ,output_padding=(1,0)
                              ),
            # in CHANNEL1*3*5
            # out : CHANNEL0*(3*3-1)+(6-1)-2*1+1*1*5
            nn.ReLU(),
            nn.ConvTranspose2d(in_channels=CHANNEL1,out_channels=CHANNEL0,kernel_size=(5,1),stride = (2,1)
                              ,output_padding=(1,0)
                              ,padding=(2,0)
                              ),
            # in CHANNEL1*11*5
            # out : CHANNEL0*22*5
            nn.ReLU(),
            nn.ConvTranspose2d(in_channels=CHANNEL0,out_channels=1,kernel_size=(4,1),stride = (2,1),
                              output_padding=(1,0),
                              padding=(1,0)
                              ),
        #    nn.ReLU(),
            # in: CHANNEL0*22*5
            # out: 1*45*5
        #    nn.ConvTranspose2d(in_channels=CHANNEL0,out_channels=1,kernel_size=(5,1),stride = (2,1)
        #                      ,output_padding=(1,0)
        #                      ,padding=(1,0)
        #                      ),
            # out: BATCH_SIZE*1*242*5
        #    nn.ReLU(),
            )



    def forward(self, x):

        x =  self.layer1(x)

        x = x.permute(2,0,1,3)
#        print(x.shape)
        x = x.reshape(x.shape[0],x.shape[1],CHANNEL0*5)

#    padding_mask=torch.sum(x.permute(1,0,2),dim=2)==0

        x = self.transformer_encoder(src = x
#                                 ,src_key_padding_mask = padding_mask
                                    )
        x = x.reshape(x.shape[0],x.shape[1],CHANNEL0,5)

        x= x.permute(1,2,0,3)
#    print('pro-trans:',x)
        #print('pro-trans.shape:', x.shape
  #,padding_mask,padding_mask.shape
#)

        #out dim: batch_size * CHANNEL1 * PADDING * 5

        x = self.layer2(x)

        x = self.layer3(x)

        mu =  self.vae1(x)

#        print(mu.get_device())
#        print('mu shape:',mu.shape)

        sigma = torch.exp(self.vae2(x))

#        print('sigma shape:',sigma.shape)

        Nsample = self.N.sample((mu.shape[0],mu.shape[1]))

#        print('Nsample:',Nsample.cuda())

        Nsample.cuda()


#        print(Nsample,Nsample.get_device())

#        if self.training:
#            encoded  = mu + sigma*Nsample.cuda()
#        else:
        encoded  = mu + sigma*Nsample.cuda()

        decoded =  self.decoder(encoded)

        return mu, encoded, decoded


CDR3 = []
with open(InFile) as file:
    while (line := file.readline().rstrip()):
        CDR3.append(line)

aa_dict_atchley=dict()
with open(aa_dict_dir,'r') as aa:
    aa_reader=csv.reader(aa)
    next(aa_reader, None)
    for rows in aa_reader:
        aa_name=rows[0]
        aa_factor=rows[1:len(rows)]
        aa_dict_atchley[aa_name]=list(aa_factor)


#random.seed(100)
#CDR3_sample = random.sample(CDR3,SAMPLE_SIZE)
#CDR3_train_pre = seqMap(CDR3_sample,aa_dict_atchley)
CDR3_train_pre = seqMap(CDR3,aa_dict_atchley)
#print(V_seq_train_pre)
CDR3_train_tensor = torch.Tensor(np.array(CDR3_train_pre,dtype=np.float32))
CDR3_train_tensor1= torch.unsqueeze(CDR3_train_tensor,dim=1)
CDR3_train_dataset = TensorDataset(CDR3_train_tensor1)
CDR3_toPredict = DataLoader(CDR3_train_dataset, 1)

#V_toPredict = DataLoader(V_dataset, 1)

autoencoder = AutoEncoder()
autoencoder.load_state_dict(checkpoint)
device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")
print(device)
autoencoder.to(device)

inseq = []
outseq = []
deseq = []
l2List =[]
with torch.no_grad():
    autoencoder.eval()
    for sequence in CDR3_toPredict:
        input_seq = sequence[0].to(device)
        mu, encoded, decoded = autoencoder(input_seq)
        inseq.append(input_seq.detach().cpu().numpy()[0,0,:,:])
        outseq.append(encoded.detach().cpu().numpy()[0,:])
        deseq.append(decoded.detach().cpu().numpy()[0,0,:,:])


encoded_V=pd.DataFrame(outseq)
encoded_V.to_csv(Output+'/encoded_CDR3.csv')


# In[24]:


index_pick=random.sample(range(len(CDR3)),10)
print(index_pick)
###select random 10 to plot out ###
for i in index_pick:
    gene_pick=i
#    print(gene_pick)
    fig,axs=plt.subplots(2,1,figsize=(50,5))
#fig,ax2=plt.subplots(figsize=(30,3))
    axs[0].set_title(str(gene_pick)+'-th gene: decoded aa seq')
    axs[1].set_title(str(gene_pick)+'-th gene: input aa seq')
    sns.heatmap(np.transpose(deseq[i]),vmax=5,vmin=-5,square=True,cmap='PiYG',center=0,linewidths=0.1,ax=axs[0])
#plt.show()
    sns.heatmap(np.transpose(inseq[i]),vmax=5,vmin=-5,square=True,cmap='PiYG',center=0,linewidths=0.1,ax=axs[1])
    plt.savefig(Output+'/heatmaps/heatmap_'+str(i)+'_CDR3.png')
