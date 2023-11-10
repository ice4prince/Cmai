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
#from config import Conf, config
torch.cuda.empty_cache()

torch.cuda.is_available()

parser = argparse.ArgumentParser(add_help=True)
parser.add_argument('--input','-i',type=str
                    #,required=True
                    ,help='The input file.')
parser.add_argument('--output','-o',type=str
                    #,required=True
                    ,help='The folder to save the results.')
parser.add_argument('--type','-t',type=str,default='sequence',help='which type of input it is, sequence or gene id. The defualt is sequence')
parser.add_argument('--species','-s',type=str,default = 'human', help = 'when input file is gene ids, species must be defined: human or mouse, in which human is defualt.')
#args=parser.parse_args()
args, unknown = parser.parse_known_args()
#parser.add_argument('--batch',type=int,default=1000,help ='The batch size input')
#parser.add_argument('--lr',type=float,default=0.001,help = 'the learning rate' )
#parser.add_argument('--epoch',type=int,required=500,help = 'the maximum of epochs')
#parser.add_argument('--mode','-m',type = 'str',required = True,help = 'somatic or germline')
#parser.add_argument('--geneID','-g',type=str,default=NA, help = 'A file lising the gene ids of the V genes to input.')


# Hyper Parameters
#EPOCH = 500
CHANNEL0 =3*20
CHANNEL1 = 5*20
CHANNEL2 = 50*20
CHANNEL3 = 25*20
CHANNEL4 = 5*20
CHANNEL5 = 15
PADDING = 105
N_TRANS = 1
#PADDING = 105
#SAMPLE_SIZE = 100000
#BATCH_SIZE = 2000
#VAL_SIZE = 5000
#LR = 0.001

device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")
print(device)


def preprocess(V_seq):
    prcentNA=V_seq.count('NA')/len(V_seq)
    if prcentNA > 0:
        print('NA sequence is deleted:',prcentNA)
    V_seq = [seq for seq in V_seq if str(seq) != 'NA']
    for i in range(len(V_seq)):
      #print(dfseq[i])
      if V_seq[i].count('X')>5 or V_seq[i].count('*')>5:
            print('The '+str(i)+'th sequence: '+ V_seq[i]+' has more than 5 stop codons or undefiend nucleotides, and it is to be deleted!')
            V_seq[i] = 'NA'
    if V_seq.count('NA') > 0:
        print(V_seq.count('NA')/len(V_seq)*100,'% Bad sequence to delete.')
    V_seq = [seq for seq in V_seq if str(seq) != 'NA']
    return(V_seq)


def aamapping(peptideSeq,aa_dict,backward=False):#Transform aa seqs to Atchley's factors.
    peptideList = []
    peptideSeq=peptideSeq.upper()
    if len(peptideSeq)>PADDING:
        print('Length: '+str(len(peptideSeq))+' over bound!')
        if backward:
            peptideSeq=peptieSeq[-PADDING:]
        else:
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
def seqMap(dataset,aa_dict,backward=False):#Wrapper of aamapping
    all_seq = []
    for seq in dataset:
        newLine = aamapping(seq,aa_dict,backward=backward)
        all_seq.append(newLine)
#    print('aaMap done!')
    return all_seq


class AutoEncoder(nn.Module):
    def __init__(self):
        super(AutoEncoder, self).__init__()
        self.layer1 = nn.Sequential(
            nn.Conv2d(in_channels=1,out_channels = CHANNEL0,kernel_size=(3,1),stride=(2,1)),
                    # input: BATCH_SIZE*1*99*5
                    # output: BATCH_SIZE*CHANNEL1*(105-3+1)/2*5
            nn.ReLU()
        )
            # nn.BatchNorm2d(30),
        self.layer2 = nn.Sequential(
            nn.MaxPool2d((2,1)),
                    # in: BATCH_SIZE*CHANNEL1*52*5
                    # out: BATCH_SIZE*CHANNEL1*26*5
            nn.Conv2d(in_channels=CHANNEL0,out_channels=CHANNEL1,kernel_size=(3,1),stride=(2,1)),
            nn.ReLU(),
                    # in: BATCH_SIZE*CHANNEL0*26*5
                    # out: BATCH_SIZE*CHANNEL1*12*5
            nn.Conv2d(in_channels=CHANNEL1,out_channels=CHANNEL2,kernel_size=(5,1),stride=(2,1)),
            nn.ReLU(),
                # in: BATCH_SIZE*CHANNEL1*12*5
                # out: BATCH_SIZE*CHANEL2*4*5
        )
        self.layer3 = nn.Sequential(
            nn.Flatten(1,3),
            # out: batch_size * CHANNEL2*12*5
            nn.Linear(CHANNEL2*4*5, CHANNEL4),
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
            nn.Linear(in_features=CHANNEL4, out_features=CHANNEL3*4*5),
            nn.ReLU(),
                # nn.Linear(in_features=int(CHANNEL*N_CONV/2), out_features=CHANNEL*N_CONV),
                # nn.ReLU(),
            nn.Unflatten(1,(CHANNEL3,4,5)), # out dim: CHANNEL2 * 4 * 5
            nn.ConvTranspose2d(in_channels=CHANNEL3,out_channels=CHANNEL2,kernel_size=(5,1),stride =(3,1)
                                ,padding=(1,0)
                               ,output_padding=(1,0)
            ),
            # in CHANNEL2*4*5
            # out : CHANNEL2*13*5 (14=(4 *3(stride)-2) + (5(kernel)-1) -2*1(padding)+1(output_padding))
            nn.ReLU(),
            nn.ConvTranspose2d(in_channels=CHANNEL2,out_channels=CHANNEL1,kernel_size=(3,1),stride = (2,1)
        #                      ,output_padding=(1,0)
        #                      ,padding=(1,0)
            ),
             # out : CHANNEL1*27*5
            nn.ReLU(),
            nn.ConvTranspose2d(in_channels=CHANNEL1,out_channels=CHANNEL0,kernel_size=(3,1),stride = (2,1),
         #                     output_padding=(1,0),
                              padding=(1,0)
            ),
            # out: CHANNEL0*53*5
            nn.ReLU(),
            nn.ConvTranspose2d(in_channels=CHANNEL0,out_channels=1,kernel_size=(3,1),stride = (2,1),
         #                     output_padding=(1,0),
                              padding=(1,0)
            ),
            # out: 1*105*5
    )


    def forward(self, x):

        x =  self.layer1(x)

        x = x.permute(2,0,1,3)

        x = x.reshape(x.shape[0],x.shape[1],CHANNEL0*5)

#        padding_mask=torch.sum(x.permute(1,0,2),dim=2)==0

        x = self.transformer_encoder(src = x
#                                 ,src_key_padding_mask = padding_mask
                                    )
        x = x.reshape(x.shape[0],x.shape[1],CHANNEL0,5)

        x= x.permute(1,2,0,3)

        #out dim: batch_size * CHANNEL1 * PADDING * 5

        x = self.layer2(x)

        x = self.layer3(x)

        mu =  self.vae1(x)

#        print(mu.get_device())
#        print('mu shape:',mu.shape)

        sigma = torch.exp(self.vae2(x))

#        print('sigma shape:',sigma.shape)

        Nsample = self.N.sample((mu.shape[0],mu.shape[1])).to(device)

#        print('Nsample:',Nsample.cuda())

#        Nsample.cuda()


#        print(Nsample,Nsample.get_device())
#         if precise:
#             encoded  = mu
#         else:
        encoded  = mu + sigma*Nsample#.cuda()
        decoded =  self.decoder(encoded)

        return mu, encoded, decoded

def InputV(InFile,type='sequence',species='human'):
    if type == 'sequence':
        V_seq = []
        with open(InFile) as file:
            while (line := file.readline().rstrip()):
                V_seq.append(line)
        return(V_seq)
    elif type == 'gene_id':
        gene_ids = []
        with open(InFile) as file:
            while (line := file.readline().rstrip()):
                gene_ids.append(line)
        if species == 'human':
            df = pd.read_csv(dir+'/V_human.csv',index_col=0)
        elif species == 'mouse':
            df = pd.read_csv(dir+'/V_mouse.csv',index_col=0)
        seq_dict = dict(zip(df.gene_id,df.sequence))
        try:
            V_seq =[seq_dict[item] for item in gene_ids]
        except KeyError as e:
            print('There is no gene named '+str(e)+' in the dictionay or this gene has a bad sequence with too many stop codon or undefined Nucleotides. So it is removed from the input list.')
            gene_ids.remove(e.args[0])
            V_seq =[seq_dict[item] for item in gene_ids]
        return(V_seq,gene_ids)

def embedV(V_seq, precise = False,backward =False):

    #BATCH_SIZE = args.batch
    #LR = args.lr
    #EPOCH = args.epoch
    #MODE = 'somatic'
    #SEED = 100

    dir = os.path.realpath(os.path.dirname(__file__))
    aa_dict_dir=os.path.join(dir,'Atchley_factors.csv')
    checkpoint = torch.load(dir+'/model.pth')

    V_seq = preprocess(V_seq)
    aa_dict_atchley=dict()
    with open(aa_dict_dir,'r') as aa:
        aa_reader=csv.reader(aa)
        next(aa_reader, None)
        for rows in aa_reader:
            aa_name=rows[0]
            aa_factor=rows[1:len(rows)]
    #        aa_dict_atchley[aa_name]=np.asarray(aa_factor,dtype='float')
            aa_dict_atchley[aa_name]=list(aa_factor)

    V_seq_pre = seqMap(V_seq,aa_dict_atchley,backward=backward)
    V_tensor = torch.Tensor(np.array(V_seq_pre,dtype=np.float32))
    V_tensor1= torch.unsqueeze(V_tensor,dim=1)
    V_dataset = TensorDataset(V_tensor1)
    V_toPredict = DataLoader(V_dataset, 1)

    # device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")
    # print(device)
    autoencoder = AutoEncoder()
    autoencoder.load_state_dict(checkpoint)
    autoencoder.to(device)

    inseq = []
    outseq = []
    deseq = []

    with torch.no_grad():
        autoencoder.eval()
        for sequence in V_toPredict:
            input_seq = sequence[0].to(device)
            mu,encoded, decoded = autoencoder(input_seq)
            inseq.append(input_seq.detach().cpu().numpy()[0,0,:,:])
            if precise:
                outseq.append(mu.detach().cpu().numpy()[0,:])
            else:
                outseq.append(encoded.detach().cpu().numpy()[0,:])
#            outseq.append(encoded.detach().cpu().numpy()[0,:])
            deseq.append(decoded.detach().cpu().numpy()[0,0,:,:])

    return(outseq,deseq,inseq)

def OutV(V_seq,outseq,deseq,inseq,Output,gene_ids=[],type='sequence',species='human'):
    if not os.path.exists(Output):
        os.makedirs(Output)
    #    os.makedirs(os.path.join(Output,'pths'))
    if not os.path.exists(os.path.join(Output,'heatmaps_V')):
        os.makedirs(os.path.join(Output,'heatmaps_V'))
    encoded_V=pd.DataFrame(outseq)
    if type == 'sequence':
        encoded_V['sequence'] = V_seq
    elif type == 'gene_id':
        encoded_V['gene_id'] = gene_ids
    ###MARK HERE: When bad sequences were deleted, the respective gene ids ALSO NEED to be deleted.
    encoded_V.to_csv(Output+'/encodedV.csv')

    if type == 'sequence':
        index_pick=random.sample(range(len(V_seq)),min(10,len(V_seq)))
        print(index_pick)
        ###select random 10 to plot out ##
        for i in index_pick:
            gene_pick=V_seq[i]
            # print(gene_pick)
            fig,axs=plt.subplots(2,1,figsize=(50,5))
            axs[0].set_title(str(gene_pick)+': decoded aa seq')
            axs[1].set_title(str(gene_pick)+': input aa seq')
            fig.tight_layout(pad=3.0)
            sns.heatmap(np.transpose(deseq[i]),vmax=5,vmin=-5,square=True,cmap='PiYG',center=0,linewidths=0.1,ax=axs[0])
            sns.heatmap(np.transpose(inseq[i]),vmax=5,vmin=-5,square=True,cmap='PiYG',center=0,linewidths=0.1,ax=axs[1])
            plt.savefig(Output+'/heatmaps_V/heatmap_'+str(i)+'th_V.png')
    elif type == 'gene_id':
        gene_ids = gene_ids[0]
        index_pick=random.sample(range(len(gene_ids)),10)
        print(index_pick)
        ###select random 10 to plot out ###
        for i in index_pick:
            gene_pick=gene_ids[i]
            # print(gene_pick)
            fig,axs=plt.subplots(2,1,figsize=(50,5))
            axs[0].set_title(gene_pick+': decoded aa seq')
            axs[1].set_title(gene_pick+': input aa seq')
            fig.tight_layout(pad=3.0)
            sns.heatmap(np.transpose(deseq[i]),vmax=5,vmin=-5,square=True,cmap='PiYG',center=0,linewidths=0.1,ax=axs[0])
            sns.heatmap(np.transpose(inseq[i]),vmax=5,vmin=-5,square=True,cmap='PiYG',center=0,linewidths=0.1,ax=axs[1])
            plt.savefig(Output+'/heatmaps_V/heatmap_'+str(i)+'th_V.png')

if __name__ == "__main__":
    try:
        V_input,_=InputV(InFile = args.input,type = args.type,species = args.species)
        V_output,V_deseq,V_inseq = embedV(V_input)
        OutV(V_input,V_output,V_deseq,V_inseq,Output = args.output,gene_ids=_,type = args.type,species = args.species)
    except ValueError:
        V_input = InputV(InFile = args.input,type = args.type,species = args.species)
#    print(len(V_input))
        V_output,V_deseq,V_inseq = embedV(V_input)
        OutV(V_input,V_output,V_deseq,V_inseq,Output = args.output,type = args.type,species = args.species)
