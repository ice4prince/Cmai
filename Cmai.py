#!/usr/bin/env python
# coding: utf-8

# In[1]:


import subprocess
import os
import pandas as pd
import argparse
import shutil
import glob


from Bio import SeqIO #NOTE: proprocess need some prerequiste
# Define the function to run the preprocess.py
# def run_preprocess():
#     subprocess.check_call(['conda', 'activate', 'runBC', '&&', 'python', 'scripts/Check_Proprocess.py']+ args)

# In[ ]:


parser = argparse.ArgumentParser(description='Parameters for the interface script.')
parser.add_argument('--code', type=str, help='the CLAnO directory',default = '/project/DPDS/Wang_lab/shared/BCR_antigen/code/CLAnO')
parser.add_argument('--input',type = str, help = 'the input files in csv which should include Antigen_id,BCR_Vh,BCR_CDR3h',default = 'data/example/binary_example.csv')
parser.add_argument('--out',type = str, help = 'the directory for output files',default = '/project/DPDS/Wang_lab/shared/BCR_antigen/code/CLAnO/data/example/output')
parser.add_argument('--env_path', help='directory of the runEmbed Conda environment',default = '/home2/s205236/.conda/envs/runEmbed/bin/python')
parser.add_argument('--rf_data',type = str, help = 'the database folder for RoseTTAFold',default= '/project/DPDS/Wang_lab/shared/BCR_antigen/data')
parser.add_argument('--fasta',type = str, help = 'The fasta file entering runEbed. When no sequence included in the input, the separate fasta file of antigens is required',default =None)
parser.add_argument('--pre_dir',type = str, help='the directory to save the preprocessed data.',default = 'data/intermediates')
parser.add_argument('--seed', type=int, help='the seed for the first 100 background BCRs. To use the prepared embeded 100 BCRs, keep the seed to default 1',default = 1)
parser.add_argument('--subsample', type=int, help='the initial sample size of background BCRs. The default is 100',default = 100)
parser.add_argument('--bottomline', type=int, help='the maximum size for subsample of background BCRs, which should no more than 1000000. The default is 10000',default = 10000)
parser.add_argument('--continuous', action='store_true', help='swtich the mode from binary to continuous, default mode is binary.')
parser.add_argument('--rf_para',action = 'store_true',help = 'use the parameters from paras/rf_para.txt for antigen embedding. Default is False')
parser.add_argument('--gen_msa',action = 'store_true',help = 'only run generating msa and exit. Default is False')
parser.add_argument('--run_rf',action = 'store_true',help = 'skip generating msa and running RoseTTAFold. Default is False')
parser.add_argument('--skip_preprocess',action = 'store_true',help = 'skip preprocess of antigen_embedding. Default is False')
parser.add_argument('--skip_extract',action = 'store_true',help = 'skip extracting NPY for antigen embedding. Default is False')
parser.add_argument('--runEmbed',action = 'store_true',help = 'only run antigen embedding. Default is False')
parser.add_argument('--runBC',action = 'store_true',help = 'only run binding or comparing. Default is False')
parser.add_argument('--species', action='store_true', help='match the species of background BCR to the target BCR. NOTE: the species MUST BE specified and unique in the target BCR input.')
parser.add_argument('--verbose', action='store_true', help='Enable verbose output, default is False.')
args = parser.parse_args()


# In[ ]:


CODE_DIR = args.code
OUT = args.out
CONT = args.continuous

def run_preprocess(args):
    subprocess.check_call(['python', 'scripts/Check_Preprocess.py']+ args)
#     from Check_Proprocess import MODE,CODE_DIR # import the MODE after running the preprocess.py
#     return MODE,CODE_DIR

# Define the function to activate a conda environment and run a script
    # Note: you must have 'conda' in your PATH for this to work
# def run_script_in_conda_env(conda_env, script, args,working_directory):
#     subprocess.check_call(['conda', 'activate', conda_env, '&&', 'python', script] + args, cwd=working_directory)
# #
# def run_script_in_conda_env(conda_env_path, script, args):
#     subprocess.check_call([conda_env_path, script] + args)
def run_script_in_conda_env(python_path, script, args, working_directory):
    subprocess.check_call([python_path, script] + args, cwd=working_directory)

def run_script_in_outer_env(script, args):
    subprocess.check_call(['python', script] + args)
# In[10]:


if not os.path.exists(OUT):
    # If not, create the directory
    os.makedirs(OUT)
    os.makedirs(OUT+'/RFoutputs')
if CONT:
    MODE = 'continuous'
    print('mode switching ...')
else:
    MODE = 'binary'
print('Entering',MODE,'mode now...')

runEmbed_env_path = args.env_path
#runEmbed_env_path = '/home2/s205236/.conda/envs/runEmbed/bin/python'
#runBC_env_path = '/home2/s205236/.conda/envs/runBC/bin/python'
# In[ ]:


os.chdir(CODE_DIR)
# Run preprocess.py and get MODE
preprocess_args = ['--code', args.code, '--input', args.input, '--pre_dir', args.pre_dir]
if args.fasta is not None:
    preprocess_args.append('--fasta')
    preprocess_args.append('args.fasta')
if CONT:
    preprocess_args.append('--continuous')
# print('Check_Proprocess ',preprocess_args)
run_preprocess(preprocess_args)



# Activate conda runEmbed and run runEmbed.py
def run_embed(args,env):
    if args.fasta is None:
        fasta_file = 'data/intermediates/antigens.fasta'
    else:
        fasta_file = args.fasta
    embed_args = ['--fasta', fasta_file, '--path.out', args.out+'/RFoutputs','--path.data',args.rf_data,'--path.scripts',CODE_DIR+'/scripts','--path.rf','rfscripts/RoseTTAFold','--path.logs','logs']
    if args.rf_para:
        rf_para = pd.read_csv('paras/rf_para.txt',sep='\t').set_index('para')['arg'].to_dict()
    #    add_arg = []
        for key,value in rf_para.items():
            embed_args.append('--'+key)
            embed_args.append(value)
        #embed_args.append(add_arg)
    if args.gen_msa:
        embed_args.append('--gen_msa')
        embed_args.append('True')
    if args.run_rf:
        embed_args.append('--run_rf')
        embed_args.append('True')
    if args.skip_preprocess:
        embed_args.append('--skip_preprocess')
        embed_args.append('True')
    if args.skip_extract:
        embed_args.append('--skip_extract')
        embed_args.append('True')
    print('runEmbed',' '.join(embed_args))
    run_script_in_conda_env(env, 'runEmbed.py', embed_args, 'scripts')

# Switch to conda env2
def run_bc(args,mode):
    if mode == 'binary':
        bind_args = ['--code',CODE_DIR,'--input',args.pre_dir,'--out',OUT,'--seed',str(args.seed),'--subsample',str(args.subsample),'--bottomline',str(args.bottomline)]
        if args.species:
            bind_args.append('--species')
        if args.verbose:
            bind_args.append('--verbose')
        print('rinBind ',' '.join(bind_args))
        run_script_in_outer_env('scripts/runBind.py',bind_args)
    elif mode == 'continuous':
        cont_args = ['--code',CODE_DIR,'--input',args.pre_dir,'--out',OUT,'--seed',str(args.seed),'--subsample',str(args.subsample),'--bottomline',str(args.bottomline)]
        if args.species:
            cont_args.append('--species')
        if args.verbose:
            cont_args.append('--verbose')

        print('runCompare ',' '.join(cont_args))
        run_script_in_outer_env('scripts/runCompare.py',cont_args)

def move_npy(source_dir,des_dir):
    for source_file in glob.glob(source_dir + '/*.pair.npy'):
    # Use shutil.move to move the file
        shutil.copy2(source_file, des_dir)

if args.runEmbed and not args.runBC:
    run_embed(args,runEmbed_env_path)
    exit(0)
if not args.runEmbed and args.runBC:
    move_npy(args.out+'/RFoutputs/results/pred',args.pre_dir+'/NPY')
    run_bc(args,MODE)
    exit(0)
run_embed(args,'runEmbed')
move_npy(args.out+'/RFoutputs/results/pred',args.pre_dir+'/NPY')
run_bc(args,MODE)
