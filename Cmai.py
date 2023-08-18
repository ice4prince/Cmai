#!/usr/bin/env python
# coding: utf-8

# In[1]:


import subprocess
import os
import pandas as pd
import argparse
import shutil
import glob
import re
# from Bio import SeqIO #NOTE: proprocess need some prerequiste
# Define the function to run the preprocess.py
# def run_preprocess():
#     subprocess.check_call(['conda', 'activate', 'runBC', '&&', 'python', 'scripts/Check_Proprocess.py']+ args)

# In[ ]:


parser = argparse.ArgumentParser(description='Parameters for the interface script.')
parser.add_argument('--code', type=str, help='the Cmai directory',default = '/project/DPDS/Wang_lab/shared/BCR_antigen/code/Cmai')
parser.add_argument('--input',type = str, help = 'the input files in csv which should include Antigen_id,BCR_Vh,BCR_CDR3h',default = 'data/example/input.csv')
parser.add_argument('--out',type = str, help = 'the directory for output files',default = '/project/DPDS/Wang_lab/shared/BCR_antigen/code/Cmai/data/example/output')
parser.add_argument('--env_path', help='the file saving the directory of the Conda environments- python of runEmbed, python of runBind, and RoseTTAFold in order.',default = 'paras/env_path')
parser.add_argument('--rf_data',type = str, help = 'the database folder for RoseTTAFold',default= '/project/DPDS/Wang_lab/shared/BCR_antigen/data')
parser.add_argument('--fasta',type = str, help = 'The fasta file entering runEbed. When no sequence included in the input, the separate fasta file of antigens is required',default =None)
parser.add_argument('--pre_dir',type = str, help='the directory to save the preprocessed data. If not defiend, same with output directory.',default = None)
parser.add_argument('--npy_dir',type = str, help = 'the npy folder if different with preprocess folder',default = None)
parser.add_argument('--cpu',type = str, help = 'the maximum of cpus for antigen embedding. If not defined, use the value of paras/rf_para.txt',default = 8)
parser.add_argument('--mem',type = str, help = 'the maximum of memory in GB for antigen embedding. If not defined, use the value of paras/rf_para.txt',default = 8)
parser.add_argument('--use_cpu',action = 'store_true', help = 'the option to use cpu or gpu.',default = False)
parser.add_argument('--seed', type=int, help='the seed for the first 100 background BCRs. To use the prepared embeded 100 BCRs, keep the seed to default 1',default = 1)
parser.add_argument('--subsample', type=int, help='the initial sample size of background BCRs. The default is 100',default = 100)
parser.add_argument('--bottomline', type=int, help='the maximum size for subsample of background BCRs, which should no more than 1000000. The default is 10000',default = 10000)
# parser.add_argument('--continuous', action='store_true', help='swtich the mode from binary to continuous, default mode is binary.')
parser.add_argument('--rf_para',action = 'store_true',help = 'use the parameters from paras/rf_para.txt for antigen embedding. Default is False')
parser.add_argument('--gen_msa',action = 'store_true',help = 'only run generating msa and exit. Default is False')
parser.add_argument('--run_rf',action = 'store_true',help = 'skip generating msa and running embedding prediction. Default is False')
parser.add_argument('--skip_preprocess',action = 'store_true',help = 'skip preprocess of antigen_embedding. Default is False')
parser.add_argument('--skip_extract',action = 'store_true',help = 'skip extracting NPY for antigen embedding. Default is False')
parser.add_argument('--runEmbed',action = 'store_true',help = 'only run antigen embedding. Default is False')
parser.add_argument('--runBind',action = 'store_true',help = 'only run binding or comparing. Default is False')
parser.add_argument('--skip_check',action = 'store_true',help = 'skip check and preprocess of input data, only use when it has been done before. Default is False')
parser.add_argument('--species', action='store_true', help='match the species of background BCR to the target BCR. NOTE: the species MUST BE specified and unique in the target BCR input.')
parser.add_argument('--suffix', action='store_true', help='Adding suffix to antigen id. Only use to distinguish same-name antigens. The default is False.')
parser.add_argument('--verbose', action='store_true', help='Enable verbose output, default is False.')
parser.add_argument('--merge', action='store_true', help='Enable merging output to input, default is False.')

args = parser.parse_args()


# In[ ]:
CODE_DIR = args.code
OUT = args.out
if args.pre_dir is None:
    args.pre_dir = OUT

if args.npy_dir is not None:
    NPY_DIR = args.npy_dir
else:
    NPY_DIR = args.pre_dir+'/NPY'
# CONT = args.continuous


#     from Check_Proprocess import MODE,CODE_DIR # import the MODE after running the preprocess.py
#     return MODE,CODE_DIR

# Define the function to activate a conda environment and run a script
    # Note: you must have 'conda' in your PATH for this to work
def run_script_in_conda_env(python_path, script, args, working_directory= os.path.dirname(os.path.abspath(__file__))):
    subprocess.check_call([python_path, script] + args, cwd=working_directory)
# #
# def run_script_in_conda_env(conda_env_path, script, args):
#     subprocess.check_call([conda_env_path, script] + args)
# def run_script_in_conda_env(python_path, script, args, working_directory):
#     subprocess.check_call([python_path, script] + args, cwd=working_directory)

# def run_script_in_conda_env(conda_env, script, args, working_directory=CODE_DIR):
#     subprocess.check_call(['conda', 'run', '-n', conda_env, 'python', script]+ args, cwd=working_directory)

#     subprocess.check_call([python_path, script] + args, cwd=working_directory)
#    subprocess.run(['conda', 'run', '-n', conda_env, 'python', script] + args, cwd=working_directory)
#
# def run_script_in_outer_env(script, args):
#     subprocess.check_call(['python', script] + args)
def run_preprocess(conda_env,args):
    preprocess_args = ['--code', args.code, '--input', args.input, '--pre_dir', args.pre_dir]
    if args.fasta is not None:
        preprocess_args.append('--fasta')
        preprocess_args.append(args.fasta)
    # if CONT:
    #     preprocess_args.append('--continuous')
#        print('rinBind ',bind_args))

    print('Check_Proprocess ',' '.join(preprocess_args))
    run_script_in_conda_env(conda_env,'scripts/Check_Preprocess.py',preprocess_args)

def run_embed(conda_env,args,path_rf_env):
    if args.fasta is None:
        fasta_file = args.pre_dir+'/antigens.fasta'
    else:
        fasta_file = args.fasta
    embed_args = ['--in-fasta', fasta_file,
                 '--env.RF_RUNTIME_BASE',args.out+'/RFoutputs',
                 '--env.RF_DATA_BASE',args.rf_data,
                 '--env.RF-BASE',CODE_DIR+'/scripts/rfold/RoseTTAFold'
                ]
    if args.rf_para:
        rf_para = pd.read_csv('paras/rf_para.txt',sep='\t').set_index('para')['arg'].to_dict()
        # if args.use_cpu is not None and args.use_cpu != rf_para['use_cpu']: rf_para['use_cpu'] =  args.use_cpu
        if args.cpu != rf_para['cpu']: rf_para['cpu'] =  args.cpu
        if args.mem != rf_para['mem']: rf_para['mem'] =  args.mem
    #    add_arg = []
        for key,value in rf_para.items():
            embed_args.append('--'+key)
            embed_args.append(value)
    else:
        args_to_add = ['--runtime.cpu',str(args.cpu),'--runtime.mem',str(args.mem),
                         '--env.CSBLAST-DATA','rfold/RoseTTAFold/csblast-2.2.3/data',
                         '--env.PSIPRED-DATA',path_rf_env+'/share/psipred_4.01/data',
                         '--env.BLASTMAT',path_rf_env+'/share/blast-2.2.26/data',
                         '--exe.hhsuit-path','rfold/deps/hhsuite/bin',
                         '--exe.psipred-path',path_rf_env+'/bin',
                         '--exe.csblast-path','rfold/RoseTTAFold/csblast-2.2.3/bin',
                         '--exe.blast-path',path_rf_env+'/bin'
                         ]
        embed_args.extend(args_to_add)
        # if args.use_cpu is not None:
        #     embed_args.append('--use-cpu')
        #     embed_args.append(args.use_cpu)
        #embed_args.append(add_arg)
    if args.use_cpu:
        embed_args.append('--runtime.use-cpu')
    if args.gen_msa:
        embed_args.append('--gen-msa')
    if args.run_rf:
        embed_args.append('--run-rf')
    if args.skip_preprocess:
        embed_args.append('--skip-preprocess')
    if args.skip_extract:
        embed_args.append('--skip-extract')
    if args.verbose:
        embed_args.append('--verbose')
    if args.suffix:
        embed_args.append('--suffix')
    # print('runEmbed',embed_args)
    print('runEmbed',' '.join(embed_args))
    # subprocess.run(['conda', 'run', '-n', 'runEmbed', 'python', 'runEmbed.py']+ embed_args, cwd='scripts',capture_output=False)
    run_script_in_conda_env(conda_env,'runEmbed.py',embed_args,working_directory='scripts')
#    run_script_in_conda_env('runEmbed', 'runEmbed.py', embed_args, 'scripts')
#    run_script_in_conda_env(env, 'runEmbed.py', embed_args, 'scripts')

# Switch to conda env2
def run_binding(conda_env,args):
    # if mode == 'binary':
    bind_args = ['--code',CODE_DIR,'--input',args.pre_dir,'--out',OUT,'--seed',str(args.seed),'--subsample',str(args.subsample),'--bottomline',str(args.bottomline)]
    if args.npy_dir is not None:
        bind_args.append('--npy_dir')
        bind_args.append(args.npy_dir)
    if args.species:
        bind_args.append('--species')
    if args.verbose:
        bind_args.append('--verbose')
    if args.merge:
        bind_args.append('--merge')
    print('rinBind ',' '.join(bind_args))

    # subprocess.run(['conda', 'run', '-n', 'runBind', 'python', 'scripts/runBind.py'] + bind_args,capture_output=False)
#    subprocess.check_call(['conda', 'run', '-n', 'runBind', 'python', 'scripts/runBind.py']+ bind_args)
#    run_script_in_outer_env('scripts/runBind.py',bind_args)
#    run_script_in_conda_env('runBind','scripts/runBind.py',bind_args)
    run_script_in_conda_env(conda_env,'scripts/runBind.py',bind_args)

    # elif mode == 'continuous':
    #     cont_args = ['--code',CODE_DIR,'--input',args.pre_dir,'--out',OUT,'--seed',str(args.seed),'--subsample',str(args.subsample),'--bottomline',str(args.bottomline)]
    #     if args.species:
    #         cont_args.append('--species')
    #     if args.verbose:
    #         cont_args.append('--verbose')
    #     print('runCompare ',' '.join(cont_args))
    #     run_script_in_outer_env('scripts/runCompare.py',cont_args)

def move_npy(source_dir,des_dir):
    if not os.path.exists(des_dir):
        os.makedirs(des_dir)
    for source_file in glob.glob(source_dir + '/*.pair.npy'):
        des_file = re.sub(r'_\d+', '', source_file.split('/')[-1])
    # Use shutil.move to move the file
        shutil.copy2(source_file, des_dir+'/'+des_file)


if not os.path.exists(OUT):
    # If not, create the directory
    os.makedirs(OUT)
    os.makedirs(OUT+'/RFoutputs')

if not os.path.exists(args.pre_dir):
    # If not, create the directory
    os.makedirs(args.pre_dir)
    os.makedirs(args.pre_dir+'/NPY')
# if CONT:
#     MODE = 'continuous'
#     print('mode switching ...')
# else:
#     MODE = 'binary'
# print('Entering',MODE,'mode now...')


os.chdir(CODE_DIR)
with open(args.env_path, 'r') as file:
    lines = file.readlines()
# Strip any leading/trailing whitespace or newline characters
path_embed = lines[0].strip()
path_bind = lines[1].strip()
path_rf = lines[2].strip()

print(f'path_embed: {path_embed}')
print(f'path_bind: {path_bind}')
# Run preprocess.py and get MODE
if not args.skip_check:
    run_preprocess(path_bind,args)

if args.runEmbed and not args.runBind:
    run_embed(path_embed,args,path_rf)
    move_npy(args.out+'/RFoutputs/pred',NPY_DIR)
    exit(0)
if not args.runEmbed and args.runBind:
    run_binding(path_bind,args)
    exit(0)
# run_embed(args,runEmbed_env_path)
# move_npy(args.out+'/RFoutputs/results/pred',NPY_DIR)
# run_binding(args)
run_embed(path_embed,args,path_rf)
move_npy(args.out+'/RFoutputs/pred',NPY_DIR)
run_binding(path_bind,args)
