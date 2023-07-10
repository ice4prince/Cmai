#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import subprocess
import os
from Bio import SeqIO #NOTE: proprocess need some prerequiste
# Define the function to run the preprocess.py
def run_preprocess():
    subprocess.check_call(['python', 'Check_Proprocess.py']+ args)
    from Check_Proprocess import MODE,CODE_DIR # import the MODE after running the preprocess.py
    return MODE,CODE_DIR

# Define the function to activate a conda environment and run a script
    # Note: you must have 'conda' in your PATH for this to work
def run_script_in_conda_env(conda_env, script, args):
    subprocess.check_call(['conda', 'activate', conda_env, '&&', 'python', script] + args)



# In[ ]:


parser = argparse.ArgumentParser(description='Parameters for the interface script.')
parser.add_argument('--code', type=str, help='the code directory', default='/some/default/path')
parser.add_argument('--input', type=str, help='the input file', default='default_input.csv')
args = parser.parse_args()


# In[ ]:


# Run preprocess.py and get MODE
MODE,CODE_DIR = run_preprocess()
os.chdir(CODE_DIR)

# Activate conda env1 and run script1.py
run_script_in_conda_env('DeLAnO', 'runEmbed.py', ['--code', args.code, '--input', args.input])

# Switch to conda env2
if MODE == 'binary':
    run_script_in_conda_env('torch_test', 'runBind.py')
elif MODE == 'continuous':
    run_script_in_conda_env('torch_test', 'runCompare.py')

