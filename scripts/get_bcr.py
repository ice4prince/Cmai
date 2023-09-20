import pandas as pd
import argparse
import subprocess
import os

parser = argparse.ArgumentParser(description='Parameters for the getting bcr embeddings.')

parser.add_argument('--input',type = str, help = 'the input folder for the preprocessed input',default = '/project/DPDS/Wang_lab/shared/BCR_antigen/code/Cmai/data/example/output')
parser.add_argument('--out',type = str, help = 'the directory for output files',default = '/project/DPDS/Wang_lab/shared/BCR_antigen/code/Cmai/data/example/output')
# parser.add_argument('--env_path',type = str, help = 'the runBind environment path',default = '/home2/s205236/.conda/envs/runBind/bin/python')
parser.add_argument('--verbose',action = 'store_true',help = 'export full outputs of bcr embedding including the heatmaps. Default is False')
args = parser.parse_args()

INPUT_DIR = args.input
OUT_DIR = args.out
INPUT = INPUT_DIR+'/processed_input.csv'
# path_bind = args.env_path
input_file= pd.read_csv(INPUT)

V_seq = input_file['BCR_Vh'].tolist()
CDR3_seq = input_file['BCR_CDR3h'].tolist()

if not args.verbose:
    from wrapV.Vwrap import embedV
    from wrapCDR3.CDR3wrap import embedCDR3
    V_output,*_= embedV(V_seq)
    CDR3_output,*_= embedCDR3(CDR3_seq)
    encoded_V=pd.DataFrame(V_output)
    encoded_V['sequence'] = V_seq
    encoded_CDR3=pd.DataFrame(CDR3_output)
    encoded_CDR3['sequence'] = CDR3_seq
    encoded_V.to_csv(OUT_DIR+'/encoded_V.csv')
    encoded_CDR3.to_csv(OUT_DIR+'/encoded_CDR3.csv')
else:
    def run_script_outer_env(script, args, working_directory= os.path.dirname(os.path.abspath(__file__))):
        # subprocess.check_call([python_path, script] + args, cwd=working_directory)
        subprocess.check_call(['python', script] + args, cwd=working_directory)
    with open(OUT_DIR+'/v_seq.txt', 'w') as file:
        file.write('\n'.join(V_seq))
    with open(OUT_DIR+'/cdr3_seq.txt', 'w') as file:
        file.write('\n'.join(CDR3_seq))
    v_args=['--input', OUT_DIR+'/v_seq.txt','--output',OUT_DIR]
    cdr3_args = ['--input', OUT_DIR+'/cdr3_seq.txt','--output',OUT_DIR]
    # Strip any leading/trailing whitespace or newline characters
    run_script_outer_env('wrapV/Vwrap.py',v_args)
    run_script_outer_env('wrapCDR3/CDR3wrap.py',cdr3_args)
