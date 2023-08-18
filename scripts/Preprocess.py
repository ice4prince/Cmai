
#from optparse import OptionParser
import os
import sys

current = os.path.dirname(os.path.realpath(__file__))
parent = os.path.dirname(current)
sys.path.append(parent)

from Bio import SeqIO
from colorama import Fore, Style
from rfold.rconfig import Conf, config

#parser = OptionParser()
#parser.add_option("-f", "--fasta",
#                  action="store", type="string", dest="fasta")
#parser.add_option("-o","--output",
#                    action = "store",type = "string",dest = "output")
#parser.add_option("-q", "--quiet",
#                  action="store_false", dest="verbose", default=True,
#                  help="don't print get-though messages.")

#(options, args) = parser.parse_args()
# In[3]:

badLetters= ['b', 'j', 'o', 'u', 'x', 'z','*']
good_sequences = []
remove_list = []
def check_fasta(conf: Conf):
    INfasta = conf.in_fasta
    output = conf.env.RF_RUNTIME_BASE
    verbose = conf.verbose
#    print('fasta: %s' %INfasta)
    for fasta in SeqIO.parse(open(INfasta),'fasta'):
        if len([x for x in badLetters if x in fasta.seq.lower()])>0:
            print(Fore.RED +'Antigen %s contains special or unknown amino acid(s) and has been removed!' % fasta.id)
            remove_list.append(fasta.id)
        else:
            if verbose:
                print(Fore.GREEN + 'Antigen %s get through!' %fasta.id)
            good_sequences.append(fasta)
    print(Style.RESET_ALL)
    if '/' in INfasta:
        infasta = INfasta.rsplit('/', 1)[1]
    else:
        infasta = INfasta
    if not os.path.exists(output):
        os.mkdir(output)
        print('Directory %s created! ' %output)
    outputP = os.path.join(output, 'preprocess')
    if not os.path.exists(outputP):
        os.mkdir(outputP)
    outfasta = outputP + '/' + infasta
#    print(good_sequences)
#    print('infasta: %s' %infasta)
#    print('outfasta: %s' %outfasta)
    SeqIO.write(good_sequences, outfasta, "fasta")
    listFile = outputP+'/remove_list.txt'
    with open(listFile, 'w') as fp:
        for item in remove_list:
        # write each item on a new line
            fp.write("%s\n" % item)
    print('Preprocess is completed!')

if __name__ == "__main__":
    check_fasta(config)
