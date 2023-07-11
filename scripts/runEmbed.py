import os
from glob import iglob

from config import Conf, config
from scripts.Preprocess import preprocess
from scripts.NPZtoPair import exPair
from rfscripts.rfpd import run_msa
from rfscripts.rfpd import run_rfpd
#from rfscripts.config import Conf, config
#parser = OptionParser()
#parser.add_option("-f", "--fasta",
#                  action="store", type="string", dest="fasta")
#parser.add_option("-o","--output",
#                    action = "store",type = "string",dest = "out")
#parser.add_option("-d", "--data",
#                  action="store", type="string", dest="data")
#parser.add_option("-s","--script",
#                    action = "store",type = "string",dest = "scripts")
#parser.add_option("-q", "--quiet",
#                  action="store_false", dest="verbose", default=True,
#                  help="don't print get-though messages.")


def Pre(conf: Conf) -> None:
    print('Preprocessing...')
    preprocess(conf)

def runEmbed(conf: Conf) -> None:
    print('Embedding...')
    if not os.path.exists(conf.path.out):
        os.mkdir(conf.path.out)
        print('results directory is created!')
    print(conf)
    if conf.skip_gen_msa:
        print('skipping generating msa...')
    else:
        run_msa(conf)
    if conf.skip_run_rf:
        print('skipping RoseTTAFold...')
    else:
        run_rfpd(conf)

def extract(conf: Conf) -> None:
    print('Extracting...')
    npzFile = conf.path.out + '/results/pred/*feature.npz'
    for npz in iglob(npzFile):
        exPair(npz)
        if conf.verbose:
            print('embedding of %s is extracted!' %npz.rsplit('/',1)[1])

#    cmd = f"./extract.sh {conf.path.out} {conf.path.scripts}"
#    print(cmd)
#    os.system(cmd)

if __name__ == "__main__":
    print(config)
    if config.skip_preprocess:
        print('skipping preprocessing...')
    else:
        Pre(config)
    runEmbed(config)
    if config.skip_extract:
        print('skipping extracting...')
    else:
        extract(config)
    #runEmbed(config)
    #extract(config)
