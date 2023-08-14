import os
from glob import iglob

from rfold.rconfig import Conf#, config
from Preprocess import check_fasta
from NPZtoPair import exPair
from rfold.gen_msa import gen_msa
from rfold.gen_embed import gen_embed
import tyro
# from rfscripts.rfpd import run_msa
# from rfscripts.rfpd import run_rfpd
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
    check_fasta(conf)

def runEmbed(conf: Conf) -> None:
    print('Embedding...')
    if not os.path.exists(conf.env.RF_RUNTIME_BASE):
        os.mkdir(conf.env.RF_RUNTIME_BASE)
        print('results directory is created!')
    # print(conf)
    if conf.gen_msa and not conf.run_rf:
        print('only generating msa...')
        gen_msa(conf)
        exit()
    if not conf.gen_msa and conf.run_rf:
        print('only running RoseTTAFold..')
        gen_embed(conf)
        return
    gen_msa(conf)
    gen_embed(conf)
    return


def extract(conf: Conf) -> None:
    print('Extracting...')
    npzFile = conf.env.RF_RUNTIME_BASE + '/pred/*feature.npz'
    for npz in iglob(npzFile):
        exPair(npz)
        if conf.verbose:
            print('embedding of %s is extracted!' %npz.rsplit('/',1)[1])

#    cmd = f"./extract.sh {conf.path.out} {conf.path.scripts}"
#    print(cmd)
#    os.system(cmd)

if __name__ == "__main__":
    config = tyro.cli(Conf)
    # print(config)
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
