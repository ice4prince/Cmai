import sys
import os

# getting the name of the directory
# where the this file is present.
current = os.path.dirname(os.path.realpath(__file__))

# Getting the parent directory name
# where the current directory is present.
parent = os.path.dirname(current)

# adding the parent directory to
# the sys.path.
sys.path.append(parent)

# importing
import numpy as np
# from .rfold.rconfig import Conf, config
from glob import iglob
#from optparse import OptionParser

def exPair(npzFile):
    features = np.load(npzFile)
    pair = features["pair"]
    np.save(npzFile.replace('feature.npz','pair.npy'),pair)

if __name__ == "__main__":
    from .rfold.rconfig import Conf, config
    npzF = conf.env.RF_RUNTIME_BASE + '/pred/*feature.npz'
    for npz in iglob(npzF):
        exPair(npz)
        print('embedding of %s is extracted!' %npz.rsplit('/',1)[1])
