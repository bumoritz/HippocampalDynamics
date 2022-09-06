# -*- coding: utf-8 -*-
"""
Created on Sun Apr 25 14:34:50 2021
@author: Moritz
"""


### IMPORT PYTHON PACKAGES ###

import os, sys
import numpy as np
import scipy.io as sio



##############################
### ----- PARAMETERS ----- ### 
##############################


animal      = sys.argv[1] #'Carlo'
date        = sys.argv[2] #'20210315'
base_path   = sys.argv[3] #r'E:\Data'

# options
do_clone_ops            = 1
do_clone_stat           = 1


##############################
### ---------------------- ### 
##############################



### DATA FOLDER ###
folder = os.path.join(base_path,date[0:4],date[0:4]+'-'+date[4:6],date[0:4]+'-'+date[4:6]+'-'+date[6:8],animal,'Imaging','suite2p','plane0')
print('--- Cloning ops and stat... [convert_ops_and_stat.py]')


### CLONE OPS AND STAT ###

if do_clone_ops:
    ops = np.load(os.path.join(folder,'ops.npy'),allow_pickle=True).item()
    sio.savemat(os.path.join(folder,'_ops.mat'),{'ops':ops})
if do_clone_stat:
    stat = np.load(os.path.join(folder,'stat.npy'),allow_pickle=True)
    sio.savemat(os.path.join(folder,'_stat.mat'),{'stat':stat})






