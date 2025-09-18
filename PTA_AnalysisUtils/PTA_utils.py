#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 10 12:55:29 2023

@author: kgrunthal
"""

import os, glob
import subprocess
import time
import numpy as np

from enterprise.pulsar import Pulsar





#------------------------------------------------------------------------------
def get_psrname(file,name_sep='.'):
    return file.split('/')[-1].split(name_sep)[0]
#------------------------------------------------------------------------------



#------------------------------------------------------------------------------
def load_pulsars(pardir, timdir, psr_list):
    pars = sorted(glob.glob(pardir + '/*.par'))
    parfiles = [f for f in pars if get_psrname(f) in psr_list]
    tims = sorted(glob.glob(timdir +'/*.tim' ) )
    timfiles = [f for f in tims if get_psrname(f) in psr_list]
    
    PSRs = []
    for p,t in zip(parfiles,timfiles):
        psr = Pulsar(p, t, ephem='DE440')
        PSRs.append(psr)
    
    return PSRs
#------------------------------------------------------------------------------



#------------------------------------------------------------------------------
def save_complex_matrix(filename, matrix):
    np.savetxt(filename, matrix)
    if os.path.isfile(filename):
        subprocess.run("sed -i s/[)(]//g {}".format(filename).split(' ') )
    else:
        time.sleep(1)
    return None
#------------------------------------------------------------------------------
