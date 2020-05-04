#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 18 15:26:14 2019

@author: martinm
"""


import pandas as pd
import numpy as np
import os

import forfsmapping
from maxentpy import maxent

martin='/net/mraid08/export/genie/Runs/Martin/'
#%% combine coverage files


'''
Run mapping to library variants in batches using the code provided in the file mapfsbins_splicinganddels
'''

lib300=pd.read_csv('./design/fslibraryall.csv', index_col='Unnamed: 0')

cov=pd.Series('',index=lib300.index)

for filename in os.listdir('./mapping/RNA'):
    if ('coveragePYTHON-' in filename):
        splitcov = pd.read_pickle('./mapping/RNA/' + filename)
        cov = cov.add(splitcov)

# mapping, should be parallelized
fsm1rnavars=forfsmapping.map_intronsdels_fs_difflibupdown(cov, 3, 35790, 53598, kind='RNA')

#%%

rnacounts=pd.DataFrame()
rnacounts['numberreads_cum']=fsm1rnavars.groupby('libindex').numberreads.sum()

fsm1rnavars['fraction_numberreads']=fsm1rnavars.index.map(lambda x: float(fsm1rnavars.numberreads[x])/rnacounts.loc[fsm1rnavars.libindex[x],'numberreads_cum'])

def classify_variant(variant):
    if variant=='wt':
        return 'wt'
    elif '_' in variant:
        return 'multiple'
    elif 'del' in variant:
        return str(int(variant.split('|')[-1])-int(variant.split('|')[0]))
    else:
        return 'pointmut'

fsm1rnavars['variantclass']=fsm1rnavars.variant.apply(lambda x: classify_variant(x))

fsm1rnavars['intronlength']=fsm1rnavars.variantclass.apply(lambda x: int(x) if x not in ['wt','multiple','pointmut'] else np.nan)


fsm1rnavars['intronseq']=fsm1rnavars.index.map(lambda x: lib300.varseq[fsm1rnavars.libindex[x]][int(fsm1rnavars.variant[x].split('|')[0])-5: int(fsm1rnavars.variant[x].split('|')[-1])+5] if fsm1rnavars.intronlength[x]>20 else np.nan)

fsm1rnavars['donorstrength']=fsm1rnavars.intronseq.map(lambda x: 
    np.max([maxent.score5(x[i-3:i+6]) for i in range(3,8)]) if (len(str(x))>28) else np.nan)
fsm1rnavars['acceptorstrength']=fsm1rnavars.intronseq.map(lambda x: 
    np.max([maxent.score3(x[-i-20:-i+3]) for i in range(4,8)]) if (len(str(x))>28) else np.nan)


rnacounts['numbersplicedreads_cum']=fsm1rnavars[fsm1rnavars.intronlength>20].groupby('libindex').numberreads.sum()
rnacounts['fractionnumbersplicedreads']=rnacounts.index.map(lambda x: rnacounts.numbersplicedreads_cum[x]/rnacounts.numberreads_cum[x])
rnacounts['fractionnumbersplicedreads'].replace(to_replace=np.nan, value=0, inplace=True)

rnacounts.to_pickle('./mapping/RNA/rnacounts_min3_dlud_minlength10.pkl')

fsm1rnavars.to_pickle('./mapping/RNA/fsm1rnavars_min3_dlud_minlength10.pkl')
