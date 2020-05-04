#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Aug  9 12:45:55 2018

@author: martinm
"""

         
#%%
#%%
import pandas as pd
from sklearn.cross_validation import train_test_split

#%%
#%% 
fsdf=pd.read_pickle('./dataframes/fsdf.pkl')

fsfilt=fsdf[(fsdf.numberreadsm_wt>=20)&(fsdf.gfpm_wt<25)]

# drop indexes for duplicated sequences

indexes=[]
samevarseq=[]
for var, group in fsfilt.groupby(by='varseq162'):
    if (len(group)>1):
        samevarseq.append(var)
        indexes+=list(fsfilt[fsfilt.varseq162==var].sample(n=len(fsfilt[fsfilt.varseq162==var])-1, random_state=77).index)
        
fsforml=fsfilt.drop(indexes) 


for i in fsforml.index:
    if (fsforml.varseq162[i].upper().find('N')>-1)|(fsforml.varseq162[i].upper().find('R')>-1)|(fsforml.varseq162[i].upper().find('Y')>-1):
        fsforml.drop(i, inplace=True)
        print(i)

fsdesforml=fsforml[fsforml.library=='fs_designed']
fsdesforml['shifting']=fsdesforml.gfpm_wt.apply(lambda x: 1 if x>1.3 else 0)

hivdbvarforml=fsforml[(fsforml.subset=='hivdbvar')]

hivdbvarforml['validseq']=hivdbvarforml.varseq162.apply(lambda s: all(i in 'ACTG' for i in s.upper()))
hivdbvarforml['shifting']=hivdbvarforml.gfpm_wt.apply(lambda x: 1 if x>1.3 else 0)

fsdesforml.to_pickle('./dataframes/ml/fsdesforml.pkl')
hivdbvarforml.to_pickle('./dataframes/ml/hivdbvarforml.pkl')

#####

fsdesmlx,fsdestestx, fsdesmly,fsdestesty=train_test_split(fsdesforml, 
                                     fsdesforml['gfpm_wt'], \
                                                test_size=0.1, random_state=0)

fsdesmlx.to_pickle('./dataframes/ml/fsdesmlx.pkl')
fsdestestx.to_pickle('./dataframes/ml/fsdestestx.pkl')
fsdesmly.to_pickle('./dataframes/ml/fsdesmly.pkl')
fsdestesty.to_pickle('./dataframes/ml/fsdestesty.pkl')
  

#####

hivmlx,hivtestx, hivmly,hivtesty=train_test_split(hivdbvarforml, hivdbvarforml['gfpm_wt'], \
                                                test_size=0.2, random_state=0)

hivmlx.to_pickle('./dataframes/ml/hivmlx.pkl')
hivtestx.to_pickle('./dataframes/ml/hivtestx.pkl')
hivmly.to_pickle('./dataframes/ml/hivmly.pkl')
hivtesty.to_pickle('./dataframes/ml/hivtesty.pkl')
      
#####


