#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 12 12:51:06 2019

@author: martinm
"""


import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os

import forfsmapping

#%% combine coverage files

'''
Run mapping to library variants in batches using the code provided in the file mapfsbins_splicinganddels
'''

lib300=pd.read_csv('./design/fslibraryall.csv', index_col='Unnamed: 0')

cov=pd.Series('',index=lib300.index)

for filename in os.listdir('./mapping/plus1'):
    if ('coveragePYTHON-' in filename):
        splitcov = pd.read_pickle('./mapping/plus1/' + filename)
        cov = cov.add(splitcov)

cov.to_pickle('./mapping/plus1/fscov_FSp1.pkl')



# mapping, should be parallelized
fsp1dnavarsdfud3, fsp1bincovdfud3=forfsmapping.map_intronsdels_fs_difflibupdown(cov, 3, 35790, 53598, kind='bins', fstype='p1')

fsp1dnavarsdfud3.to_pickle('./mapping/plus1/fscov_FSp1_DNAvariants_min3_difflibupdown_minlength10.pkl')
fsp1bincovdfud3.to_pickle('./mapping/plus1/fscov_FSp1_bincoverage_min3_difflibupdown_minlength10.pkl')


#%%

#%%

def apply_filters_p1(bincov, dnavars, name):
    forfsmapping.check_binprimercorrelation(bincov)
    covcomb=forfsmapping.combine_primers(bincov)
       
    covf2=forfsmapping.FilterMinBinCoverage(covcomb,2)    
    
    covf3=forfsmapping.FilterExtremeBins(covf2)                                                   
    
    covf4=forfsmapping.FilterIsolated(covf3)    
    

#%% check the relative number of reads from each bin and compare with the percentage of cells
    
    binreads300=covf4.sum()/covf4.sum().sum()*100
    
    binpercentages=pd.read_excel('./mapping/plus1/percentofcellssorted_fs2.xlsx', index_col='bin') 
    binpercentages.index=binreads300.index
    
    plt.figure(figsize=(4,4))
    plt.plot(binreads300,binpercentages,'.')
    plt.xlabel('number of reads')
    plt.ylabel('percent of cells')
    for x in binreads300.index:
        plt.annotate(x, xy=(binreads300.loc[x], binpercentages.loc[x]), fontsize=12)
    
    
    #%% correct for number of reads per bin and percent sorted
    
    
    correctionfactor = binpercentages.percentcells/(binreads300/binreads300.sum()*100)
    
    covf5 = forfsmapping.NormalizeByBinPercentages(covf4, correctionfactor)
    
    #sanity check
    binreads300norm=covf5.sum()/covf5.sum().sum()*100
    
    plt.figure(figsize=(4,4))
    plt.plot(binreads300norm,binpercentages,'.')
    plt.xlabel('number of reads')
    plt.ylabel('percent of cells')
    
    #%% expression estimates
    
    ratiocalcbins=pd.read_pickle('./mapping/plus1/ratiocalcbins1.pkl') 
    
    xval=ratiocalcbins['median'].apply(lambda x: np.log2(x))
    
    
    covnormp1variants=forfsmapping.NormalizeToReadNumbers(covf5)
    covnormp1variants.to_pickle('./mapping/plus1/'+name+'norm_p1variants.pkl')
    exprp1variants = forfsmapping.ExpCalc(covnormp1variants, xval)
    peaksp1variants=forfsmapping.findpeaks(covnormp1variants, xval)
    
    return exprp1variants.join(peaksp1variants).join(dnavars)

    #%%

#%%


exprp1variants3dlud=apply_filters_p1(fsp1bincovdfud3, fsp1dnavarsdfud3, 'fsp1_difflibupdown_min3reads_minlength10')

forp1dlud3=forfsmapping.classify_variants(exprp1variants3dlud, lib300, 1)
exprp1variants3dlud=exprp1variants3dlud.join(forp1dlud3)
exprp1variants3dlud.to_pickle('/net/mraid08/export/genie/Runs/Martin/fs_rev/Sample_FS1_v0_s/expr_fsm1_difflibupdown_min3reads_minlength10.pkl')

cov_wtmut_dlud3, vars_wtmut_dlud3=forfsmapping.create_wt_combmuts_bincov(exprp1variants3dlud, fsp1bincovdfud3)

exprp1comb=apply_filters_p1(cov_wtmut_dlud3, vars_wtmut_dlud3, 'fsp1_mutscombined_difflibupdown_min3reads_minlength10')

exprp1comb.to_pickle('./mapping/plus1/exprfsp1_mutscombined_difflibupdown_min3reads_minlength10variants.pkl')
