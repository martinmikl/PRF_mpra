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

for filename in os.listdir('./mapping/minus1'):
    if ('coveragePYTHON-' in filename):
        splitcov = pd.read_pickle('./mapping/minus1/' + filename)
        cov = cov.add(splitcov)

cov.to_pickle('./mapping/minus1/fscov_FSm1.pkl')


# mapping, should be parallelized
fsm1dnavarsdfud3, fsm1bincovdfud3=forfsmapping.map_intronsdels_fs_difflibupdown(cov, 3, 35790, 53598, kind='bins')

fsm1dnavarsdfud3.to_pickle('./mapping/minus1/fscov_FSm1_DNAvariants_min3_difflibupdown_minlength10.pkl')
fsm1bincovdfud3.to_pickle('./mapping/minus1/fscov_FSm1_bincoverage_min3_difflibupdown_minlength10.pkl')


#%%
#%%

def apply_filters_m1(bincov, dnavars, name):
    forfsmapping.check_binprimercorrelation(bincov)
    covcomb=forfsmapping.combine_primers(bincov)
       
    covf2=forfsmapping.FilterMinBinCoverage(covcomb,2)    
    
    covf3=forfsmapping.FilterExtremeBins(covf2)                                                   
    
    covf4=forfsmapping.FilterIsolated(covf3)    
    

#%% check the relative number of reads from each bin and compare with the percentage of cells
    
    binreads300=covf4.sum()/covf4.sum().sum()*100
    
    binpercentages=pd.read_excel('./mapping/minus1/percentofcellssorted_fs1.xlsx', index_col='bin') 
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
    
    ratiocalcbins=pd.read_pickle('./mapping/minus1/ratiocalcbins1.pkl') 
    
    xval=ratiocalcbins['median'].apply(lambda x: np.log2(x))
    
    
    covnormm1variants=forfsmapping.NormalizeToReadNumbers(covf5)
    covnormm1variants.to_pickle('./mapping/minus1/'+name+'norm_m1variants.pkl')
    exprm1variants = forfsmapping.ExpCalc(covnormm1variants, xval)
    peaksm1variants=forfsmapping.findpeaks(covnormm1variants, xval)
    
    return exprm1variants.join(peaksm1variants).join(dnavars)

    #%%
    
exprm1variants3dlud=apply_filters_m1(fsm1bincovdfud3, fsm1dnavarsdfud3, 'fsm1_difflibupdown_min3reads_minlength10')

form1dlud3=forfsmapping.classify_variants(exprm1variants3dlud, lib300, -1)
exprm1variants3dlud=exprm1variants3dlud.join(form1dlud3)
exprm1variants3dlud.to_pickle('./mapping/minus1/expr_fsm1_difflibupdown_min3reads_minlength10.pkl')

cov_wtmut_dlud3, vars_wtmut_dlud3=forfsmapping.create_wt_combmuts_bincov(exprm1variants3dlud, fsm1bincovdfud3)
#vars_wtmut_dlud3.to_pickle('./mapping/minus1/fscov_FSm1_mutscombined_DNAvariants_min3_difflibupdown_minlength10.pkl')
#cov_wtmut_dlud3.to_pickle('./mapping/minus1/fscov_FSm1_mutscombined_bincoverage_min3_difflibupdown_minlength10.pkl')

exprm1comb=apply_filters_m1(cov_wtmut_dlud3, vars_wtmut_dlud3, 'fsm1_mutscombined_difflibupdown_min3reads_minlength10')

exprm1comb.to_pickle('./mapping/minus1/exprfsm1_mutscombined_difflibupdown_min3reads_minlength10variants.pkl')

