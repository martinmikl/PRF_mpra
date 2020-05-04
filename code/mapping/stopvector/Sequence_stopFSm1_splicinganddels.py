#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 12 12:51:06 2019

@author: martinm
"""


import pandas as pd
import matplotlib.pyplot as plt
import os

import forfsmapping

#%% combine coverage files


'''
Run mapping to library variants in batches using the code provided in the file mapfsbins_splicinganddels
'''

lib300=pd.read_csv('./design/fslibraryall.csv', index_col='Unnamed: 0')

cov=pd.Series('',index=lib300.index)

for filename in os.listdir('./mapping/stopvector'):
    if ('coveragePYTHON-' in filename):
        splitcov = pd.read_pickle('./mapping/stopvector/' + filename)
        cov = cov.add(splitcov)

cov.to_pickle('./mapping/stopvector/fscov_stopFS.pkl')


# mapping, should be parallelized
fsstopm1dnavars3, fsstopm1bincov3=forfsmapping.map_intronsdels_fs_difflibupdown(cov, 3, 35790, 53598, kind='stop')

fsstopm1dnavars3.to_pickle('./mapping/stopvector/fscov_FSstop_DNAvariants_min3_difflibupdown.pkl')
fsstopm1bincov3.to_pickle('./mapping/stopvector/fscov_FSstop_bincoverage_min3_difflibupdown.pkl')


#%%

#%%

def apply_filters(bincov, dnavars, name):
       
    covf2=forfsmapping.FilterMinBinCoverage(bincov,2)    
    
    covf3=forfsmapping.FilterExtremeBins(covf2)                                                   
    
    covf4=forfsmapping.FilterIsolated12bins(covf3)    
    

#%% check the relative number of reads from each bin and compare with the percentage of cells
    
    binreads300=covf4.sum()/covf4.sum().sum()*100
    
    binpercentages=pd.read_excel('./mapping/stopvector/percentofcellssorted_fsstop.xlsx', index_col='bin') 
    binpercentages.index=binreads300.index
    
    plt.figure(figsize=(4,4))
    plt.plot(binreads300,binpercentages,'.')
    plt.xlabel('number of reads')
    plt.ylabel('percent of cells')
    for x in binreads300.index:
        plt.annotate(x, xy=(binreads300.loc[x], binpercentages.loc[x]), fontsize=12)
    
    
    #%% correct for number of reads per bin and percent sorted
    
    
    correctionfactor = binpercentages.fraction/(binreads300/binreads300.sum()*100)
    
    covf5 = forfsmapping.NormalizeByBinPercentages(covf4, correctionfactor)
    
    #sanity check
    binreads300norm=covf5.sum()/covf5.sum().sum()*100
    
    plt.figure(figsize=(4,4))
    plt.plot(binreads300norm,binpercentages,'.')
    plt.xlabel('number of reads')
    plt.ylabel('percent of cells')
    
    #%% expression estimates
    
    xval=pd.read_pickle('./mapping/stopvector/xvals_corrected.pkl') 
    
    covnormstopvariants=forfsmapping.NormalizeToReadNumbers12bins(covf5)
    covnormstopvariants.to_pickle('./mapping/stopvector/'+name+'norm_m1variants.pkl')
    exprm1variants = forfsmapping.ExpCalc(covnormstopvariants, xval)
    peaksm1variants=forfsmapping.findpeaks(covnormstopvariants, xval)
    
    return exprm1variants.join(peaksm1variants).join(dnavars)

    #%%
    

exprstopvariants3dlud=apply_filters(fsstopm1bincov3, fsstopm1dnavars3, 'fsstop_difflibupdown_min3reads_minlength10')

form1dlud3=forfsmapping.classify_variants(exprstopvariants3dlud, lib300, -1)
exprstopvariants3dlud=exprstopvariants3dlud.join(form1dlud3)
exprstopvariants3dlud.to_pickle('./mapping/stopvector/exprfsstop_difflibupdown_min3reads_minlength10variants.pkl')

cov_wtmut_dlud3, vars_wtmut_dlud3=forfsmapping.create_wt_combmuts_bincov(exprstopvariants3dlud, fsstopm1bincov3)
vars_wtmut_dlud3.to_pickle('./mapping/stopvector/fscov_FSstop_mutscombined_DNAvariants_min3_difflibupdown_minlength10.pkl')
cov_wtmut_dlud3.to_pickle('./mapping/stopvector/fscov_FSstop_mutscombined_bincoverage_min3_difflibupdown_minlength10.pkl')

exprstopcomb=apply_filters(cov_wtmut_dlud3, vars_wtmut_dlud3, 'fsstop_mutscombined_difflibupdown_min3reads_minlength10')

exprstopcomb.to_pickle('./mapping/stopvector/exprfsstop_mutscombined_difflibupdown_min3reads_minlength10variants.pkl')
    
    
