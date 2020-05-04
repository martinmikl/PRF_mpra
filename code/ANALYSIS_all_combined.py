#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 20 12:32:09 2019

@author: martinm
"""



import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import RNA
from Bio.Seq import Seq

sns.set_context('talk')


ratiocalcbins=pd.read_pickle('./mapping/minus1/ratiocalcbins1.pkl') 

xval=ratiocalcbins['median'].apply(lambda x: np.log2(x))

xvalstop=pd.read_pickle('./mapping/stopvector/xvals_corrected.pkl') 


# load library

fs=pd.read_csv('./design/fslibraryall.csv', index_col='Unnamed: 0')

fs['varseq162']=fs.varseq.apply(lambda x: x[30:-18].upper())


# new mapping - FSm1

exprm1comb=pd.read_pickle('./mapping/minus1/exprfsm1_mutscombined_difflibupdown_min3reads_minlength10variants.pkl')
covm1comb=pd.read_pickle('./mapping/minus1/fsm1_mutscombined_difflibupdown_min3reads_minlength10norm_m1variants.pkl')

exprm1vars=pd.read_pickle('./mapping/minus1/exprfsm1_difflibupdown_min3reads_minlength10variants.pkl')

# new mapping - FSp1

exprp1comb=pd.read_pickle('./mapping/plus1/exprfsp1_mutscombined_difflibupdown_min3reads_minlength10variants.pkl')
covp1comb=pd.read_pickle('./mapping/plus1/fsp1_mutscombined_difflibupdown_min3reads_minlength10norm_p1variants.pkl')

exprp1vars=pd.read_pickle('./mapping/plus1/exprfsp1_difflibupdown_min3reads_minlength10variants.pkl')

# new mapping - stopFSm1

exprstopcomb=pd.read_pickle('./mapping/stopvector/exprfsstop_mutscombined_difflibupdown_min3reads_minlength10variants.pkl')
covstopcomb=pd.read_pickle('./mapping/stopvector/fsstop_mutscombined_difflibupdown_min3reads_minlength10norm_m1variants.pkl')

exprstopvars=pd.read_pickle('./mapping/stopvector/exprfsstop_difflibupdown_min3reads_minlength10variants.pkl')


# RNA

rnacounts=pd.read_pickle('./mapping/RNA/rnacounts_min3_dlud_minlength10.pkl')

fsm1rnavars=pd.read_pickle('./mapping/RNA/fsm1rnavars_min3_dlud_minlength10.pkl')



#%%

def isbonafide(var):
    istrue=False
    if '_' not in var:
        if 'del' in var:
            coord=var.split('|del|')
            if (int(coord[1])-int(coord[0])==2):
                istrue=True
                
    return istrue
        
exprm1vars['bonafide']=exprm1vars.variant.apply(lambda x: isbonafide(x))

        
        
#%% make one df


wavm=exprm1comb.pivot(index='libindex', columns='variant', values='wav')
wavm.columns=['wavm_'+str(x) for x in wavm.columns]
peaksm=exprm1comb.pivot(index='libindex', columns='variant', values='smoothedpeaks')
peaksm.columns=['peaksm_'+str(x) for x in peaksm.columns]
readsm=exprm1comb.pivot(index='libindex', columns='variant', values='numberreads')
readsm.columns=['numberreadsm_'+str(x) for x in readsm.columns]

wavp=exprp1comb.pivot(index='libindex', columns='variant', values='wav')
wavp.columns=['wavp_'+str(x) for x in wavp.columns]
peaksp=exprp1comb.pivot(index='libindex', columns='variant', values='smoothedpeaks')
peaksp.columns=['peaksp_'+str(x) for x in peaksp.columns]
readsp=exprm1comb.pivot(index='libindex', columns='variant', values='numberreads')
readsp.columns=['numberreadsp_'+str(x) for x in readsp.columns]


wavstop=exprstopcomb.pivot(index='libindex', columns='variant', values='wav')
wavstop.columns=['wavstop_'+str(x) for x in wavstop.columns]
peaksstop=exprstopcomb.pivot(index='libindex', columns='variant', values='smoothedpeaks')
peaksstop.columns=['peaksstop_'+str(x) for x in peaksstop.columns]
readsstop=exprm1comb.pivot(index='libindex', columns='variant', values='numberreads')
readsstop.columns=['numberreadsstop_'+str(x) for x in readsstop.columns]

rnaexp=fsm1rnavars[(fsm1rnavars.variant=='wt')][['numberreads','libindex']].set_index('libindex')

fsdf=pd.concat([fs, wavm, wavp, wavstop, peaksm, peaksp, peaksstop, readsm, readsp, readsstop, rnaexp, rnacounts], axis=1)

fsdf['rnaperdna']=np.log2(fsdf.numberreads.replace(to_replace=np.nan, value=0)/(fsdf.numberreadsp_wt.replace(to_replace=np.nan, value=0) +fsdf.numberreadsp_wt.replace(to_replace=np.nan, value=0)))
fsdf['rnaperdna']=fsdf['rnaperdna'].replace(to_replace=np.inf, value=np.nan).replace(to_replace=-np.inf, value=np.nan)

fsdf['wavm_byrna']=np.log2(np.exp2(fsdf['wavm_wt'])/np.exp2(fsdf['rnaperdna']))
fsdf['wavm_byrna'].replace(to_replace=np.inf, value=np.nan).replace(to_replace=-np.inf, value=np.nan).hist(bins=100)
fsdf['wavm_wt'].hist(bins=100)

fsdf['fsevent']=fsdf.index.map(lambda x: ' - '.join(pd.Series([fsdf.loc[x,'organism'], fsdf.loc[x,'gene']]).dropna()))
fsdf.rename(columns={'sequence ID':'seqid'}, inplace=True)

shortenfsevent={'PLRV luteovirus - RDRP and coat':'PLRV luteovirus',
'herpes simplex virus - thymidine kinase wt':'herpes simplex - thymidine kinase',
'human - peg10 human homolog of edr':'human - PEG10',
'influenza a virus +1':'influenza a virus',
'human t cell lymphotropic virus':'human T-lymphotropic virus',
'human ccr5':'human CCR5'}

fsdf.fsevent.replace(shortenfsevent, inplace=True)


absmin=fsdf[['wavm_wt','wavp_wt','wavstop_wt']].min().min()
absmax=fsdf[['wavm_wt','wavp_wt','wavstop_wt']].max().max()

for col in ['wavm_wt','wavm_del','wavm_ins','wavm_del_withmultiple','wavm_ins_withmultiple',
            'wavp_wt','wavp_del','wavp_ins','wavp_del_withmultiple','wavp_ins_withmultiple',
            'wavstop_wt','wavstop_del','wavstop_ins','wavstop_del_withmultiple','wavstop_ins_withmultiple']:
    fsdf['gfp'+col[3:]]=fsdf[col].apply(lambda x: (2**x-2**absmin)/2**absmax*100)     
    
    
    
exprm1vars.columns=[str(x)+'_m' for x in exprm1vars.columns]
exprp1vars.columns=[str(x)+'_p' for x in exprp1vars.columns]
exprstopvars.columns=[str(x)+'_stop' for x in exprstopvars.columns]

fsvarsdf=pd.concat([exprm1vars, exprp1vars, exprstopvars, fsm1rnavars], axis=1)
fsvarsdf['rnaperdna']=np.log2(fsvarsdf.numberreads.replace(to_replace=np.nan, value=0)/(fsvarsdf.numberreads_p.replace(to_replace=np.nan, value=0) +fsvarsdf.numberreads_m.replace(to_replace=np.nan, value=0)))


#%% calculate wt values


wtvalues=pd.DataFrame()
for var in fsdf[fsdf.subset=='fsmain'].fsevent.unique():
    vrsq=fsdf[(fsdf.subset=='fsmain')&(fsdf.fsevent==var)].varseq162.values[0]
    indexes=list(fsdf[(fsdf.wavm_wt<12)&(fsdf.varseq162==vrsq)&(fsdf.peaksm_wt==1)&(fsdf.numberreadsm_wt>=20)].index)  
    wtvalues.loc[var,'wavm']=fsdf.loc[indexes,'wavm_wt'].median()
    indexes=list(fsdf[(fsdf.wavp_wt<12)&(fsdf.varseq162==vrsq)&(fsdf.peaksp_wt==1)&(fsdf.numberreadsp_wt>=20)].index)  
    wtvalues.loc[var,'wavp']=fsdf.loc[indexes,'wavp_wt'].median()
    indexes=list(fsdf[(fsdf.wavstop_wt<12)&(fsdf.varseq162==vrsq)&(fsdf.peaksstop_wt==1)&(fsdf.numberreadsstop_wt>=20)].index)  
    wtvalues.loc[var,'wavstop']=fsdf.loc[indexes,'wavstop_wt'].median()

    indexes=list(fsdf[(fsdf.wavm_wt<12)&(fsdf.varseq162==vrsq)&(fsdf.peaksm_wt==1)&(fsdf.numberreadsm_wt>=20)].index)  
    wtvalues.loc[var,'gfpm']=fsdf.loc[indexes,'gfpm_wt'].median()
    indexes=list(fsdf[(fsdf.wavp_wt<12)&(fsdf.varseq162==vrsq)&(fsdf.peaksp_wt==1)&(fsdf.numberreadsp_wt>=20)].index)  
    wtvalues.loc[var,'gfpp']=fsdf.loc[indexes,'gfpp_wt'].median()
    indexes=list(fsdf[(fsdf.wavstop_wt<12)&(fsdf.varseq162==vrsq)&(fsdf.peaksstop_wt==1)&(fsdf.numberreadsstop_wt>=20)].index)  
    wtvalues.loc[var,'gfpstop']=fsdf.loc[indexes,'gfpstop_wt'].median()
    
wtvalues.to_pickle('./dataframes/wtvalues.pkl')
   
fsdf['normwavm_wt']=fsdf.index.map(lambda x: fsdf.wavm_wt[x] - wtvalues.wavm[fsdf.fsevent[x]] if fsdf.fsevent[x] in wtvalues.index else np.nan)
fsdf['normwavp_wt']=fsdf.index.map(lambda x: fsdf.wavp_wt[x] - wtvalues.wavp[fsdf.fsevent[x]] if fsdf.fsevent[x] in wtvalues.index else np.nan)
    

fsdf['normgfpm_wt']=fsdf.index.map(lambda x: fsdf.gfpm_wt[x]-wtvalues.loc[fsdf.fsevent[x],'gfpm'] if 
    fsdf.fsevent[x] in wtvalues.index else np.nan)
fsdf['normgfpp_wt']=fsdf.index.map(lambda x: fsdf.gfpp_wt[x]-wtvalues.loc[fsdf.fsevent[x],'gfpp'] if 
    fsdf.fsevent[x] in wtvalues.index else np.nan)


fsdf['percgfpm_wt']=fsdf.index.map(lambda x: (fsdf.gfpm_wt[x])/(wtvalues.gfpm[fsdf.fsevent[x]])*100 if fsdf.fsevent[x] in wtvalues.index else np.nan)
fsdf['percgfpp_wt']=fsdf.index.map(lambda x: (fsdf.gfpp_wt[x])/(wtvalues.gfpp[fsdf.fsevent[x]])*100 if fsdf.fsevent[x] in wtvalues.index else np.nan)


# Secondary structure

fsdf['dg_downstream']=fsdf.varseq162.apply(lambda x: RNA.fold(x[42:])[1])

for length in range(0,70,10):
    fsdf['dgup'+str(length)]=fsdf.varseq.apply(lambda x: RNA.fold(x[length:30+42])[1])
for length in range(20,130,10):
    fsdf['dg'+str(length)]=fsdf.varseq.apply(lambda x: RNA.fold(x[30+42:30+42+length])[1])
    
# endogenous slipperysite

def endoslip(idx):
    try:
        if fsdf.varseq162[idx][35:42]==fsdf[(fsdf.subset=='fsmain')&(fsdf.fsevent==fsdf.fsevent[idx])].varseq162.values[0][35:42]:
            return True
        else:
            return False
    except:
        return np.nan
        
        
fsdf['endoslippery']=fsdf.index.map(lambda x: endoslip(x))
  

fsdf['laststopinframem1']=fsdf.varseq162.apply(lambda x: str(Seq(x[2:]).translate()).rfind('*'))
fsdf['laststopinframep1']=fsdf.varseq162.apply(lambda x: str(Seq(x[1:]).translate()).rfind('*'))

# Save Dataframes
   
fsdf.to_pickle('./dataframes/fsdf.pkl')
fsvarsdf.to_pickle('./dataframes/fsvarsdf.pkl')

   
fsdf.to_csv('./dataframes/fsdf.csv')
fsvarsdf.to_csv('./dataframes/fsvarsdf.csv')
    
    
# Save Supplementary Data 
  
fsdf[['subset','fsevent','changes',
      'wavm_wt','gfpm_wt','percgfpm_wt','peaksm_wt','numberreadsm_wt', 
      'wavp_wt','gfpp_wt','percgfpp_wt','peaksp_wt','numberreadsp_wt',
      'wavstop_wt','gfpstop_wt','peaksstop_wt','numberreadsstop_wt',
      'numberreads_cum', 'numbersplicedreads_cum','fractionnumbersplicedreads',
      'varseq','varseq162']].to_csv('./tables/SupplementaryData2.csv')
    
    
    
#%% Compare readouts - minus 1 vs plus 1

f=plt.figure(figsize=(4,4))
plt.scatter(fsdf[(fsdf.gfpm_wt<25)&(fsdf.gfpp_wt<25)&(fsdf.numberreadsm_wt>=10)&(fsdf.peaksm_wt==1)&(fsdf.numberreadsp_wt>=10)&(fsdf.peaksp_wt==1)].gfpm_wt, 
            fsdf[(fsdf.gfpm_wt<25)&(fsdf.gfpp_wt<25)&(fsdf.numberreadsm_wt>=10)&(fsdf.peaksm_wt==1)&(fsdf.numberreadsp_wt>=10)&(fsdf.peaksp_wt==1)].gfpp_wt, alpha=0.1)
plt.xlim(-0.5,25.5)
plt.ylim(-0.5,25.5)
plt.axhspan(-0.5, 1.3, color=sns.xkcd_rgb['medium blue'], alpha=0.2)
plt.axvspan(-0.5, 1.3, color=sns.xkcd_rgb['medium blue'], alpha=0.2)
f.savefig('./figures/gfpm_wt_vs_gfpp_wt.png', \
      dpi = 300, format='png', bbox_inches='tight', frameon=True)

f=plt.figure(figsize=(4,4))
plt.scatter(fsdf[(fsdf.gfpm_wt<25)&(fsdf.gfpp_wt<25)&(fsdf.numberreadsm_wt>=10)&(fsdf.peaksm_wt==1)&(fsdf.numberreadsp_wt>=10)&(fsdf.peaksp_wt==1)].gfpm_wt, 
            fsdf[(fsdf.gfpm_wt<25)&(fsdf.gfpp_wt<25)&(fsdf.numberreadsm_wt>=10)&(fsdf.peaksm_wt==1)&(fsdf.numberreadsp_wt>=10)&(fsdf.peaksp_wt==1)].gfpp_wt, alpha=0.1)
plt.xlim(-0.2,5.2)
plt.ylim(-0.2,5.2)
plt.axhspan(-0.2, 1.3, color=sns.xkcd_rgb['medium blue'], alpha=0.2)
plt.axvspan(-0.2, 1.3, color=sns.xkcd_rgb['medium blue'], alpha=0.2)
f.savefig('./figures/gfpm_wt_vs_gfpp_wt_max5.png', \
      dpi = 300, format='png', bbox_inches='tight', frameon=True)


# experimental vs stop vector


f=plt.figure(figsize=(4,4))
plt.scatter(fsdf[(fsdf.gfpm_wt<25)&(fsdf.gfpstop_wt<25)&(fsdf.numberreadsm_wt>=10)&(fsdf.peaksm_wt==1)&(fsdf.numberreadsstop_wt>=10)&(fsdf.peaksstop_wt==1)].gfpm_wt, 
            fsdf[(fsdf.gfpm_wt<25)&(fsdf.gfpstop_wt<25)&(fsdf.numberreadsm_wt>=10)&(fsdf.peaksm_wt==1)&(fsdf.numberreadsstop_wt>=10)&(fsdf.peaksstop_wt==1)].gfpstop_wt, alpha=0.1)
plt.xlim(-0.5,25.5)
plt.ylim(-0.5,25.5)
plt.axhspan(-0.5, 1.3, color=sns.xkcd_rgb['medium blue'], alpha=0.2)
plt.axvspan(-0.5, 1.3, color=sns.xkcd_rgb['medium blue'], alpha=0.2)
f.savefig('./figures/gfpm_wt_vs_gfpstop_wt.png', \
      dpi = 300, format='png', bbox_inches='tight', frameon=True)

f=plt.figure(figsize=(4,4))
plt.scatter(fsdf[(fsdf.gfpm_wt<25)&(fsdf.gfpstop_wt<25)&(fsdf.numberreadsm_wt>=10)&(fsdf.peaksm_wt==1)&(fsdf.numberreadsstop_wt>=10)&(fsdf.peaksstop_wt==1)].gfpm_wt, 
            fsdf[(fsdf.gfpm_wt<25)&(fsdf.gfpstop_wt<25)&(fsdf.numberreadsm_wt>=10)&(fsdf.peaksm_wt==1)&(fsdf.numberreadsstop_wt>=10)&(fsdf.peaksstop_wt==1)].gfpstop_wt, alpha=0.1)
plt.xlim(-0.2,5.2)
plt.ylim(-0.2,5.2)
plt.axhspan(-0.2, 1.3, color=sns.xkcd_rgb['medium blue'], alpha=0.2)
plt.axvspan(-0.2, 1.3, color=sns.xkcd_rgb['medium blue'], alpha=0.2)
f.savefig('./figures/gfpm_wt_vs_gfpstop_wt_max5.png', \
      dpi = 300, format='png', bbox_inches='tight', frameon=True)


# test for splicing

fsdf['donorstrength']=fsdf.index.map(lambda x: fsm1rnavars[fsm1rnavars.libindex==x].donorstrength.max())
fsdf['acceptorstrength']=fsdf.index.map(lambda x: fsm1rnavars[fsm1rnavars.libindex==x].acceptorstrength.max())
fsdf['donorandacceptorstrength']=fsdf['donorstrength']+fsdf['acceptorstrength']


def countsplicedreads(var):
    rnavrs=list(fsm1rnavars[(fsm1rnavars.libindex==var)&(fsm1rnavars.intronlength>20)].variant)
    m1vrs=list(exprm1vars[exprm1vars.libindex==var].variant)
    p1vrs=list(exprp1vars[exprp1vars.libindex==var].variant)
    readcount=0
    for i in rnavrs:
        if i not in m1vrs+p1vrs:
            readcount+=fsm1rnavars.numberreads[str(var)+'_'+i]
    return readcount

fsdf['numbersplicedreadsnotindna_cum']=fsdf.index.map(lambda x: countsplicedreads(x)) 
fsdf['fractionnumbersplicedreadsnotindna']=fsdf.index.map(lambda x: fsdf.numbersplicedreadsnotindna_cum[x]/rnacounts.numberreads_cum[float(x)]
    if float(x) in rnacounts.index else np.nan)
fsdf['fractionnumbersplicedreadsnotindna'].replace(to_replace=np.nan, value=0, inplace=True)


def comparernaanddna(var):
    rnavrs=list(fsm1rnavars[(fsm1rnavars.libindex==var)&(fsm1rnavars.intronlength>20)].variant)
    m1vrs=list(exprm1vars[exprm1vars.libindex==var].variant)
    p1vrs=list(exprp1vars[exprp1vars.libindex==var].variant)
    isindna=True
    for i in rnavrs:
        if i not in m1vrs+p1vrs:
            isindna=False
    return isindna
               
fsdf['rnavariantindna']=fsdf.index.map(lambda x: comparernaanddna(x))

f=plt.figure(figsize=(4,4))
plt.scatter(fsdf[(fsdf.gfpm_wt<25)&(fsdf.numberreadsm_wt>=20)&(fsdf.peaksm_wt==1)].gfpm_wt, 
            fsdf[(fsdf.gfpm_wt<25)&(fsdf.numberreadsm_wt>=20)&(fsdf.peaksm_wt==1)].fractionnumbersplicedreadsnotindna, alpha=0.1)
plt.scatter(fsdf[((fsdf.donorstrength>0)|(fsdf.acceptorstrength>0))&(fsdf.gfpm_wt<25)&(fsdf.numberreadsm_wt>=20)&(fsdf.peaksm_wt==1)].gfpm_wt, 
            fsdf[((fsdf.donorstrength>0)|(fsdf.acceptorstrength>0))&(fsdf.gfpm_wt<25)&(fsdf.numberreadsm_wt>=20)&(fsdf.peaksm_wt==1)].fractionnumbersplicedreadsnotindna, 
                 color=sns.xkcd_rgb['red'], alpha=0.3)
plt.axvspan(-0.5, 1.3, color=sns.xkcd_rgb['medium blue'], alpha=0.2)
plt.xlim(-0.5,25.5)
plt.ylim(-0.03,1.03)
f.savefig('./figures/gfpm_wt_vs_fractionnumbersplicedreads_withssstrength.png', \
      dpi = 300, format='png', bbox_inches='tight', frameon=True)


f=plt.figure(figsize=(4,4))
plt.scatter(fsdf[(fsdf.gfpp_wt<25)&(fsdf.numberreadsp_wt>=20)&(fsdf.peaksp_wt==1)].gfpp_wt, 
            fsdf[(fsdf.gfpp_wt<25)&(fsdf.numberreadsp_wt>=20)&(fsdf.peaksp_wt==1)].fractionnumbersplicedreadsnotindna, alpha=0.1)
plt.scatter(fsdf[((fsdf.donorstrength>0)|(fsdf.acceptorstrength>0))&(fsdf.gfpp_wt<25)&(fsdf.numberreadsp_wt>=20)&(fsdf.peaksp_wt==1)].gfpp_wt, 
            fsdf[((fsdf.donorstrength>0)|(fsdf.acceptorstrength>0))&(fsdf.gfpp_wt<25)&(fsdf.numberreadsp_wt>=20)&(fsdf.peaksp_wt==1)].fractionnumbersplicedreadsnotindna,
                 color=sns.xkcd_rgb['red'], alpha=0.3)
plt.axvspan(-0.5, 1.3, color=sns.xkcd_rgb['medium blue'], alpha=0.2)
plt.xlim(-0.5,25.5)
plt.ylim(-0.03,1.03)
f.savefig('./figures/gfpp_wt_vs_fractionnumbersplicedreads_withssstrength.png', \
      dpi = 300, format='png', bbox_inches='tight', frameon=True)

fsdf[(fsdf.numberreadsm_wt>=20)&(fsdf.peaksm_wt==1)].fractionnumbersplicedreadsnotindna.dropna().apply(lambda x: 1 if x>0.01 else 0).sum()/float(fsdf[(fsdf.gfpm_wt<25)&(fsdf.numberreadsm_wt>=20)&(fsdf.peaksm_wt==1)].gfpm_wt.dropna().shape[0])
# 0.029542920847268672
fsdf[(fsdf.numberreadsm_wt>=20)&(fsdf.peaksm_wt==1)].fractionnumbersplicedreadsnotindna.dropna().apply(lambda x: 1 if x>0.01 else 0).sum()
# 163/6599
fsdf[(fsdf.numberreadsp_wt>=20)&(fsdf.peaksp_wt==1)].fractionnumbersplicedreadsnotindna.dropna().apply(lambda x: 1 if x>0.01 else 0).sum()/float(fsdf[(fsdf.gfpp_wt<25)&(fsdf.numberreadsp_wt>=20)&(fsdf.peaksp_wt==1)].gfpp_wt.dropna().shape[0])
# 0.048065650644783117 = 82/1706
fsdf[((fsdf.donorstrength>0)|(fsdf.acceptorstrength>0))&(fsdf.numberreadsm_wt>=20)&(fsdf.peaksm_wt==1)].fractionnumbersplicedreadsnotindna.dropna().apply(lambda x: 1 if x>0.01 else 0).sum()/float(fsdf[(fsdf.gfpm_wt<25)&(fsdf.numberreadsm_wt>=20)&(fsdf.peaksm_wt==1)].gfpm_wt.dropna().shape[0])
# 0.006688963210702341 = 36/5382
fsdf[((fsdf.donorstrength>0)|(fsdf.acceptorstrength>0))&(fsdf.numberreadsp_wt>=20)&(fsdf.peaksp_wt==1)].fractionnumbersplicedreadsnotindna.dropna().apply(lambda x: 1 if x>0.01 else 0).sum()/float(fsdf[(fsdf.gfpp_wt<25)&(fsdf.numberreadsp_wt>=20)&(fsdf.peaksp_wt==1)].gfpp_wt.dropna().shape[0])
# 0.0070339976553341153


f=plt.figure(figsize=(4,4))
plt.scatter(fsdf[(fsdf.gfpm_wt<25)&(fsdf.numberreadsm_wt>=20)&(fsdf.peaksm_wt==1)].gfpm_wt, 
            fsdf[(fsdf.gfpm_wt<25)&(fsdf.numberreadsm_wt>=20)&(fsdf.peaksm_wt==1)].fractionnumbersplicedreadsnotindna, alpha=0.2)
plt.scatter(fsdf[((fsdf.donorstrength>0)|(fsdf.acceptorstrength>0))&(fsdf.gfpm_wt<25)&(fsdf.numberreadsm_wt>=20)&(fsdf.peaksm_wt==1)].gfpm_wt, 
            fsdf[((fsdf.donorstrength>0)|(fsdf.acceptorstrength>0))&(fsdf.gfpm_wt<25)&(fsdf.numberreadsm_wt>=20)&(fsdf.peaksm_wt==1)].fractionnumbersplicedreadsnotindna, 
                 color=sns.xkcd_rgb['red'], alpha=0.3)
plt.axvspan(-0.5, 1.3, color=sns.xkcd_rgb['medium blue'], alpha=0.2)
plt.xlim(-0.1,5.1)
plt.ylim(-0.03,1.03)
f.savefig('./figures/gfpm_wt_vs_fractionnumbersplicedreads_max5_withssstrength.png', \
      dpi = 300, format='png', bbox_inches='tight', frameon=True)


f=plt.figure(figsize=(4,4))
plt.scatter(fsdf[(fsdf.gfpp_wt<25)&(fsdf.numberreadsp_wt>=20)&(fsdf.peaksp_wt==1)].gfpp_wt, 
            fsdf[(fsdf.gfpp_wt<25)&(fsdf.numberreadsp_wt>=20)&(fsdf.peaksp_wt==1)].fractionnumbersplicedreadsnotindna, alpha=0.1)
plt.scatter(fsdf[((fsdf.donorstrength>0)|(fsdf.acceptorstrength>0))&(fsdf.gfpp_wt<25)&(fsdf.numberreadsp_wt>=20)&(fsdf.peaksp_wt==1)].gfpp_wt, 
            fsdf[((fsdf.donorstrength>0)|(fsdf.acceptorstrength>0))&(fsdf.gfpp_wt<25)&(fsdf.numberreadsp_wt>=20)&(fsdf.peaksp_wt==1)].fractionnumbersplicedreadsnotindna,
                 color=sns.xkcd_rgb['red'], alpha=0.3)
plt.axvspan(-0.5, 1.3, color=sns.xkcd_rgb['medium blue'], alpha=0.2)
plt.xlim(-0.1,5.1)
plt.ylim(-0.03,1.03)
f.savefig('./figures/gfpp_wt_vs_fractionnumbersplicedreads_max5_withssstrength.png', \
      dpi = 300, format='png', bbox_inches='tight', frameon=True)



#%%

with open('./downstreamsequences.fa', 'w') as f:
    for i in fsdf[fsdf.subset=='fs_downstrvar'].index:
        f.write('>'+str(i)+'\n')
        f.write(str(fsdf.varseq162[i])+'\n')


with open('./downstreamsequences_onlydownstr.fa', 'w') as f:
    for i in fsdf[fsdf.subset=='fs_downstrvar'].index:
        f.write('>'+str(i)+'\n')
        f.write(str(fsdf.varseq162[i][42:])+'\n')


