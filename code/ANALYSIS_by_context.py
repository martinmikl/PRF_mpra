#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 20 10:28:09 2019

@author: martinm
"""



import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
import seaborn as sns
import scipy.stats
import RNA

sns.set_context('talk')


ratiocalcbins=pd.read_pickle('./mapping/minus1/ratiocalcbins1.pkl') 

xval=ratiocalcbins['median'].apply(lambda x: np.log2(x))

xvalstop=pd.read_pickle('./mapping/stopvector/xvals_corrected.pkl') 

# combined dataframes

fsdf=pd.read_pickle('./dataframes/fsdf.pkl')
fsvarsdf=pd.read_pickle('./dataframes/fsvarsdf.pkl')

# new mapping - FS,1

covm1comb=pd.read_pickle('/net/mraid08/export/genie/Runs/Martin/fs_rev/Sample_FS1_v0_s/fsm1_mutscombined_difflibupdown_min3reads_minlength10norm_m1variants.pkl')

covm1vars=pd.read_pickle('/net/mraid08/export/genie/Runs/Martin/fs_rev/Sample_FS1_v0_s/fsm1_difflibupdown_min3reads_minlength10norm_m1variants.pkl')

# new mapping - FSp1
covp1comb=pd.read_pickle('/net/mraid08/export/genie/Runs/Martin/fs_rev/Sample_FS_v0_s/fsp1_mutscombined_difflibupdown_min3reads_minlength10norm_p1variants.pkl')

covp1vars=pd.read_pickle('/net/mraid08/export/genie/Runs/Martin/fs_rev/Sample_FS_v0_s/fsp1_difflibupdown_min3reads_minlength10norm_p1variants.pkl')

# new mapping - stopFSm1

covstopvariants=pd.read_pickle('/net/mraid08/export/genie/Runs/Martin/fs_rev/Sample_stopFS1_v0_s/fscov_FSstop_mutscombined_bincoverage_min3_difflibupdown_minlength10.pkl')


# wt values

wtvalues=pd.read_pickle('./dataframes/wtvalues.pkl')

# define background cutoff

thrshld=np.percentile(fsdf[fsdf.numberreadsstop_wt>=20].wavstop_wt.dropna(), 95) # 8.2

eventlist=list(wtvalues[wtvalues.wavm>thrshld].index.values)                     
plist=list(wtvalues[wtvalues.wavp>thrshld].index.values)                     
nofslist=list(wtvalues[(wtvalues.wavm<thrshld)&(wtvalues.wavp<thrshld)].index.values)                     
  
thrshldgfp=np.percentile(fsdf[fsdf.numberreadsstop_wt>=20].gfpstop_wt.dropna(), 95) # 1.296

                    
fsdf['shifting']=fsdf.wavm_wt.apply(lambda x: 1 if x>thrshld else 0)
fsdf['shifting_p1']=fsdf.wavp_wt.apply(lambda x: 1 if x>thrshld else 0)

#%%
eventcolor={'HIV HXB2':'blue',
 'human - HERV-K10':'cyan',
 'simian srv1 retrovirus':'green',
 'SIVmac239':'purple',
 'PLRV luteovirus':'orange',
 'SARS coronavirus':'red',
 'human T-lymphotropic virus':'olive'}
  

#%% ANALYSIS BY CONTEXT

#event='HIV HXB2'
#event='human - HERV-K10'
os.mkdir('./figures/byevent')
eventfolder='./figures/byevent/'

f=plt.figure(figsize=(3,2))
fsdf[(((fsdf.fsevent=='human - OAZ1')&(fsdf.endoslippery==True))|(fsdf.slippery=='oaz'))&(fsdf.peaksp_wt==1)&(fsdf.numberreadsp_wt>=20)&(fsdf.gfpp_wt<35)&(fsdf.laststopinframep1<13)].gfpp_wt.hist(bins=50, linewidth=0)
plt.axvspan(0,1.3, color=sns.xkcd_rgb['medium blue'], alpha=0.2)
f.savefig('./figures/byevent/oaz_gfpp_hist.png', \
              dpi = 300, format='png', bbox_inches='tight', frameon=True)

f=plt.figure(figsize=(3,2))
fsdf[(((fsdf.fsevent=='human - OAZ1')&(fsdf.endoslippery==True))|(fsdf.slippery=='oaz'))&(fsdf.peaksm_wt==1)&(fsdf.numberreadsm_wt>=20)&(fsdf.gfpm_wt<35)&(fsdf.laststopinframem1<13)].gfpm_wt.hist(bins=30, linewidth=0)
plt.axvspan(0,1.3, color=sns.xkcd_rgb['medium blue'], alpha=0.2)
f.savefig('./figures/byevent/oaz_gfpm_hist.png', \
              dpi = 300, format='png', bbox_inches='tight', frameon=True)



f=plt.figure(figsize=(3,2))
fsdf[[x in eventlist for x in fsdf.fsevent]&(fsdf.endoslippery==True)&(fsdf.peaksp_wt==1)&
     (fsdf.numberreadsp_wt>=20)&(fsdf.gfpp_wt<10)].gfpp_wt.hist(bins=50, linewidth=0)
plt.axvspan(0,1.3, color=sns.xkcd_rgb['medium blue'], alpha=0.2)
f.savefig('./figures/byevent/eventlist_gfpp_wt_hist.png', \
              dpi = 300, format='png', bbox_inches='tight', frameon=True)

f=plt.figure(figsize=(3,2))
fsdf[[x in eventlist for x in fsdf.fsevent]&(fsdf.endoslippery==True)&(fsdf.peaksm_wt==1)&
     (fsdf.numberreadsm_wt>=20)&(fsdf.gfpm_wt<10)&(fsdf.laststopinframem1<13)].gfpm_wt.hist(bins=50, linewidth=0)
plt.axvspan(0,1.3, color=sns.xkcd_rgb['medium blue'], alpha=0.2)
plt.xlim(0,10)
f.savefig('./figures/byevent/eventlist_gfpm_wt_hist.png', \
              dpi = 300, format='png', bbox_inches='tight', frameon=True)


f=plt.figure(figsize=(3,2))
fsdf[[x in plist for x in fsdf.fsevent]&(fsdf.endoslippery==True)&(fsdf.peaksp_wt==1)&
     (fsdf.numberreadsp_wt>=20)&(fsdf.gfpp_wt<10)&(fsdf.laststopinframep1<13)].gfpp_wt.hist(bins=50, linewidth=0)
plt.axvspan(0,1.3, color=sns.xkcd_rgb['medium blue'], alpha=0.2)
f.savefig('./figures/byevent/plist_gfpp_wt_hist.png', \
              dpi = 300, format='png', bbox_inches='tight', frameon=True)

f=plt.figure(figsize=(3,2))
fsdf[[x in plist for x in fsdf.fsevent]&(fsdf.endoslippery==True)&(fsdf.peaksm_wt==1)&
     (fsdf.numberreadsm_wt>=20)&(fsdf.gfpm_wt<10)&(fsdf.laststopinframem1<13)].gfpm_wt.hist(bins=50, linewidth=0)
plt.axvspan(0,1.3, color=sns.xkcd_rgb['medium blue'], alpha=0.2)
plt.xlim(0,10)
f.savefig('./figures/byevent/plist_gfpm_wt_hist.png', \
              dpi = 300, format='png', bbox_inches='tight', frameon=True)


f=plt.figure(figsize=(3,2))
fsdf[[x in nofslist for x in fsdf.fsevent]&(fsdf.endoslippery==True)&(fsdf.peaksp_wt==1)&
     (fsdf.numberreadsp_wt>=20)&(fsdf.gfpp_wt<10)].gfpp_wt.hist(bins=50, linewidth=0)
plt.axvspan(0,1.3, color=sns.xkcd_rgb['medium blue'], alpha=0.2)
f.savefig('./figures/byevent/nofslist_gfpp_wt_hist.png', \
              dpi = 300, format='png', bbox_inches='tight', frameon=True)

f=plt.figure(figsize=(3,2))
fsdf[[x in nofslist for x in fsdf.fsevent]&(fsdf.endoslippery==True)&(fsdf.peaksm_wt==1)&
     (fsdf.numberreadsm_wt>=20)&(fsdf.gfpm_wt<10)&(fsdf.laststopinframem1<13)].gfpm_wt.hist(bins=50, linewidth=0)
plt.axvspan(0,1.3, color=sns.xkcd_rgb['medium blue'], alpha=0.2)
plt.xlim(0,10)
f.savefig('./figures/byevent/nofslist_gfpm_wt_hist.png', \
              dpi = 300, format='png', bbox_inches='tight', frameon=True)


######################


df=fsdf[(fsdf.laststopinframem1<12)&(fsdf.peaksm_wt==1)&(fsdf.numberreadsm_wt>=20)&(fsdf.gfpm_wt<25)&[x in eventlist+nofslist+plist for x in fsdf.fsevent]]

wtseqs=pd.Series()

for event in fsdf[fsdf.subset=='fsmain'].fsevent.unique():
    wtseqs.loc[event]=fsdf[(fsdf.subset=='fsmain')&(fsdf.fsevent==event)].varseq162.values[0]
        
df['endoslip']=df.index.map(lambda x: bool(df.varseq162[x][33:42]==wtseqs[df.fsevent[x]][33:42]))
df['endoup']=df.index.map(lambda x: bool(df.varseq162[x][:33]==wtseqs[df.fsevent[x]][:33]))
df['endodown']=df.index.map(lambda x: bool(df.varseq162[x][42:]==wtseqs[df.fsevent[x]][42:]))

    
overview_normgfp=pd.DataFrame()
for i in eventlist:
    overview_normgfp.loc[i,'upstream']=df[(df.fsevent==i)&
            (df.endoup==False)&(df.endoslip==True)&(df.endodown==True)].percgfpm_wt.dropna().median()
    overview_normgfp.loc[i,'slippery']=df[(df.fsevent==i)&
            (df.endoup==True)&(df.endoslip==False)&(df.endodown==True)].percgfpm_wt.dropna().median()
    overview_normgfp.loc[i,'downstream']=df[(df.fsevent==i)&
            (df.endoup==True)&(df.endoslip==True)&(df.endodown==False)].percgfpm_wt.dropna().median()

cg=sns.clustermap(data=overview_normgfp.dropna(), figsize=(3,4), annot=True, 
                  cbar_kws={'ticks':[-100,0,100]},
               annot_kws={'fontsize':12}, vmin=0, fmt='.0f', col_cluster=False)
plt.setp(cg.ax_heatmap.yaxis.get_majorticklabels(), rotation=0)
cg.savefig('./figures/byevent/clustermap_overview_percgfpm_median_eventlist_peaks1minreads20.png', \
              dpi = 300, format='png', bbox_inches='tight', frameon=True)


### Plot mean difference per position for each event

for event in fsdf[(fsdf.subset=='fsmain')].fsevent.unique():    
    df=fsdf[(fsdf.fsevent==event)&(fsdf.numberreadsm_wt>=10)&(fsdf.peaksm_wt==1)&(fsdf.gfpm_wt<25)]
    df=df.dropna(axis=1, how='all')
    endo=fsdf[(fsdf.fsevent==event)&(fsdf.subset=='fsmain')].varseq162.unique()[0]

    ident=df.varseq162.apply(lambda x: pd.Series([x[i]==endo[i] for i in range(len(x))]))

    meandiff=pd.Series()
    for pos in ident.columns:
        meandiff.loc[pos]=df[ident[pos]==False].gfpm_wt.mean()-df[ident[pos]==True].gfpm_wt.mean()
    
    meandiffp=pd.Series()
    for pos in ident.columns:
        meandiffp.loc[pos]=scipy.stats.mannwhitneyu(df[ident[pos]==True].gfpm_wt.dropna(),df[ident[pos]==False].gfpm_wt.dropna())[1]
    
    f=plt.figure(figsize=(6,3))
    ax=meandiff.plot(color=sns.xkcd_rgb['light blue'])
    meandiff.rolling(5, center=True).median().plot(color=sns.xkcd_rgb['medium blue'])
    plt.axvspan(35,42,alpha=0.2)
    plt.xlabel('position along the RNA')
    plt.ylabel('difference in mean\nGFP fluorescence')
    plt.axhline(y=0, linewidth=2, c='gray', alpha=0.6)
    f.savefig('./figures/byevent/identicaltoendogenous_perposition_diffinmeanratio_gfpm_wt_'+'_'.join(event.split(' ')) +'.png', \
              dpi = 300, format='png', bbox_inches='tight', frameon=True)
    
    
    f=plt.figure(figsize=(6,3))
    ax=np.log10(meandiffp).plot(color=sns.xkcd_rgb['light blue'])
    np.log10(meandiffp).rolling(5, center=True).median().plot(color=sns.xkcd_rgb['medium blue'])
    plt.axvspan(35,42,alpha=0.2)
    plt.xlabel('position along the RNA')
    plt.ylabel('Mann-Whitney U\np-value [log10]')
    f.savefig('./figures/byevent/identicaltoendogenous_perposition_diffinmeanratio_gfpm_wt_pval_'+'_'.join(event.split(' ')) +'.png', \
              dpi = 300, format='png', bbox_inches='tight', frameon=True)
    


for event in fsdf[(fsdf.subset=='fsmain')].fsevent.unique():    
    df=fsdf[(fsdf.fsevent==event)&(fsdf.numberreadsp_wt>=10)&(fsdf.peaksp_wt==1)&(fsdf.gfpp_wt<35)]
    df=df.dropna(axis=1, how='all')
    endo=fsdf[(fsdf.fsevent==event)&(fsdf.subset=='fsmain')].varseq162.unique()[0]

    ident=df.varseq162.apply(lambda x: pd.Series([x[i]==endo[i] for i in range(len(x))]))

    meandiff=pd.Series()
    for pos in ident.columns:
        meandiff.loc[pos]=df[ident[pos]==False].gfpp_wt.mean()-df[ident[pos]==True].gfpp_wt.mean()
    
    meandiffp=pd.Series()
    for pos in ident.columns:
        meandiffp.loc[pos]=scipy.stats.mannwhitneyu(df[ident[pos]==True].gfpp_wt.dropna(),df[ident[pos]==False].gfpp_wt.dropna())[1]
    
    f=plt.figure(figsize=(6,3))
    ax=meandiff.plot(color=sns.xkcd_rgb['light blue'])
    meandiff.rolling(5, center=True).median().plot(color=sns.xkcd_rgb['medium blue'])
    plt.axvspan(35,42,alpha=0.2)
    plt.xlabel('position along the RNA')
    plt.ylabel('difference in mean\nGFP fluorescence')
    plt.axhline(y=0, linewidth=2, c='gray', alpha=0.6)
    f.savefig('./figures/byevent/identicaltoendogenous_perposition_diffinmeanratio_gfpp_wt_'+'_'.join(event.split(' ')) +'.png', \
              dpi = 300, format='png', bbox_inches='tight', frameon=True)
    
    
    f=plt.figure(figsize=(6,3))
    ax=np.log10(meandiffp).plot(color=sns.xkcd_rgb['light blue'])
    np.log10(meandiffp).rolling(5, center=True).median().plot(color=sns.xkcd_rgb['medium blue'])
    plt.axvspan(35,42,alpha=0.2)
    plt.xlabel('position along the RNA')
    plt.ylabel('Mann-Whitney U\np-value [log10]')
    f.savefig('./figures/byevent/identicaltoendogenous_perposition_diffinmeanratio_gfpp_wt_pval_'+'_'.join(event.split(' ')) +'.png', \
              dpi = 300, format='png', bbox_inches='tight', frameon=True)
    

# with wild-type slippery site

for event in fsdf[(fsdf.subset=='fsmain')].fsevent.unique():    
    wtslippery=fsdf[(fsdf.subset=='fsmain')&(fsdf.fsevent==event)].varseq162.values[0][35:42]
    df=fsdf[(fsdf.fsevent==event)&(fsdf.numberreadsm_wt>=10)&(fsdf.peaksm_wt==1)&(fsdf.gfpm_wt<25)&[x[35:42]==wtslippery for x in fsdf.varseq162]]
    df=df.dropna(axis=1, how='all')
    endo=fsdf[(fsdf.fsevent==event)&(fsdf.subset=='fsmain')].varseq162.unique()[0]

    ident=df.varseq162.apply(lambda x: pd.Series([x[i]==endo[i] for i in range(len(x))]))

    meandiff=pd.Series()
    for pos in ident.columns:
        meandiff.loc[pos]=df[ident[pos]==False].gfpm_wt.mean()-df[ident[pos]==True].gfpm_wt.mean()
    
    meandiffp=pd.Series()
    for pos in ident.columns:
        meandiffp.loc[pos]=scipy.stats.mannwhitneyu(df[ident[pos]==True].gfpm_wt.dropna(),df[ident[pos]==False].gfpm_wt.dropna())[1]
    
    f=plt.figure(figsize=(6,3))
    ax=meandiff.plot(color=sns.xkcd_rgb['light blue'])
    meandiff.rolling(5, center=True).median().plot(color=sns.xkcd_rgb['medium blue'])
    plt.axvspan(35,42,alpha=0.2)
    plt.xlabel('position along the RNA')
    plt.ylabel('difference in mean\nGFP fluorescence')
    plt.axhline(y=0, linewidth=2, c='gray', alpha=0.6)
    f.savefig('./figures/byevent/identicaltoendogenous_perposition_diffinmeanratio_gfpm_wt_'+'_'.join(event.split(' ')) +'.png', \
              dpi = 300, format='png', bbox_inches='tight', frameon=True)
    
    
    f=plt.figure(figsize=(6,3))
    ax=np.log10(meandiffp).plot(color=sns.xkcd_rgb['light blue'])
    np.log10(meandiffp).rolling(5, center=True).median().plot(color=sns.xkcd_rgb['medium blue'])
    plt.axvspan(35,42,alpha=0.2)
    plt.xlabel('position along the RNA')
    plt.ylabel('Mann-Whitney U\np-value [log10]')
    f.savefig('./figures/byevent/identicaltoendogenous_perposition_endoslippery_diffinmeanratio_gfpm_wt_pval_'+'_'.join(event.split(' ')) +'.png', \
              dpi = 300, format='png', bbox_inches='tight', frameon=True)
    


for event in fsdf[(fsdf.subset=='fsmain')].fsevent.unique():    
    wtslippery=fsdf[(fsdf.subset=='fsmain')&(fsdf.fsevent==event)].varseq162.values[0][35:42]
    df=fsdf[(fsdf.fsevent==event)&(fsdf.numberreadsp_wt>=10)&(fsdf.peaksp_wt==1)&(fsdf.gfpp_wt<35)&[x[35:42]==wtslippery for x in fsdf.varseq162]]
    df=df.dropna(axis=1, how='all')
    endo=fsdf[(fsdf.fsevent==event)&(fsdf.subset=='fsmain')].varseq162.unique()[0]

    ident=df.varseq162.apply(lambda x: pd.Series([x[i]==endo[i] for i in range(len(x))]))

    meandiff=pd.Series()
    for pos in ident.columns:
        meandiff.loc[pos]=df[ident[pos]==False].gfpp_wt.mean()-df[ident[pos]==True].gfpp_wt.mean()
    
    meandiffp=pd.Series()
    for pos in ident.columns:
        meandiffp.loc[pos]=scipy.stats.mannwhitneyu(df[ident[pos]==True].gfpp_wt.dropna(),df[ident[pos]==False].gfpp_wt.dropna())[1]
    
    f=plt.figure(figsize=(6,3))
    ax=meandiff.plot(color=sns.xkcd_rgb['light blue'])
    meandiff.rolling(5, center=True).median().plot(color=sns.xkcd_rgb['medium blue'])
    plt.axvspan(35,42,alpha=0.2)
    plt.xlabel('position along the RNA')
    plt.ylabel('difference in mean\nGFP fluorescence')
    plt.axhline(y=0, linewidth=2, c='gray', alpha=0.6)
    f.savefig('./figures/byevent/identicaltoendogenous_perposition_diffinmeanratio_gfpp_wt_'+'_'.join(event.split(' ')) +'.png', \
              dpi = 300, format='png', bbox_inches='tight', frameon=True)
    
    
    f=plt.figure(figsize=(6,3))
    ax=np.log10(meandiffp).plot(color=sns.xkcd_rgb['light blue'])
    np.log10(meandiffp).rolling(5, center=True).median().plot(color=sns.xkcd_rgb['medium blue'])
    plt.axvspan(35,42,alpha=0.2)
    plt.xlabel('position along the RNA')
    plt.ylabel('Mann-Whitney U\np-value [log10]')
    f.savefig('./figures/byevent/identicaltoendogenous_perposition_endoslippery_diffinmeanratio_gfpp_wt_pval_'+'_'.join(event.split(' ')) +'.png', \
              dpi = 300, format='png', bbox_inches='tight', frameon=True)
    
    
    
# All PRF events together - normalized PRF rates   
    
dfs=pd.DataFrame()
endo=pd.Series()
for event in eventlist:   
    df=fsdf[(fsdf.fsevent==event)&(fsdf.peaksm_wt==1)&(fsdf.numberreadsm_wt>=10)&(fsdf.laststopinframem1<13)&(fsdf.gfpm_wt<25)]
    df=df.dropna(axis=1, how='all')
    endo.loc[event]=fsdf[(fsdf.fsevent==event)&(fsdf.subset=='fsmain')].varseq162.unique()[0]
    dfs=dfs.append(df)
    
dfs['libindex']=dfs.index.map(lambda x: int(x))

ident=dfs.libindex.apply(lambda x: pd.Series([dfs.varseq162[x][i]==endo.loc[dfs.fsevent[x]][i] for i in range(162)]))

meandiff=pd.Series()
for pos in ident.columns:
    meandiff.loc[pos]=dfs[ident[pos]==False].gfpm_wt.mean()-dfs[ident[pos]==True].gfpm_wt.mean()

meandiffp=pd.Series()
for pos in ident.columns:
    meandiffp.loc[pos]=scipy.stats.mannwhitneyu(dfs[ident[pos]==True].gfpm_wt.dropna(),dfs[ident[pos]==False].gfpm_wt.dropna())[1]

f=plt.figure(figsize=(6,3))
ax=meandiff.plot(color=sns.xkcd_rgb['light blue'])
meandiff.rolling(5, center=True).median().plot(color=sns.xkcd_rgb['medium blue'])
plt.axvspan(35,42,alpha=0.2)
plt.xlabel('position along the RNA')
plt.ylabel('difference in mean\nGFP fluorescence')
plt.axhline(y=0, linewidth=2, c='gray', alpha=0.6)
f.savefig('./figures/byevent/identicaltoendogenous_perposition_diffinmeanratio_eventlist.png', \
          dpi = 300, format='png', bbox_inches='tight', frameon=True)


f=plt.figure(figsize=(6,3))
ax=np.log10(meandiffp).plot(color=sns.xkcd_rgb['light blue'])
np.log10(meandiffp).rolling(5, center=True).median().plot(color=sns.xkcd_rgb['medium blue'])
plt.axvspan(35,42,alpha=0.2)
plt.xlabel('position along the RNA')
plt.ylabel('Mann-Whitney U\np-value [log10]')
f.savefig('./figures/byevent/identicaltoendogenous_perposition_diffinmeanratio_pval_eventlist.png', \
          dpi = 300, format='png', bbox_inches='tight', frameon=True)



# All PRF events together - percent PRF rates

dfs=pd.DataFrame()
endo=pd.Series()
for event in eventlist:   
    df=fsdf[(fsdf.fsevent==event)&(fsdf.peaksm_wt==1)&(fsdf.numberreadsm_wt>=10)&(fsdf.laststopinframem1<13)&(fsdf.gfpm_wt<25)]
    df=df.dropna(axis=1, how='all')
    endo.loc[event]=fsdf[(fsdf.fsevent==event)&(fsdf.subset=='fsmain')].varseq162.unique()[0]
    dfs=dfs.append(df)
    
dfs['libindex']=dfs.index.map(lambda x: int(x))

ident=dfs.libindex.apply(lambda x: pd.Series([dfs.varseq162[x][i]==endo.loc[dfs.fsevent[x]][i] for i in range(162)]))

meandiff=pd.Series()
for pos in ident.columns:
    meandiff.loc[pos]=dfs[ident[pos]==False].percgfpm_wt.mean()-dfs[ident[pos]==True].percgfpm_wt.mean()

meandiffp=pd.Series()
for pos in ident.columns:
    meandiffp.loc[pos]=scipy.stats.mannwhitneyu(dfs[ident[pos]==True].percgfpm_wt.dropna(),dfs[ident[pos]==False].percgfpm_wt.dropna())[1]

f=plt.figure(figsize=(6,3))
ax=meandiff.plot(color=sns.xkcd_rgb['light blue'])
meandiff.rolling(5, center=True).median().plot(color=sns.xkcd_rgb['medium blue'])
plt.axvspan(35,42,alpha=0.2)
plt.xlabel('position along the RNA')
plt.ylabel('difference in\n% GFP fluorescence')
plt.axhline(y=0, linewidth=2, c='gray', alpha=0.6)
f.savefig('./figures/byevent/identicaltoendogenous_perposition_diffinmeanperc_eventlist.png', \
          dpi = 300, format='png', bbox_inches='tight', frameon=True)


f=plt.figure(figsize=(6,3))
ax=np.log10(meandiffp).plot(color=sns.xkcd_rgb['light blue'])
np.log10(meandiffp).rolling(5, center=True).median().plot(color=sns.xkcd_rgb['medium blue'])
plt.axvspan(35,42,alpha=0.2)
plt.xlabel('position along the RNA')
plt.ylabel('Mann-Whitney U\np-value [log10]')
f.savefig('./figures/byevent/identicaltoendogenous_perposition_diffinmeanperc_pval_eventlist.png', \
          dpi = 300, format='png', bbox_inches='tight', frameon=True)



# For PRF events tested negative in our assay
    
dfs=pd.DataFrame()
endo=pd.Series()
for event in nofslist:   
    df=fsdf[(fsdf.fsevent==event)&(fsdf.peaksm_wt==1)&(fsdf.numberreadsm_wt>=10)&(fsdf.laststopinframem1<13)&(fsdf.gfpm_wt<25)]
    df=df.dropna(axis=1, how='all')
    endo.loc[event]=fsdf[(fsdf.fsevent==event)&(fsdf.subset=='fsmain')].varseq162.unique()[0]
    dfs=dfs.append(df)
    
dfs['libindex']=dfs.index.map(lambda x: int(x))

ident=dfs.libindex.apply(lambda x: pd.Series([dfs.varseq162[x][i]==endo.loc[dfs.fsevent[x]][i] for i in range(162)]))

meandiff=pd.Series()
for pos in ident.columns:
    meandiff.loc[pos]=dfs[ident[pos]==False].gfpm_wt.mean()-dfs[ident[pos]==True].gfpm_wt.mean()

meandiffp=pd.Series()
for pos in ident.columns:
    meandiffp.loc[pos]=scipy.stats.mannwhitneyu(dfs[ident[pos]==True].gfpm_wt.dropna(),dfs[ident[pos]==False].gfpm_wt.dropna())[1]

f=plt.figure(figsize=(6,3))
ax=meandiff.plot(color=sns.xkcd_rgb['light blue'])
meandiff.rolling(5, center=True).median().plot(color=sns.xkcd_rgb['medium blue'])
plt.axvspan(35,42,alpha=0.2)
plt.xlabel('position along the RNA')
plt.ylabel('difference in mean\nGFP fluorescence')
plt.axhline(y=0, linewidth=2, c='gray', alpha=0.6)
f.savefig('./figures/byevent/identicaltoendogenous_perposition_diffinmeanratio_nofslist.png', \
          dpi = 300, format='png', bbox_inches='tight', frameon=True)


f=plt.figure(figsize=(6,3))
ax=np.log10(meandiffp).plot(color=sns.xkcd_rgb['light blue'])
np.log10(meandiffp).rolling(5, center=True).median().plot(color=sns.xkcd_rgb['medium blue'])
plt.axvspan(35,42,alpha=0.2)
plt.xlabel('position along the RNA')
plt.ylabel('Mann-Whitney U\np-value [log10]')
f.savefig('./figures/byevent/identicaltoendogenous_perposition_diffinmeanratio_pval_nofslist.png', \
          dpi = 300, format='png', bbox_inches='tight', frameon=True)


#### plot together


df=pd.DataFrame()
endo=pd.Series()
for event in eventlist:   
    dfs=fsdf[(fsdf.fsevent==event)&(fsdf.peaksm_wt==1)&(fsdf.numberreadsm_wt>=10)&(fsdf.laststopinframem1<13)&(fsdf.gfpm_wt<25)]
    endo.loc[event]=fsdf[(fsdf.fsevent==event)&(fsdf.subset=='fsmain')].varseq162.unique()[0]
    df=df.append(dfs)
    
df['libindex']=df.index.map(lambda x: int(x))

ident=df.libindex.apply(lambda x: pd.Series([df.varseq162[x][i]==endo.loc[df.fsevent[x]][i] for i in range(162)]))

eventcolor={'HIV HXB2':'blue',
 'human - HERV-K10':'cyan',
 'simian srv1 retrovirus':'green',
 'SIVmac239':'purple',
 'PLRV luteovirus':'orange',
 'SARS coronavirus':'red',
 'human T-lymphotropic virus':'olive'}
 
nofscolor={'human - PEG10':'blue',
 'human - CCR5':'red',
 'herpes simplex - thymidine kinase':'green',
 'SIVmac239':'purple',
 'Rous sarcoma virus':'orange'}
 

f=plt.figure(figsize=(6,3))
#for event in ['HIV HXB2','SARS coronavirus','SIVmac239']:
for event in eventcolor.keys():
    if event!='human T-lymphotropic virus':
        meandiff=pd.Series()
        for pos in ident.columns:
            meandiff.loc[pos]=df[(ident[pos]==False)&(df.fsevent==event)].percgfpm_wt.median()-df[(ident[pos]==True)&(df.fsevent==event)].percgfpm_wt.median()
        
     #   ax=meandiff.plot(color=sns.xkcd_rgb['light blue'])
        meandiff.rolling(5, center=True).median().plot(color=eventcolor[event])
        plt.axvspan(35,42,alpha=0.05)
        plt.ylim(-50,100)
        plt.xlabel('position along the RNA')
        plt.ylabel('difference in mean\nGFP fluorescence')
        plt.axhline(y=0, linewidth=2, c='gray', alpha=0.6)
f.savefig('./figures/byevent/identicaltoendogenous_perposition_diffinmeanratio_eventlist_plottedtogether.png', \
              dpi = 300, format='png', bbox_inches='tight', frameon=True)
    
f=plt.figure(figsize=(6,3))
for event in eventcolor.keys():
    if event!='human T-lymphotropic virus':
        meandiffp=pd.Series()
        for pos in ident.columns:
            meandiffp.loc[pos]=scipy.stats.mannwhitneyu(df[(ident[pos]==False)&(df.fsevent==event)].percgfpm_wt.dropna(),df[(ident[pos]==True)&(df.fsevent==event)].percgfpm_wt.dropna())[1]
        
    #    ax=np.log10(meandiffp).plot(color=sns.xkcd_rgb['light blue'])
        np.log10(meandiffp).rolling(5, center=True).median().plot(color=eventcolor[event])
        plt.axvspan(35,42,alpha=0.05)
        plt.xlabel('position along the RNA')
        plt.ylabel('Mann-Whitney U\np-value [log10]')
f.savefig('./figures/byevent/identicaltoendogenous_perposition_diffinmeanratio_eventlist_plottedtogether_pvals.png', \
                  dpi = 300, format='png', bbox_inches='tight', frameon=True)
    


f=plt.figure(figsize=(6,3))
#for event in ['HIV HXB2','SARS coronavirus','SIVmac239']:
for event in eventcolor.keys():
    if event!='human T-lymphotropic virus':
        meandiff=pd.Series()
        for pos in ident.columns:
            meandiff.loc[pos]=df[(ident[pos]==False)&(df.fsevent==event)].percgfpm_wt.mean()
        
     #   ax=meandiff.plot(color=sns.xkcd_rgb['light blue'])
        meandiff.rolling(5, center=True).median().plot(color=eventcolor[event])
        plt.axvspan(35,42,alpha=0.05)
        plt.ylim(0,130)
        plt.xlabel('position along the RNA')
        plt.ylabel('mean % wild-type\nGFP fluorescence')
        plt.axhline(y=100, linewidth=2, c='gray', alpha=0.6)
f.savefig('./figures/byevent/identicaltoendogenous_perposition_meanpercgfpwhenmutated_eventlist_plottedtogether.png', \
              dpi = 300, format='png', bbox_inches='tight', frameon=True)

f=plt.figure(figsize=(6,3))
#for event in ['HIV HXB2','SARS coronavirus','SIVmac239']:
for event in eventcolor.keys():
    if event!='human T-lymphotropic virus':
        meandiff=pd.Series()
        for pos in ident.columns:
            meandiff.loc[pos]=df[(ident[pos]==False)&(df.fsevent==event)].percgfpm_wt.median()
        
     #   ax=meandiff.plot(color=sns.xkcd_rgb['light blue'])
        meandiff.rolling(5, center=True).median().plot(color=eventcolor[event])
        plt.axvspan(35,42,alpha=0.05)
        plt.ylim(0,130)
        plt.xlabel('position along the RNA')
        plt.ylabel('median % wild-type\nGFP fluorescence')
        plt.axhline(y=100, linewidth=2, c='gray', alpha=0.6)
f.savefig('./figures/byevent/identicaltoendogenous_perposition_medianpercgfpwhenmutated_eventlist_plottedtogether.png', \
              dpi = 300, format='png', bbox_inches='tight', frameon=True)

    
f=plt.figure(figsize=(6,3))
for event in eventcolor.keys():
    if event!='human T-lymphotropic virus':
        meandiffp=pd.Series()
        for pos in ident.columns:
            meandiffp.loc[pos]=scipy.stats.wilcoxon(df[(ident[pos]==False)&(df.fsevent==event)].percgfpm_wt.dropna()-100)[1]
        
    #    ax=np.log10(meandiffp).plot(color=sns.xkcd_rgb['light blue'])
        np.log10(meandiffp).rolling(5, center=True).median().plot(color=eventcolor[event])
        plt.axvspan(35,42,alpha=0.05)
        plt.xlabel('position along the RNA')
        plt.ylabel('Wilcoxon signed-rank test\np-value [log10]')
f.savefig('./figures/byevent/identicaltoendogenous_perposition_meanpercgfpwhenmutated_eventlist_plottedtogether_pvals.png', \
                  dpi = 300, format='png', bbox_inches='tight', frameon=True)
    

# Plot together for PRF events tested negative in our assay

df=pd.DataFrame()
endo=pd.Series()
for event in nofslist:   
    dfs=fsdf[(fsdf.fsevent==event)&(fsdf.peaksm_wt==1)&(fsdf.numberreadsm_wt>=10)&(fsdf.laststopinframem1<13)&(fsdf.gfpm_wt<25)]
    endo.loc[event]=fsdf[(fsdf.fsevent==event)&(fsdf.subset=='fsmain')].varseq162.unique()[0]
    df=df.append(dfs)
    
df['libindex']=df.index.map(lambda x: int(x))

ident=df.libindex.apply(lambda x: pd.Series([df.varseq162[x][i]==endo.loc[df.fsevent[x]][i] for i in range(162)]))

f=plt.figure(figsize=(6,3))
#for event in ['HIV HXB2','SARS coronavirus','SIVmac239']:
for event in nofslist:
    meandiff=pd.Series()
    for pos in ident.columns:
        meandiff.loc[pos]=df[(ident[pos]==False)&(df.fsevent==event)].percgfpm_wt.mean()
    
 #   ax=meandiff.plot(color=sns.xkcd_rgb['light blue'])
    meandiff.rolling(5, center=True).median().plot(color=nofscolor[event])
    plt.axvspan(35,42,alpha=0.05)
    plt.xlabel('position along the RNA')
    plt.ylabel('mean % wild-type\nGFP fluorescence')
    plt.axhline(y=100, linewidth=2, c='gray', alpha=0.6)
f.savefig('./figures/byevent/identicaltoendogenous_perposition_meanpercgfpwhenmutated_nofslist_plottedtogether.png', \
              dpi = 300, format='png', bbox_inches='tight', frameon=True)

f=plt.figure(figsize=(6,3))
#for event in ['HIV HXB2','SARS coronavirus','SIVmac239']:
for event in nofslist:
    meandiff=pd.Series()
    for pos in ident.columns:
        meandiff.loc[pos]=df[(ident[pos]==False)&(df.fsevent==event)].percgfpm_wt.median()
    
 #   ax=meandiff.plot(color=sns.xkcd_rgb['light blue'])
    meandiff.rolling(5, center=True).median().plot(color=nofscolor[event])
    plt.axvspan(35,42,alpha=0.05)
    plt.xlabel('position along the RNA')
    plt.ylabel('median % wild-type\nGFP fluorescence')
    plt.axhline(y=100, linewidth=2, c='gray', alpha=0.6)
f.savefig('./figures/byevent/identicaltoendogenous_perposition_medianpercgfpwhenmutated_nofslist_plottedtogether.png', \
              dpi = 300, format='png', bbox_inches='tight', frameon=True)

    
f=plt.figure(figsize=(6,3))
for event in nofslist:
    meandiffp=pd.Series()
    for pos in ident.columns:
        meandiffp.loc[pos]=scipy.stats.wilcoxon(df[(ident[pos]==False)&(df.fsevent==event)].percgfpm_wt.dropna()-100)[1]
    
#    ax=np.log10(meandiffp).plot(color=sns.xkcd_rgb['light blue'])
    np.log10(meandiffp).rolling(5, center=True).median().plot(color=nofscolor[event])
    plt.axvspan(35,42,alpha=0.05)
    plt.xlabel('position along the RNA')
    plt.ylabel('Wilcoxon signed-rank test\np-value [log10]')
f.savefig('./figures/byevent/identicaltoendogenous_perposition_meanpercgfpwhenmutated_nofslist_plottedtogether_pvals.png', \
                  dpi = 300, format='png', bbox_inches='tight', frameon=True)
    


###    
# Plot histograms for variants with mutations in upstream, slippery, or downstream region separately

for event in fsdf[(fsdf.subset=='fsmain')].fsevent.unique():    
    df=fsdf[(fsdf.fsevent==event)&(fsdf.peaksm_wt==1)&(fsdf.numberreadsm_wt>=20)&(fsdf.laststopinframem1<13)&(fsdf.gfpm_wt<25)]
    df=df.dropna(axis=1, how='all')
    endo=fsdf[(fsdf.fsevent==event)&(fsdf.subset=='fsmain')].varseq162.unique()[0]
    
    df['upstream']=df.varseq162.apply(lambda x: x[:33])
    df['slippery']=df.varseq162.apply(lambda x: x[33:42])
    df['downstream']=df.varseq162.apply(lambda x: x[42:])
    maxval=np.ceil(df[(df.gfpm_wt<25)].gfpm_wt.max())
    
    f=plt.figure(figsize=(4,3))
    df[(df.gfpm_wt<25)].gfpm_wt.hist(bins=50)
    plt.axvline(x=wtvalues.loc[event,'gfpm'], linewidth=2, color=sns.xkcd_rgb['dark red'])
    plt.xlim(0,maxval)
    plt.axvspan(0,1.3,alpha=0.2)
    plt.xlabel('% GFP expression')
    plt.ylabel('counts')
    f.savefig('./figures/byevent/hist_allvars_'+'_'.join(event.split(' ')) +'.png', \
              dpi = 300, format='png', bbox_inches='tight', frameon=True)
    
    f=plt.figure(figsize=(4,3))
    df[(df.gfpm_wt<25)&(df.slippery!=endo[33:42])&(df.downstream==endo[42:])&(df.upstream==endo[:33])].gfpm_wt.hist(bins=50)
    plt.axvline(x=wtvalues.loc[event,'gfpm'], linewidth=2, color=sns.xkcd_rgb['dark red'], alpha=0.5)
    plt.xlim(0,maxval)
    plt.axvspan(0,1.3,alpha=0.2)
    plt.xlabel('% GFP expression')
    plt.ylabel('counts')
    f.savefig('./figures/byevent/hist_mutslippery_'+'_'.join(event.split(' ')) +'.png', \
              dpi = 300, format='png', bbox_inches='tight', frameon=True)
    
    f=plt.figure(figsize=(4,3))
    df[(df.gfpm_wt<25)&(df.peaksm_wt==1)&(df.numberreadsm_wt>=20)&(df.downstream!=endo[42:])&(df.slippery==endo[33:42])&(df.upstream==endo[:33])].gfpm_wt.hist(bins=50)
    plt.axvline(x=wtvalues.loc[event,'gfpm'], linewidth=2, color=sns.xkcd_rgb['dark red'], alpha=0.5)
    plt.xlim(0,maxval)
    plt.axvspan(0,1.3,alpha=0.2)
    plt.xlabel('% GFP expression')
    plt.ylabel('counts')
    f.savefig('./figures/byevent/hist_mutdownstream_'+'_'.join(event.split(' ')) +'.png', \
              dpi = 300, format='png', bbox_inches='tight', frameon=True)
    
    f=plt.figure(figsize=(4,3))
    df[(df.gfpm_wt<25)&(df.peaksm_wt==1)&(df.numberreadsm_wt>=20)&(df.upstream!=endo[:33])&(df.downstream==endo[42:])&(df.slippery==endo[33:42])].gfpm_wt.hist(bins=50)
    plt.axvline(x=wtvalues.loc[event,'gfpm'], linewidth=2, color=sns.xkcd_rgb['dark red'], alpha=0.5)
    plt.xlim(0,maxval)
    plt.axvspan(0,1.3,alpha=0.2)
    plt.xlabel('% GFP expression')
    plt.ylabel('counts')
    f.savefig('./figures/byevent/hist_mutupstream_'+'_'.join(event.split(' ')) +'.png', \
              dpi = 300, format='png', bbox_inches='tight', frameon=True)

for event in fsdf[(fsdf.subset=='fsmain')].fsevent.unique():    
    df=fsdf[(fsdf.fsevent==event)&(fsdf.peaksm_wt==1)&(fsdf.numberreadsm_wt>=20)&(fsdf.laststopinframem1<13)&(fsdf.gfpm_wt<25)]
    df=df.dropna(axis=1, how='all')
    endo=fsdf[(fsdf.fsevent==event)&(fsdf.subset=='fsmain')].varseq162.unique()[0]
    
    df['upstream']=df.varseq162.apply(lambda x: x[:33])
    df['slippery']=df.varseq162.apply(lambda x: x[33:42])
    df['downstream']=df.varseq162.apply(lambda x: x[42:])
    print
    print event
    print scipy.stats.wilcoxon(df[(df.gfpm_wt<25)&(df.peaksm_wt==1)&(df.numberreadsm_wt>=20)&
                            (df.upstream!=endo[:33])&(df.downstream==endo[42:])&(df.slippery==endo[33:42])].gfpm_wt-wtvalues.loc[event,'gfpm'])
    print scipy.stats.wilcoxon(df[(df.gfpm_wt<25)&(df.peaksm_wt==1)&(df.numberreadsm_wt>=20)&
                            (df.slippery!=endo[33:42])&(df.downstream==endo[42:])&(df.upstream==endo[:33])].gfpm_wt-wtvalues.loc[event,'gfpm'])
    print scipy.stats.wilcoxon(df[(df.gfpm_wt<25)&(df.peaksm_wt==1)&(df.numberreadsm_wt>=20)&
                            (df.upstream==endo[:33])&(df.downstream!=endo[42:])&(df.slippery==endo[33:42])].gfpm_wt-wtvalues.loc[event,'gfpm'])

'''
results: upstreammut, slipperymut, downstreammut

human - CCR5
WilcoxonResult(statistic=40.0, pvalue=0.70070341225168198)
WilcoxonResult(statistic=100.0, pvalue=0.85192459481464233)
WilcoxonResult(statistic=188.0, pvalue=0.15528204380540125)

SARS coronavirus
WilcoxonResult(statistic=319.0, pvalue=0.32178087953675816)
WilcoxonResult(statistic=137.0, pvalue=1.1363463258731764e-12)
WilcoxonResult(statistic=227.0, pvalue=1.3432608890718677e-24)

SIVmac239
WilcoxonResult(statistic=788.0, pvalue=0.93498798684581397)
WilcoxonResult(statistic=290.0, pvalue=9.7744503453816537e-05)
WilcoxonResult(statistic=701.0, pvalue=1.969521892383602e-13)

simian srv1 retrovirus
WilcoxonResult(statistic=172.0, pvalue=0.0023486394938546282)
WilcoxonResult(statistic=0.0, pvalue=1.7106552178285207e-15)
WilcoxonResult(statistic=42.0, pvalue=7.8306892144285228e-26)

PLRV luteovirus
WilcoxonResult(statistic=1006.0, pvalue=0.30752503023611466)
WilcoxonResult(statistic=759.0, pvalue=2.4235271204097168e-13)
WilcoxonResult(statistic=1747.0, pvalue=3.3851171828298306e-13)

human - OAZ1
WilcoxonResult(statistic=0.0, pvalue=nan)
WilcoxonResult(statistic=0.0, pvalue=nan)
WilcoxonResult(statistic=95.0, pvalue=2.8996750933304209e-11)

human - PEG10
WilcoxonResult(statistic=91.0, pvalue=0.0020933562871036549)
WilcoxonResult(statistic=1199.0, pvalue=0.02017785002091806)
WilcoxonResult(statistic=2533.0, pvalue=6.1442295609633017e-05)

Rous sarcoma virus
WilcoxonResult(statistic=1734.0, pvalue=0.96740654285245498)
WilcoxonResult(statistic=1119.0, pvalue=0.00076613126385520632)
WilcoxonResult(statistic=2428.0, pvalue=0.0078494213555923755)

HCV - F protein
WilcoxonResult(statistic=0.0, pvalue=nan)
WilcoxonResult(statistic=0.0, pvalue=nan)
WilcoxonResult(statistic=229.0, pvalue=0.51281362625336624)

human - HERV-K10
WilcoxonResult(statistic=434.0, pvalue=0.00039871776839621703)
WilcoxonResult(statistic=6.0, pvalue=2.1227659050061804e-14)
WilcoxonResult(statistic=253.0, pvalue=2.21014836195947e-24)

influenza a virus
WilcoxonResult(statistic=0.0, pvalue=nan)
WilcoxonResult(statistic=0.0, pvalue=nan)
WilcoxonResult(statistic=152.0, pvalue=0.0075896370375709635)

human T-lymphotropic virus
WilcoxonResult(statistic=389.0, pvalue=0.43452083346948633)
WilcoxonResult(statistic=268.0, pvalue=0.82326477762683115)
WilcoxonResult(statistic=1139.0, pvalue=0.24493780831422007)

HIV HXB2
WilcoxonResult(statistic=856.0, pvalue=0.21850960321275303)
WilcoxonResult(statistic=38.0, pvalue=1.0874474301356334e-21)
WilcoxonResult(statistic=1398.0, pvalue=8.5928844585964466e-13)

west nile virus
WilcoxonResult(statistic=36.0, pvalue=0.17284834913542624)
WilcoxonResult(statistic=140.0, pvalue=0.36725440950812072)
WilcoxonResult(statistic=392.0, pvalue=8.7070957352798552e-14)
'''

for event in fsdf[(fsdf.subset=='fsmain')].fsevent.unique():    
    df=fsdf[(fsdf.fsevent==event)&(fsdf.peaksm_wt==1)&(fsdf.numberreadsm_wt>=20)&(fsdf.laststopinframem1<13)&(fsdf.gfpm_wt<25)]
    df=df.dropna(axis=1, how='all')
    endo=fsdf[(fsdf.fsevent==event)&(fsdf.subset=='fsmain')].varseq162.unique()[0]
    
    df['upstream']=df.varseq162.apply(lambda x: x[:33])
    df['slippery']=df.varseq162.apply(lambda x: x[33:42])
    df['downstream']=df.varseq162.apply(lambda x: x[42:])
    
    f=plt.figure(figsize=(4,3))
    df[(df.gfpm_wt<25)].gfpm_wt.hist(bins=50)
    plt.xlabel('% GFP expression')
    plt.ylabel('counts')
    f.savefig('./figures/byevent/hist_allvars_'+'_'.join(event.split(' ')) +'.png', \
              dpi = 300, format='png', bbox_inches='tight', frameon=True)
    
    f=plt.figure(figsize=(4,3))
    df[(df.gfpm_wt<25)&(df.slippery==endo[33:42])].gfpm_wt.hist(bins=50)
    plt.xlabel('% GFP expression')
    plt.ylabel('counts')
    f.savefig('./figures/byevent/hist_endoslippery_'+'_'.join(event.split(' ')) +'.png', \
              dpi = 300, format='png', bbox_inches='tight', frameon=True)
    
    f=plt.figure(figsize=(4,3))
    df[(df.gfpm_wt<25)&(df.downstream==endo[42:])].gfpm_wt.hist(bins=50)
    plt.xlabel('% GFP expression')
    plt.ylabel('counts')
    f.savefig('./figures/byevent/hist_endodownstream_'+'_'.join(event.split(' ')) +'.png', \
              dpi = 300, format='png', bbox_inches='tight', frameon=True)
    
    f=plt.figure(figsize=(4,3))
    df[(df.gfpm_wt<25)&(df.upstream==endo[:33])].gfpm_wt.hist(bins=50)
    plt.xlabel('% GFP expression')
    plt.ylabel('counts')
    f.savefig('./figures/byevent/hist_endoupstream_'+'_'.join(event.split(' ')) +'.png', \
              dpi = 300, format='png', bbox_inches='tight', frameon=True)
        
    sec=df.varseq.apply(lambda x: pd.Series(list(RNA.fold(x)[0])))
    sec=sec.replace({'.':False,'(':True,')':True})
    
    meandiff=pd.Series()
    for pos in sec.columns[30:-18]:
        meandiff.loc[pos-30]=df[sec[pos]==True].gfpm_wt.mean()-df[sec[pos]==False].gfpm_wt.mean()
    
    
    meandiffp=pd.Series()
    for pos in sec.columns[30:-18]:
        meandiffp.loc[pos-30]=scipy.stats.mannwhitneyu(df[sec[pos]==True].gfpm_wt.dropna(),df[sec[pos]==False].gfpm_wt.dropna())[1]
    
    f=plt.figure(figsize=(6,3))
    ax=meandiff.plot(color=sns.xkcd_rgb['light blue'])
    meandiff.rolling(5, center=True).median().plot(color=sns.xkcd_rgb['medium blue'])
    plt.axvspan(35,42,alpha=0.2)
    plt.xlabel('position along the RNA')
    plt.ylabel('difference in mean\nGFP fluorescence')
    plt.axhline(y=0, linewidth=2, c='gray', alpha=0.6)
    f.savefig('./figures/byevent/paired_perposition_diffinmeanratio_'+'_'.join(event.split(' ')) +'.png', \
              dpi = 300, format='png', bbox_inches='tight', frameon=True)
    
    
    f=plt.figure(figsize=(6,3))
    ax=np.log10(meandiffp).plot(color=sns.xkcd_rgb['light blue'])
    np.log10(meandiffp).rolling(5, center=True).median().plot(color=sns.xkcd_rgb['medium blue'])
    plt.axvspan(35,42,alpha=0.2)
    plt.xlabel('position along the RNA')
    plt.ylabel('Mann-Whitney U\np-value [log10]')
    f.savefig('./figures/byevent/paired_perposition_diffinmeanratio_pval_'+'_'.join(event.split(' ')) +'.png', \
              dpi = 300, format='png', bbox_inches='tight', frameon=True)
    
    
    

for event in fsdf[(fsdf.subset=='fsmain')].fsevent.unique():   
    df=fsdf[(fsdf.fsevent==event)&(fsdf.endoslippery==True)&(fsdf.peaksm_wt==1)&(fsdf.numberreadsm_wt>=10)&(fsdf.laststopinframem1<13)&(fsdf.gfpm_wt<25)]
    df=df.dropna(axis=1, how='all')
    
    df['upstream']=df.varseq162.apply(lambda x: x[:33])
    df['slippery']=df.varseq162.apply(lambda x: x[33:42])
    df['downstream']=df.varseq162.apply(lambda x: x[42:])
    
    
    sec=df.varseq.apply(lambda x: pd.Series(list(RNA.fold(x)[0])))
    sec=sec.replace({'.':False,'(':True,')':True})
    
    
    meandiff=pd.Series()
    for pos in sec.columns[30:-18]:
        meandiff.loc[pos-30]=df[sec[pos]==True].gfpm_wt.mean()-df[sec[pos]==False].gfpm_wt.mean()
    
    
    meandiffp=pd.Series()
    for pos in sec.columns[30:-18]:
        meandiffp.loc[pos-30]=scipy.stats.mannwhitneyu(df[sec[pos]==True].gfpm_wt.dropna(),df[sec[pos]==False].gfpm_wt.dropna())[1]
    
    f=plt.figure(figsize=(6,3))
    ax=meandiff.plot(color=sns.xkcd_rgb['light blue'])
    meandiff.rolling(5, center=True).median().plot(color=sns.xkcd_rgb['medium blue'])
    plt.axvspan(35,42,alpha=0.2)
    plt.xlabel('position along the RNA')
    plt.ylabel('difference in mean\nGFP fluorescence')
    plt.axhline(y=0, linewidth=2, c='gray', alpha=0.6)
    f.savefig('./figures/byevent/paired_perposition_slipperyendog_diffinmeanratio_'+'_'.join(event.split(' ')) +'.png', \
              dpi = 300, format='png', bbox_inches='tight', frameon=True)
    
    
    f=plt.figure(figsize=(6,3))
    ax=np.log10(meandiffp).plot(color=sns.xkcd_rgb['light blue'])
    np.log10(meandiffp).rolling(5, center=True).median().plot(color=sns.xkcd_rgb['medium blue'])
    plt.axvspan(35,42,alpha=0.2)
    plt.xlabel('position along the RNA')
    plt.ylabel('Mann-Whitney U\np-value [log10]')
    f.savefig('./figures/byevent/paired_perposition_slipperyendog_diffinmeanratio_pval_'+'_'.join(event.split(' ')) +'.png', \
              dpi = 300, format='png', bbox_inches='tight', frameon=True)
    



    
for event in fsdf[(fsdf.subset=='fsmain')].fsevent.unique():    
    df=fsdf[(fsdf.fsevent==event)&(fsdf.peaksp_wt==1)&(fsdf.numberreadsp_wt>=20)&(fsdf.laststopinframep1<13)&(fsdf.gfpp_wt<25)]
    df=df.dropna(axis=1, how='all')
    endo=fsdf[(fsdf.fsevent==event)&(fsdf.subset=='fsmain')].varseq162.unique()[0]
    
    df['upstream']=df.varseq162.apply(lambda x: x[:33])
    df['slippery']=df.varseq162.apply(lambda x: x[33:42])
    df['downstream']=df.varseq162.apply(lambda x: x[42:])
    
    f=plt.figure(figsize=(4,3))
    df[(df.gfpp_wt<25)].gfpp_wt.hist(bins=50)
    plt.xlabel('% GFP expression')
    plt.ylabel('counts')
    f.savefig('./figures/byevent/hist_allvars_gfpp_wt_'+'_'.join(event.split(' ')) +'.png', \
              dpi = 300, format='png', bbox_inches='tight', frameon=True)
    
    f=plt.figure(figsize=(4,3))
    df[(df.gfpp_wt<25)&(df.slippery==endo[33:42])].gfpp_wt.hist(bins=50)
    plt.xlabel('% GFP expression')
    plt.ylabel('counts')
    f.savefig('./figures/byevent/hist_endoslippery_gfpp_wt_'+'_'.join(event.split(' ')) +'.png', \
              dpi = 300, format='png', bbox_inches='tight', frameon=True)
    
    f=plt.figure(figsize=(4,3))
    df[(df.gfpp_wt<25)&(df.downstream==endo[42:])].gfpp_wt.hist(bins=50)
    plt.xlabel('% GFP expression')
    plt.ylabel('counts')
    f.savefig('./figures/byevent/hist_endodownstream_gfpp_wt_'+'_'.join(event.split(' ')) +'.png', \
              dpi = 300, format='png', bbox_inches='tight', frameon=True)
    
    f=plt.figure(figsize=(4,3))
    df[(df.gfpp_wt<25)&(df.upstream==endo[:33])].gfpp_wt.hist(bins=50)
    plt.xlabel('% GFP expression')
    plt.ylabel('counts')
    f.savefig('./figures/byevent/hist_endoupstream_gfpp_wt_'+'_'.join(event.split(' ')) +'.png', \
              dpi = 300, format='png', bbox_inches='tight', frameon=True)
        
    sec=df.varseq.apply(lambda x: pd.Series(list(RNA.fold(x)[0])))
    sec=sec.replace({'.':False,'(':True,')':True})
    
    meandiff=pd.Series()
    for pos in sec.columns[30:-18]:
        meandiff.loc[pos-30]=df[sec[pos]==True].gfpp_wt.mean()-df[sec[pos]==False].gfpp_wt.mean()
    
    
    meandiffp=pd.Series()
    for pos in sec.columns[30:-18]:
        meandiffp.loc[pos-30]=scipy.stats.mannwhitneyu(df[sec[pos]==True].gfpp_wt.dropna(),df[sec[pos]==False].gfpp_wt.dropna())[1]
    
    f=plt.figure(figsize=(6,3))
    ax=meandiff.plot(color=sns.xkcd_rgb['light blue'])
    meandiff.rolling(5, center=True).median().plot(color=sns.xkcd_rgb['medium blue'])
    plt.axvspan(35,42,alpha=0.2)
    plt.xlabel('position along the RNA')
    plt.ylabel('difference in mean\nGFP fluorescence')
    plt.axhline(y=0, linewidth=2, c='gray', alpha=0.6)
    f.savefig('./figures/byevent/paired_perposition_diffinmeanratio_gfpp_wt_'+'_'.join(event.split(' ')) +'.png', \
              dpi = 300, format='png', bbox_inches='tight', frameon=True)
    
    
    f=plt.figure(figsize=(6,3))
    ax=np.log10(meandiffp).plot(color=sns.xkcd_rgb['light blue'])
    np.log10(meandiffp).rolling(5, center=True).median().plot(color=sns.xkcd_rgb['medium blue'])
    plt.axvspan(35,42,alpha=0.2)
    plt.xlabel('position along the RNA')
    plt.ylabel('Mann-Whitney U\np-value [log10]')
    f.savefig('./figures/byevent/paired_perposition_diffinmeanratio_gfpp_wt_pval_'+'_'.join(event.split(' ')) +'.png', \
              dpi = 300, format='png', bbox_inches='tight', frameon=True)
    
    
    

### plot preferences for secondary structure together

dfs=fsdf[(fsdf.gfpm_wt<25)&[x in eventlist for x in fsdf.fsevent]&(fsdf.endoslippery==True)]
dfs=dfs.dropna(axis=1, how='all')

sec=dfs.varseq.apply(lambda x: pd.Series(list(RNA.fold(x)[0])))
sec=sec.replace({'.':False,'(':True,')':True})

f=plt.figure(figsize=(6,3))
for event in ['HIV HXB2','SARS coronavirus','SIVmac239']:
    meandiff=pd.Series()
    df=fsdf[(fsdf.gfpm_wt<25)&(fsdf.fsevent==event)&(fsdf.endoslippery==True)]
    for pos in sec.columns[30:-18]:
        meandiff.loc[pos-30]=df[(sec[pos]==True)&(df.fsevent==event)].percgfpm_wt.mean()-df[(sec[pos]==False)&(df.fsevent==event)].percgfpm_wt.mean()
    
 #   ax=meandiff.plot(color=sns.xkcd_rgb['light blue'])
    meandiff.rolling(5, center=True).median().plot(color=eventcolor[event])
plt.axvspan(35,42,alpha=0.2)
plt.ylim(-50,50)
plt.xlabel('position along the RNA')
plt.ylabel('difference in mean\nGFP fluorescence')
plt.axhline(y=0, linewidth=2, c='gray', alpha=0.6)
f.savefig('./figures/byevent/paired_perposition_slipperyendog_diffinmeanratio_hiv_sars_siv.png', \
          dpi = 300, format='png', bbox_inches='tight', frameon=True)


f=plt.figure(figsize=(6,3))
for event in ['HIV HXB2','SARS coronavirus','simian srv1 retrovirus', 'human - HERV-K10']:
    meandiff=pd.Series()
    df=fsdf[(fsdf.gfpm_wt<25)&(fsdf.fsevent==event)&(fsdf.endoslippery==True)]
    for pos in sec.columns[30:-18]:
        meandiff.loc[pos-30]=df[(sec[pos]==True)&(df.fsevent==event)].percgfpm_wt.mean()-df[(sec[pos]==False)&(df.fsevent==event)].percgfpm_wt.mean()
    
 #   ax=meandiff.plot(color=sns.xkcd_rgb['light blue'])
    meandiff.rolling(5, center=True).median().plot(color=eventcolor[event])
plt.axvspan(35,42,alpha=0.2)
plt.ylim(-50,50)
plt.xlabel('position along the RNA')
plt.ylabel('difference in mean\nGFP fluorescence')
plt.axhline(y=0, linewidth=2, c='gray', alpha=0.6)
f.savefig('./figures/byevent/paired_perposition_slipperyendog_diffinmeanratio_hiv_sars_srv1_herv.png', \
          dpi = 300, format='png', bbox_inches='tight', frameon=True)

f=plt.figure(figsize=(6,3))
for event in ['HIV HXB2','SARS coronavirus']:
    meandiff=pd.Series()
    df=fsdf[(fsdf.gfpm_wt<25)&(fsdf.fsevent==event)&(fsdf.endoslippery==True)]
    for pos in sec.columns[30:-18]:
        meandiff.loc[pos-30]=df[(sec[pos]==True)&(df.fsevent==event)].percgfpm_wt.mean()-df[(sec[pos]==False)&(df.fsevent==event)].percgfpm_wt.mean()
    
 #   ax=meandiff.plot(color=sns.xkcd_rgb['light blue'])
    meandiff.rolling(5, center=True).median().plot(color=eventcolor[event])
plt.axvspan(35,42,alpha=0.2)
plt.ylim(-50,50)
plt.xlabel('position along the RNA')
plt.ylabel('difference in mean\nGFP fluorescence')
plt.axhline(y=0, linewidth=2, c='gray', alpha=0.6)
f.savefig('./figures/byevent/paired_perposition_slipperyendog_diffinmeanratio_hiv_sars.png', \
          dpi = 300, format='png', bbox_inches='tight', frameon=True)

f=plt.figure(figsize=(6,3))
for event in ['simian srv1 retrovirus', 'human - HERV-K10']:
    meandiff=pd.Series()
    df=fsdf[(fsdf.gfpm_wt<25)&(fsdf.fsevent==event)&(fsdf.endoslippery==True)]
    for pos in sec.columns[30:-18]:
        meandiff.loc[pos-30]=df[(sec[pos]==True)&(df.fsevent==event)].percgfpm_wt.mean()-df[(sec[pos]==False)&(df.fsevent==event)].percgfpm_wt.mean()
    
 #   ax=meandiff.plot(color=sns.xkcd_rgb['light blue'])
    meandiff.rolling(5, center=True).median().plot(color=eventcolor[event])
plt.axvspan(35,42,alpha=0.2)
plt.ylim(-50,50)
plt.xlabel('position along the RNA')
plt.ylabel('difference in mean\nGFP fluorescence')
plt.axhline(y=0, linewidth=2, c='gray', alpha=0.6)
f.savefig('./figures/byevent/paired_perposition_slipperyendog_diffinmeanratio_srv1_herv.png', \
          dpi = 300, format='png', bbox_inches='tight', frameon=True)


    
f=plt.figure(figsize=(6,3))
for event in ['HIV HXB2','SARS coronavirus','SIVmac239']:
    meandiffp=pd.Series()
    for pos in sec.columns[30:-18]:
        meandiffp.loc[pos-30]=scipy.stats.mannwhitneyu(df[(sec[pos]==True)&(df.fsevent==event)].percgfpm_wt.dropna(),df[(sec[pos]==False)&(df.fsevent==event)].percgfpm_wt.dropna())[1]
    
#    ax=np.log10(meandiffp).plot(color=sns.xkcd_rgb['light blue'])
    np.log10(meandiffp).rolling(5, center=True).median().plot()
    plt.axvspan(35,42,alpha=0.2)
    plt.xlabel('position along the RNA')
    plt.ylabel('Mann-Whitney U\np-value [log10]')
f.savefig('./figures/byevent/paired_perposition_slipperyendog_diffinmeanratio_pval_hiv_sars_siv.png', \
              dpi = 300, format='png', bbox_inches='tight', frameon=True)
