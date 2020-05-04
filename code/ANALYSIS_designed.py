#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 14 15:01:01 2019

@author: martinm
"""


import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from Bio.Seq import Seq
from scipy.stats import pearsonr
import scipy.stats
import difflib

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


#%%
eventcolor={'HIV HXB2':'blue',
 'human - HERV-K10':'cyan',
 'simian srv1 retrovirus':'green',
 'SIVmac239':'purple',
 'PLRV luteovirus':'orange',
 'SARS coronavirus':'red',
 'human T-lymphotropic virus':'olive'}
    #%% read in codon usage table into dictionary


codonusagedict={}

with open ('../additional/CodonUsageHuman_from_kazusa_dot_org.txt','r') as infile:
    table = infile.read()
    table_list = table.splitlines()
    
for row in range (0, len(table_list)):
    linearr=str.split(table_list[row], '\t')
    if len(str(linearr[1]))==3:
        codonusagedict[linearr[1]]=float(linearr[3])


#%% read in tAI table into dictionary


tAIdict={}

with open ('../additional/tAI_index_human_nar-02315.txt') as infile:
    table = infile.read()
    table_list = table.splitlines()
    
for row in range (0, len(table_list)):
    linearr=str.split(table_list[row], '\t')
    if len(str(linearr[2]))==3:
        if linearr[4]!=' ----' and linearr[4]!='----':
            tAIdict[linearr[2]]=float(linearr[4])
        else:
            tAIdict[linearr[2]]=float(0)

genomic_copynumber={}

for row in range (0, len(table_list)):
    linearr=str.split(table_list[row], '\t')
    if len(str(linearr[2]))==3:
        genomic_copynumber[linearr[2]]=int(float(linearr[3]))
        


        
#%% Plot multiple barcodes for same variant

wtseqs=pd.Series()

for event in fsdf[fsdf.subset=='fsmain'].fsevent.unique():
    wtseqs.loc[event]=fsdf[(fsdf.subset=='fsmain')&(fsdf.fsevent==event)].varseq162.values[0]

indexes=pd.Series()
samevarseq=[]
for var in eventlist:
    samevarseq.append(var)
    indexes.append(list(fsdf[fsdf.varseq162==var].index))


for event in eventlist:
    f=plt.figure(figsize=(3,2))
    indexes=list(fsdf[(fsdf.varseq162==wtseqs.loc[event])&(fsdf.peaksm_wt==1)&(fsdf.numberreadsm_wt>=20)&(fsdf.gfpm_wt<25)].index)
    for i in indexes:
        plt.plot(xval, covm1comb.loc[str(i)+'_wt'], color=sns.xkcd_rgb['medium blue'], alpha=0.2)
        plt.title(event)
    plt.xlim(7,np.ceil(np.max(xval)))
    f.savefig('./figures/multiplebarcodes_peaks1minreads20_'+('_').join(event.split(' '))+'.png', \
       dpi = 300, format='png', bbox_inches='tight', frameon=True)

for event in nofslist:
    f=plt.figure(figsize=(3,2))
    indexes=list(fsdf[(fsdf.varseq162==wtseqs.loc[event])&(fsdf.peaksm_wt==1)&(fsdf.numberreadsm_wt>=20)&(fsdf.gfpm_wt<25)].index)
    for i in indexes:
        plt.plot(xval, covm1comb.loc[str(i)+'_wt'], color=sns.xkcd_rgb['medium blue'], alpha=0.2)
        plt.title(event)
    plt.xlim(7,np.ceil(np.max(xval)))
    f.savefig('./figures/multiplebarcodes_peaks1minreads20_'+('_').join(event.split(' '))+'.png', \
       dpi = 300, format='png', bbox_inches='tight', frameon=True)

for event in plist:
    f=plt.figure(figsize=(3,2))
    indexes=list(fsdf[(fsdf.varseq162==wtseqs.loc[event])&(fsdf.peaksm_wt==1)&(fsdf.numberreadsm_wt>=20)&(fsdf.gfpm_wt<25)].index)
    for i in indexes:
        plt.plot(xval, covm1comb.loc[str(i)+'_wt'], color=sns.xkcd_rgb['medium blue'], alpha=0.2)
        plt.title(event)
    plt.xlim(7,np.ceil(np.max(xval)))
    f.savefig('./figures/multiplebarcodes_m1_peaks1minreads20_'+('_').join(event.split(' '))+'.png', \
       dpi = 300, format='png', bbox_inches='tight', frameon=True)

'''
for event in eventlist:
    f=plt.figure(figsize=(3,2))
    indexes=list(fsdf[(fsdf.varseq162==wtseqs.loc[event])&(fsdf.peaksstop_wt==1)&(fsdf.numberreadsstop_wt>=20)].index)
    for i in indexes:
        plt.plot(xvalstop, covstopvariants.loc[str(i)+'.0_wt'], color=sns.xkcd_rgb['medium blue'], alpha=0.2)
    plt.title(event)
    plt.xlim(7,np.ceil(np.max(xval)))
    f.savefig('./figures/multiplebarcodes_peaks1minreads20_'+('_').join(event.split(' '))+'.png', \
       dpi = 300, format='png', bbox_inches='tight', frameon=True)
'''   


for event in plist:
    f=plt.figure(figsize=(3,2))
    indexes=list(fsdf[(fsdf.varseq162==wtseqs.loc[event])&(fsdf.peaksp_wt==1)&(fsdf.numberreadsp_wt>=20)&(fsdf.gfpp_wt<25)].index)
    for i in indexes:
        plt.plot(xval, covp1comb.loc[str(i)+'_wt'], color=sns.xkcd_rgb['medium blue'], alpha=0.2)
        plt.title(event)
    plt.xlim(7,np.ceil(np.max(xval)))
    f.savefig('./figures/multiplebarcodes_p_peaks1minreads20_'+('_').join(event.split(' '))+'.png', \
       dpi = 300, format='png', bbox_inches='tight', frameon=True)


for event in eventlist:
    f=plt.figure(figsize=(3,2))
    indexes=list(fsdf[(fsdf.varseq162==wtseqs.loc[event])&(fsdf.peaksp_wt==1)&(fsdf.numberreadsp_wt>=20)&(fsdf.gfpp_wt<25)].index)
    for i in indexes:
        plt.plot(xval, covp1comb.loc[str(i)+'_wt'], color=sns.xkcd_rgb['medium blue'], alpha=0.2)
        plt.title(event)
    plt.xlim(7,np.ceil(np.max(xval)))
    f.savefig('./figures/multiplebarcodes_p_peaks1minreads20_'+('_').join(event.split(' '))+'.png', \
       dpi = 300, format='png', bbox_inches='tight', frameon=True)



# Plot distributions

f=plt.figure(figsize=(4,3))
fsdf[(fsdf.peaksm_wt==1)&(fsdf.numberreadsm_wt>20)&(fsdf.laststopinframem1<12)].wavm_wt.hist(bins=100, linewidth=0)
plt.xlabel('GFP expression')
plt.ylabel('# variants')
f.savefig('./figures/hist_wavm_wt_fullibrary.png', \
   dpi = 300, format='png', bbox_inches='tight', frameon=True)

f=plt.figure(figsize=(4,3))
fsdf[(fsdf.peaksm_del==1)&(fsdf.numberreadsm_del>20)].wavm_del.hist(bins=100, linewidth=0)
plt.xlabel('GFP expression')
plt.ylabel('# variants')
f.savefig('./figures/hist_wavm_del_fullibrary.png', \
   dpi = 300, format='png', bbox_inches='tight', frameon=True)

f=plt.figure(figsize=(4,3))
fsdf[(fsdf.peaksm_ins==1)&(fsdf.numberreadsm_ins>20)].wavm_ins.hist(bins=100, linewidth=0)
plt.xlabel('GFP expression')
plt.ylabel('# variants')
f.savefig('./figures/hist_wavm_ins_fullibrary.png', \
   dpi = 300, format='png', bbox_inches='tight', frameon=True)


f=plt.figure(figsize=(4,3))
fsdf[(fsdf.peaksstop_wt==1)&(fsdf.numberreadsstop_wt>20)].wavstop_wt.hist(bins=100, linewidth=0)
plt.xlabel('GFP expression')
plt.ylabel('# variants')
plt.xlim(7,14)
f.savefig('./figures/hist_wavstop_wt_fullibrary.png', \
   dpi = 300, format='png', bbox_inches='tight', frameon=True)



np.percentile(fsdf[(fsdf.peaksm_ins==1)&(fsdf.wavm_ins>10)].wavm_ins.values, 5)
# 12.05 --> 25%

np.percentile(fsdf[(fsdf.peaksp_ins==1)&(fsdf.numberreadsp_ins>20)&(fsdf.wavp_ins>thrshld)].wavp_ins.values, 5)
# 12.34



f=plt.figure(figsize=(4,3))
fsdf[(fsdf.peaksp_wt==1)&(fsdf.numberreadsp_wt>20)&(fsdf.laststopinframep1<12)].wavp_wt.hist(bins=100, linewidth=0)
plt.xlabel('GFP expression')
plt.ylabel('# variants')
f.savefig('./figures/hist_wavp_wt_fullibrary.png', \
   dpi = 300, format='png', bbox_inches='tight', frameon=True)

f=plt.figure(figsize=(4,3))
fsdf[(fsdf.peaksp_del==1)&(fsdf.numberreadsp_del>20)].wavp_del.hist(bins=100, linewidth=0)
plt.xlabel('GFP expression')
plt.ylabel('# variants')
f.savefig('./figures/hist_wavp_del_fullibrary.png', \
   dpi = 300, format='png', bbox_inches='tight', frameon=True)

f=plt.figure(figsize=(4,3))
fsdf[(fsdf.peaksp_ins==1)&(fsdf.numberreadsp_ins>20)].wavp_ins.hist(bins=100, linewidth=0)
plt.xlabel('GFP expression')
plt.ylabel('# variants')
f.savefig('./figures/hist_wavp_ins_fullibrary.png', \
   dpi = 300, format='png', bbox_inches='tight', frameon=True)






        #%%   SLIPPERY
##############################
##############################

df=fsdf[(fsdf.subset=="fs_slipperyvar")]
df=df.dropna(axis=1, how='all')
df['subgroup']=df.group.apply(lambda x: x.split('subgroup ')[1][0])
df['stopmut']=df.group.apply(lambda x: True if 'stop' in x else False)

df1=df[(df.subgroup=='1')&(df.laststopinframem1<12)]
df1['slippery']=df1.varseq162.apply(lambda x: x.upper()[35:42])
df1['slipperyfirst']=df1.slippery.apply(lambda x: x[:3])
df1['slipperysecond']=df1.slippery.apply(lambda x: x[3:6])
df1['slipperylast']=df1.slippery.apply(lambda x: x[-1:])

df1['tai_0first']=df1.slippery.apply(lambda x: tAIdict[x[1:4]])
df1['tai_0second']=df1.slippery.apply(lambda x: tAIdict[x[4:7]])
df1['tai_m1first']=df1.slippery.apply(lambda x: tAIdict[x[0:3]])
df1['tai_m1second']=df1.slippery.apply(lambda x: tAIdict[x[3:6]])

df1['cpy_0first']=df1.slippery.apply(lambda x: genomic_copynumber[x[1:4]])
df1['cpy_0second']=df1.slippery.apply(lambda x: genomic_copynumber[x[4:7]])
df1['cpy_m1first']=df1.slippery.apply(lambda x: genomic_copynumber[x[0:3]])
df1['cpy_m1second']=df1.slippery.apply(lambda x: genomic_copynumber[x[3:6]])

def gettaidiffmap(vrsq):
    tai=[]
    for i in range(len(vrsq)/3-1):
        tai.append(tAIdict[vrsq[2+3*i:5+3*i]]-tAIdict[vrsq[3+3*i:6+3*i]])
    return pd.Series(tai)

taidiffslip=df1.varseq162.apply(lambda x: gettaidiffmap(x.upper()))



############ HEATMAPS FOR SLIPPERY GROUP WISE

# first vs second    
f=plt.figure(figsize=(4,3))
sns.heatmap(data=df1[(df1.numberreadsm_wt>=20)&(df1.wavm_wt<12)&[x in eventlist for x in df1.fsevent]].pivot_table(index='slipperyfirst',
                     columns='slipperysecond',values='gfpm_wt',aggfunc=np.median), 
    annot=df1[(df1.numberreadsm_wt>=20)&(df1.wavm_wt<12)&[x in eventlist for x in df1.fsevent]].pivot_table(index='slipperyfirst',
              columns='slipperysecond',values='gfpm_wt',aggfunc=np.median), fmt='.1f')
f.savefig('./figures/slippery/subgroup1_first_vs_second_median_eventlist_gfp.png', \
   dpi = 300, format='png', bbox_inches='tight', frameon=True)

f=plt.figure(figsize=(4,3))
sns.heatmap(data=df1[(df1.numberreadsm_wt>=20)&(df1.wavm_wt<12)&[x in eventlist for x in df1.fsevent]].pivot_table(index='slipperyfirst',
                     columns='slipperysecond',values='wavm_wt',aggfunc=np.median), 
    annot=df1[(df1.numberreadsm_wt>=20)&(df1.wavm_wt<12)&[x in eventlist for x in df1.fsevent]].pivot_table(index='slipperyfirst',
              columns='slipperysecond',values='wavm_wt',aggfunc=np.median), fmt='.1f')
f.savefig('./figures/slippery/subgroup1_first_vs_second_median_eventlist_wav.png', \
   dpi = 300, format='png', bbox_inches='tight', frameon=True)

f=plt.figure(figsize=(4,3))
sns.heatmap(data=df1[(df1.numberreadsm_wt>=20)&(df1.wavm_wt<12)&[x in eventlist for x in df1.fsevent]].pivot_table(index='slipperyfirst',
                     columns='slipperysecond',values='percgfpm_wt',aggfunc=np.median), 
    annot=df1[(df1.numberreadsm_wt>=20)&(df1.wavm_wt<12)&[x in eventlist for x in df1.fsevent]].pivot_table(index='slipperyfirst',
              columns='slipperysecond',values='percgfpm_wt',aggfunc=np.median), fmt='.1f')
f.savefig('./figures/slippery/subgroup1_first_vs_second_median_eventlist_gfpperc.png', \
   dpi = 300, format='png', bbox_inches='tight', frameon=True)

# second vs last
f=plt.figure(figsize=(4,3))
sns.heatmap(data=df1[(df1.numberreadsm_wt>=20)&(df1.wavm_wt<12)&[x in eventlist for x in df1.fsevent]].pivot_table(index='slipperysecond',
                     columns='slipperylast',values='gfpm_wt',aggfunc=np.median), 
    annot=df1[(df1.numberreadsm_wt>=20)&(df1.wavm_wt<12)&[x in eventlist for x in df1.fsevent]].pivot_table(index='slipperysecond',
              columns='slipperylast',values='gfpm_wt',aggfunc=np.median), fmt='.1f')
f.savefig('./figures/slippery/subgroup1_second_vs_last_median_eventlist_gfp.png', \
   dpi = 300, format='png', bbox_inches='tight', frameon=True)

f=plt.figure(figsize=(4,3))
sns.heatmap(data=df1[(df1.numberreadsm_wt>=20)&(df1.wavm_wt<12)&[x in eventlist for x in df1.fsevent]].pivot_table(index='slipperysecond',
                     columns='slipperylast',values='wavm_wt',aggfunc=np.median), 
    annot=df1[(df1.numberreadsm_wt>=20)&(df1.wavm_wt<12)&[x in eventlist for x in df1.fsevent]].pivot_table(index='slipperysecond',
              columns='slipperylast',values='wavm_wt',aggfunc=np.median), fmt='.1f')
f.savefig('./figures/slippery/subgroup1_second_vs_last_median_eventlist_wav.png', \
   dpi = 300, format='png', bbox_inches='tight', frameon=True)

f=plt.figure(figsize=(4,3))
sns.heatmap(data=df1[(df1.numberreadsm_wt>=20)&(df1.wavm_wt<12)&[x in eventlist for x in df1.fsevent]].pivot_table(index='slipperysecond',
                     columns='slipperylast',values='percgfpm_wt',aggfunc=np.median), 
    annot=df1[(df1.numberreadsm_wt>=20)&(df1.wavm_wt<12)&[x in eventlist for x in df1.fsevent]].pivot_table(index='slipperysecond',
              columns='slipperylast',values='percgfpm_wt',aggfunc=np.median), fmt='.1f')
f.savefig('./figures/slippery/subgroup1_second_vs_last_median_eventlist_gfpperc.png', \
   dpi = 300, format='png', bbox_inches='tight', frameon=True)



f=plt.figure(figsize=(4,3))
sns.heatmap(data=df1[(df1.numberreadsm_wt>=20)&(df1.wavm_wt<12)&[x in nofslist for x in df1.fsevent]].pivot_table(index='slipperyfirst',
                     columns='slipperysecond',values='percgfpm_wt',aggfunc=np.median), 
    annot=df1[(df1.numberreadsm_wt>=20)&(df1.wavm_wt<12)&[x in nofslist for x in df1.fsevent]].pivot_table(index='slipperyfirst',
              columns='slipperysecond',values='percgfpm_wt',aggfunc=np.mean), fmt='.1f')
f.savefig('./figures/slippery/subgroup1_first_vs_second_median_nofslist_gfpperc.png', \
   dpi = 300, format='png', bbox_inches='tight', frameon=True)

f=plt.figure(figsize=(4,3))
sns.heatmap(data=df1[(df1.numberreadsm_wt>=20)&(df1.wavm_wt<12)&[x in nofslist for x in df1.fsevent]].pivot_table(index='slipperysecond',
                     columns='slipperylast',values='percgfpm_wt',aggfunc=np.median), 
    annot=df1[(df1.numberreadsm_wt>=20)&(df1.wavm_wt<12)&[x in nofslist for x in df1.fsevent]].pivot_table(index='slipperysecond',
              columns='slipperylast',values='percgfpm_wt',aggfunc=np.mean), fmt='.1f')
f.savefig('./figures/slippery/subgroup1_second_vs_last_median_nofslist_gfpperc.png', \
   dpi = 300, format='png', bbox_inches='tight', frameon=True)

f=plt.figure(figsize=(4,3))
ax=sns.boxplot(data=df1[[x in eventlist for x in df1.fsevent]&
        (df1.endoslippery==False)&(df1.gfpm_wt<25)&(df1.percgfpm_wt<200)], 
    x='fsevent', y='percgfpm_wt', order=[x for x in eventlist if x not in ['PRRSV - nsp2TF','west nile virus']], color=sns.xkcd_rgb['light blue'])
ax=sns.swarmplot(data=df1[[x in eventlist for x in df1.fsevent]&
        (df1.endoslippery==False)&(df1.gfpm_wt<25)&(df1.percgfpm_wt<200)], 
    x='fsevent', y='percgfpm_wt', order=[x for x in eventlist if x not in ['PRRSV - nsp2TF','west nile virus']], color=sns.xkcd_rgb['medium blue'])
ax.set_xticklabels([x for x in eventlist if x not in ['PRRSV - nsp2TF','west nile virus']], rotation=45, horizontalalignment='right')
plt.ylim(-5,205)
plt.xlabel('')
plt.axhline(y=100, color='grey', linewidth=2)
plt.ylabel('% wild-type GFP fluorescence')
f.savefig('./figures/optimalityslippery.png', \
   dpi = 300, format='png', bbox_inches='tight', frameon=True)


df1[[x in eventlist for x in df1.fsevent]&
        (df1.endoslippery==False)&(df1.gfpm_wt<25)&(df1.percgfpm_wt<200)].pivot_table(
    index='fsevent',values='percgfpm_wt', aggfunc=np.count_nonzero).to_csv('./n_numbers/Fig2e.csv')


# for subgroup 3
df3=df[df.subgroup=='3']
df3['slippery']=df3.varseq162.apply(lambda x: x.upper()[35:42])
df3['slipperyfirst']=df3.slippery.apply(lambda x: x[:2])
df3['slipperysecond']=df3.slippery.apply(lambda x: x[2:])


for event in df3.fsevent.unique():
    f=plt.figure(figsize=(4,3))
    sns.heatmap(data=df3[(df3.fsevent==event)&(df3.numberreadsm_wt>=20)&(df3.gfpm_wt<25)].pivot_table(index='slipperyfirst',
                         columns='slipperysecond',values='gfpm_wt',aggfunc=np.median), 
        annot=df3[(df3.fsevent==event)&(df3.numberreadsm_wt>=20)&(df3.gfpm_wt<25)].pivot_table(index='slipperyfirst',
                  columns='slipperysecond',values='gfpm_wt',aggfunc=np.median), fmt='.1f')
    f.savefig('./figures/slippery/subgroup3_first_vs_second_median_'+'_'.join(event.split(' '))+'_gfp.png', \
       dpi = 300, format='png', bbox_inches='tight', frameon=True)
    
for event in df3.fsevent.unique():
    f=plt.figure(figsize=(4,3))
    sns.heatmap(data=df3[(df3.fsevent==event)&(df3.numberreadsm_wt>=20)&(df3.gfpm_wt<25)].pivot_table(index='slipperyfirst',
                         columns='slipperysecond',values='percgfpm_wt',aggfunc=np.median), 
        annot=df3[(df3.fsevent==event)&(df3.numberreadsm_wt>=20)&(df3.gfpm_wt<25)].pivot_table(index='slipperyfirst',
                  columns='slipperysecond',values='percgfpm_wt',aggfunc=np.median), fmt='.1f')
    f.savefig('./figures/slippery/subgroup3_first_vs_second_median_percgfp_'+'_'.join(event.split(' '))+'_gfp.png', \
       dpi = 300, format='png', bbox_inches='tight', frameon=True)



# for subgroup 4
df4=df[(df.subgroup=='4')&(df.laststopinframem1<12)]
df4['slippery']=df4.varseq162.apply(lambda x: x.upper()[35:42])
df4['slipperyfirst']=df4.slippery.apply(lambda x: x[:1])
df4['slipperysecond']=df4.slippery.apply(lambda x: x[1:3])
df4['slipperylast']=df4.slippery.apply(lambda x: x[3:])


for event in df4.fsevent.unique():
    f=plt.figure(figsize=(4,3))
    sns.heatmap(data=df4[(df4.fsevent==event)&(df4.numberreadsm_wt>=20)&(df4.gfpm_wt<25)].pivot_table(index='slipperyfirst',
                         columns='slipperysecond',values='gfpm_wt',aggfunc=np.median), 
        annot=df4[(df4.fsevent==event)&(df4.numberreadsm_wt>=20)&(df4.gfpm_wt<25)].pivot_table(index='slipperyfirst',
                  columns='slipperysecond',values='gfpm_wt',aggfunc=np.median), fmt='.1f')
    f.savefig('./figures/slippery/subgroup4_first_vs_second_median_'+'_'.join(event.split(' '))+'_gfp.png', \
       dpi = 300, format='png', bbox_inches='tight', frameon=True)

    f=plt.figure(figsize=(4,3))
    sns.heatmap(data=df4[(df4.fsevent==event)&(df4.numberreadsm_wt>=20)&(df4.gfpm_wt<25)].pivot_table(index='slipperysecond',
                         columns='slipperylast',values='gfpm_wt',aggfunc=np.median), 
        annot=df4[(df4.fsevent==event)&(df4.numberreadsm_wt>=20)&(df4.gfpm_wt<25)].pivot_table(index='slipperysecond',
                  columns='slipperylast',values='gfpm_wt',aggfunc=np.median), fmt='.1f')
    f.savefig('./figures/slippery/subgroup4_second_vs_last_median_'+'_'.join(event.split(' '))+'_gfp.png', \
       dpi = 300, format='png', bbox_inches='tight', frameon=True)

for event in df4.fsevent.unique():
    f=plt.figure(figsize=(4,3))
    sns.heatmap(data=df4[(df4.fsevent==event)&(df4.numberreadsm_wt>=20)&(df4.gfpm_wt<25)].pivot_table(index='slipperyfirst',
                         columns='slipperysecond',values='percgfpm_wt',aggfunc=np.median), 
        annot=df4[(df4.fsevent==event)&(df4.numberreadsm_wt>=20)&(df4.gfpm_wt<25)].pivot_table(index='slipperyfirst',
                  columns='slipperysecond',values='percgfpm_wt',aggfunc=np.median), fmt='.1f')
    f.savefig('./figures/slippery/subgroup4_first_vs_second_median_percgfp_'+'_'.join(event.split(' '))+'_gfp.png', \
       dpi = 300, format='png', bbox_inches='tight', frameon=True)


df['wtslippery']=df.index.map(lambda x: bool(df.slippery[x]==wtseqs[df.fsevent[x]][35:42]))


f=plt.figure(figsize=(4,3))
ax=sns.boxplot(data=df[(df.subgroup=='4')&(df.gfpm_wt<10)&(df.laststopinframem1<12)], x='fsevent', y='gfpm_wt', order=df[(df.subgroup=='4')&(df.gfpm_wt<10)].pivot_table(index='fsevent', values='gfpm_wt', aggfunc=np.median).sort_values(by='gfpm_wt').index, color=sns.xkcd_rgb['light blue'])
ax=sns.swarmplot(data=df[(df.subgroup=='4')&(df.gfpm_wt<10)&(df.laststopinframem1<12)], x='fsevent', y='gfpm_wt', order=df[(df.subgroup=='4')&(df.gfpm_wt<10)].pivot_table(index='fsevent', values='gfpm_wt', aggfunc=np.median).sort_values(by='gfpm_wt').index, color=sns.xkcd_rgb['medium blue'])
ax=sns.swarmplot(data=wtvalues.reset_index(), x='index', y='gfpm', order=df[(df.subgroup=='4')&(df.gfpm_wt<10)].pivot_table(index='fsevent', values='gfpm_wt', aggfunc=np.median).sort_values(by='gfpm_wt').index, color=sns.xkcd_rgb['red'])
ax.set_xticklabels(df[(df.subgroup=='4')&(df.gfpm_wt<10)&(df.laststopinframem1<12)].pivot_table(index='fsevent', values='gfpm_wt', aggfunc=np.median).sort_values(by='gfpm_wt').index, rotation=45, horizontalalignment='right')
plt.ylim(0,10.2)
plt.axhspan(0,1.3,color=sns.xkcd_rgb['medium blue'], alpha=0.2)
f.savefig('./figures/slippery/subgroup4_byevent_gfpm.png', \
   dpi = 300, format='png', bbox_inches='tight', frameon=True)


df[(df.subgroup=='4')&(df.gfpm_wt<10)&(df.laststopinframem1<12)].pivot_table(
    index='fsevent',values='percgfpm_wt', aggfunc=np.count_nonzero).to_csv('./n_numbers/Fig2g.csv')

# for individual fsevents

for event in df1.fsevent.unique():
    f=plt.figure(figsize=(4,3))
    sns.heatmap(data=df1[(df1.numberreadsm_wt>=20)&(df1.wavm_wt<12)&(df1.fsevent==event)].pivot_table(index='slipperyfirst',
                         columns='slipperysecond',values='percgfpm_wt',aggfunc=np.median), 
        annot=df1[(df1.numberreadsm_wt>=20)&(df1.wavm_wt<12)&(df1.fsevent==event)].pivot_table(index='slipperyfirst',
                  columns='slipperysecond',values='percgfpm_wt',aggfunc=np.median), fmt='.1f', center=100, vmin=0)
    plt.title(event + ' - ' + wtseqs[event][35:42])
    f.savefig('./figures/slippery/subgroup1_first_vs_second_median_'+'_'.join(event.split(' '))+'_gfpperc.png', \
       dpi = 300, format='png', bbox_inches='tight', frameon=True)
    
    
    f=plt.figure(figsize=(4,3))
    sns.heatmap(data=df1[(df1.numberreadsm_wt>=20)&(df1.wavm_wt<12)&(df1.fsevent==event)].pivot_table(index='slipperysecond',
                         columns='slipperylast',values='percgfpm_wt',aggfunc=np.median), 
        annot=df1[(df1.numberreadsm_wt>=20)&(df1.wavm_wt<12)&(df1.fsevent==event)].pivot_table(index='slipperysecond',
                  columns='slipperylast',values='percgfpm_wt',aggfunc=np.median), fmt='.1f', center=100, vmin=0)
    plt.title(event + ' - ' + wtseqs[event][35:42])
    f.savefig('./figures/slippery/subgroup1_second_vs_last_median_'+'_'.join(event.split(' '))+'_gfpperc.png', \
       dpi = 300, format='png', bbox_inches='tight', frameon=True)


for event in df1.fsevent.unique():
    f=plt.figure(figsize=(4,3))
    sns.heatmap(data=df1[(df1.numberreadsm_wt>=20)&(df1.wavm_wt<12)&(df1.fsevent==event)].pivot_table(index='slipperyfirst',
                         columns='slipperysecond',values='gfpm_wt',aggfunc=np.median), 
        annot=df1[(df1.numberreadsm_wt>=20)&(df1.wavm_wt<12)&(df1.fsevent==event)].pivot_table(index='slipperyfirst',
                  columns='slipperysecond',values='gfpm_wt',aggfunc=np.median), fmt='.1f', vmin=0)
    plt.title(event + ' - ' + wtseqs[event][35:42])
    f.savefig('./figures/slippery/subgroup1_first_vs_second_median_'+'_'.join(event.split(' '))+'_gfp.png', \
       dpi = 300, format='png', bbox_inches='tight', frameon=True)
    
    
    f=plt.figure(figsize=(4,3))
    sns.heatmap(data=df1[(df1.numberreadsm_wt>=20)&(df1.wavm_wt<12)&(df1.fsevent==event)].pivot_table(index='slipperysecond',
                         columns='slipperylast',values='gfpm_wt',aggfunc=np.median), 
        annot=df1[(df1.numberreadsm_wt>=20)&(df1.wavm_wt<12)&(df1.fsevent==event)].pivot_table(index='slipperysecond',
                  columns='slipperylast',values='gfpm_wt',aggfunc=np.median), fmt='.1f', vmin=0)
    plt.title(event + ' - ' + wtseqs[event][35:42])
    f.savefig('./figures/slippery/subgroup1_second_vs_last_median_'+'_'.join(event.split(' '))+'_gfp.png', \
       dpi = 300, format='png', bbox_inches='tight', frameon=True)


for event in df1.fsevent.unique():
    f=plt.figure(figsize=(4,3))
    sns.heatmap(data=df1[(df1.numberreadsm_wt>=20)&(df1.wavm_wt<12)&(df1.fsevent==event)].pivot_table(index='slipperyfirst',
                         columns='slipperysecond',values='gfpm_wt',aggfunc=np.median), 
        annot=df1[(df1.numberreadsm_wt>=20)&(df1.wavm_wt<12)&(df1.fsevent==event)].pivot_table(index='slipperyfirst',
                  columns='slipperysecond',values='gfpm_wt',aggfunc=np.count_nonzero), fmt='.1f', vmin=0)
    plt.title(event + ' - ' + wtseqs[event][35:42])
    f.savefig('./figures/slippery/subgroup1_first_vs_second_median_'+'_'.join(event.split(' '))+'_gfp_withnnumbers.png', \
       dpi = 300, format='png', bbox_inches='tight', frameon=True)
    
    
    f=plt.figure(figsize=(4,3))
    sns.heatmap(data=df1[(df1.numberreadsm_wt>=20)&(df1.wavm_wt<12)&(df1.fsevent==event)].pivot_table(index='slipperysecond',
                         columns='slipperylast',values='gfpm_wt',aggfunc=np.median), 
        annot=df1[(df1.numberreadsm_wt>=20)&(df1.wavm_wt<12)&(df1.fsevent==event)].pivot_table(index='slipperysecond',
                  columns='slipperylast',values='gfpm_wt',aggfunc=np.count_nonzero), fmt='.1f', vmin=0)
    plt.title(event + ' - ' + wtseqs[event][35:42])
    f.savefig('./figures/slippery/subgroup1_second_vs_last_median_'+'_'.join(event.split(' '))+'_gfp_withnnumbers.png', \
       dpi = 300, format='png', bbox_inches='tight', frameon=True)


################
# For OAZ frameshifting site - systematic replacement of codons

df2=df[(df.subgroup=='2')&(df.numberreadsp_wt>=20)&(df.peaksp_wt==1)]
df2['slippery']=df2.varseq162.apply(lambda x: x.upper()[36:42])
df2['slipperyfirst']=df2.slippery.apply(lambda x: x[:3])
df2['slipperysecond']=df2.slippery.apply(lambda x: x[3:6])
df2['slipperyfirstaa']=df2.slipperyfirst.apply(lambda x: str(Seq(x).translate()))
df2['slipperysecondaa']=df2.slipperysecond.apply(lambda x: str(Seq(x).translate()))

df2['firstorsecond']=df2.slippery.apply(lambda x: 'first' if x[3:6]=='TGA' else 'second')


# plot TGA vs TAA and TAG at the second position - with all the others
df2[(df2.slipperysecondaa=='*')&(df2.slipperyfirst=='TCC')].index

f=plt.figure(figsize=(4,2))
plt.plot(xval, covp1comb.loc['38345.0_wt'], 'g')
plt.plot(xval, covp1comb.loc['38346.0_wt'], 'b')
plt.plot(xval, covp1comb.loc['38349.0_wt'], 'magenta')
for var in list(df2[(df2.wavp_wt<12)&(df2.slipperysecondaa!='*')&(df2.slipperyfirst=='TCC')].index.values):
    plt.plot(xval, covp1comb.loc[str(var)+'.0_wt'], 'gray', alpha=0.5)
f.savefig('./figures/slippery/oaz_slipperysecond_variants_cov.png', \
   dpi = 300, format='png', bbox_inches='tight', frameon=True)


# plot options at the second position - with all the others
f=plt.figure(figsize=(4,2))
plt.plot(xval, covp1comb.loc['38346.0_wt'], 'b')
var=df2[(df2.slipperysecond=='TGA')&(df2.slipperyfirst=='TTT')].index.values[0]
plt.plot(xval, covp1comb.loc[str(var)+'.0_wt'], 'g')
for var in list(df2[(df2.wavp_wt<12)&(df2.slipperysecond=='TGA')&(df2.slipperyfirst!='TCC')&(df2.slipperyfirst!='TTT')].index.values):
    plt.plot(xval, covp1comb.loc[str(var)+'.0_wt'], 'gray', alpha=0.5)
f.savefig('./figures/slippery/oaz_slipperyfirst_variants_cov.png', \
   dpi = 300, format='png', bbox_inches='tight', frameon=True)


#%%
df6=df[(df.subgroup=='6')&(df.numberreadsm_wt>=20)&(df.peaksm_wt==1)]
df6['slippery']=df6.varseq162.apply(lambda x: x.upper()[35:42].replace('T','U'))
df6['slipperybefore']=df6.varseq162.apply(lambda x: x.upper()[33:36])
df6['slipperyfirst']=df6.slippery.apply(lambda x: x[1:4])
df6['slipperysecond']=df6.slippery.apply(lambda x: x[4:7])
df6['slipperyshifted']=df6.slippery.apply(lambda x: x[3:6])
df6['slipperybeforeaa']=df6.slipperybefore.apply(lambda x: str(Seq(x).translate()))
df6['slipperyfirstaa']=df6.slipperyfirst.apply(lambda x: str(Seq(x).translate()))
df6['slipperysecondaa']=df6.slipperysecond.apply(lambda x: str(Seq(x).translate()))
df6['slipperyshiftedaa']=df6.slipperyshifted.apply(lambda x: str(Seq(x).translate()))

df6['secondvsshifted']=df6.index.map(lambda x: bool(df6.slipperysecond[x]==df6.slipperyshifted[x]))
df6['secondvsshiftedaa']=df6.index.map(lambda x: bool(df6.slipperysecondaa[x]==df6.slipperyshiftedaa[x]))
df6['secondvsshifted3rdpos']=df6.index.map(lambda x: bool(df6.slipperysecond[x][0:2]==df6.slipperyshifted[x][0:2]))


def findmutpos(x):
    mutpos=8
    endo=fsdf[(fsdf.fsevent==df.fsevent[x])&(fsdf.subset=='fsmain')].varseq162.unique()[0]
    for i in range(35,42):
        if df.varseq162[x][i]!=endo[i]:
            mutpos=i-34
            break
    if mutpos<8:
        return mutpos
    else:
        return 0
            

df6['mutpos']=df6.index.map(lambda x: findmutpos(x))

df1=df[df.subgroup=='1']
df1['slippery']=df1.varseq162.apply(lambda x: x.upper()[35:42].replace('T','U'))
df1['slipperylast']=df1.slippery.apply(lambda x: x[-1:])

df1['mutpos']=df1.index.map(lambda x: findmutpos(x))

df6=df6.append(df1[df1.mutpos==7])
df6['mutnt']=df6.index.map(lambda x: df6.slippery[x][df6.mutpos[x]-1])

f=plt.figure(figsize=(4,3))
sns.boxplot(data=df6[([x in eventlist for x in df6.fsevent])&(df6.percgfpm_wt<200)], x='mutpos', y='percgfpm_wt', color=sns.xkcd_rgb['light blue'], linewidth=2)
sns.swarmplot(data=df6[([x in eventlist for x in df6.fsevent])&(df6.percgfpm_wt<200)], x='mutpos', y='percgfpm_wt', color=sns.xkcd_rgb['dark blue'])
plt.axhline(y=100, linewidth=2, color='grey')
plt.ylim(0,200)
f.savefig('./figures/slippery/slipperyvar_subgroup16_eventlist_gfpperc_per_slipperyposition.png', \
          dpi = 300, format='png', bbox_inches='tight', frameon=True)

df6[([x in eventlist for x in df6.fsevent])&(df6.percgfpm_wt<200)].pivot_table(
    index='mutpos',values='percgfpm_wt', aggfunc=np.count_nonzero).to_csv('./n_numbers/Fig2c.csv')


f=plt.figure(figsize=(4,3))
sns.boxplot(data=df6[([x in nofslist for x in df6.fsevent])&(df6.percgfpm_wt<200)], x='mutpos', y='percgfpm_wt', color=sns.xkcd_rgb['light blue'], linewidth=2)
sns.swarmplot(data=df6[([x in nofslist for x in df6.fsevent])&(df6.percgfpm_wt<200)], x='mutpos', y='percgfpm_wt', color=sns.xkcd_rgb['dark blue'])
plt.axhline(y=100, linewidth=2, color='grey')
plt.ylim(0,200)
f.savefig('./figures/slippery/slipperyvar_subgroup16_nofslist_gfpperc_per_slipperyposition.png', \
          dpi = 300, format='png', bbox_inches='tight', frameon=True)

   


for i in eventlist:
    try:
        f=plt.figure(figsize=(len(df6[(df6.fsevent==i)&(df6.gfpm_wt<25)&(df6.numberreadsm_wt>=20)].slippery.unique())/2,3))
        ax=sns.swarmplot(data=df6[(df6.fsevent==i)&(df6.gfpm_wt<25)&(df6.numberreadsm_wt>=20)], x='slippery',y='gfpm_wt', color=sns.xkcd_rgb['medium blue'])
        ax.set_xticklabels(df6[(df6.fsevent==i)&(df6.gfpm_wt<25)&(df6.numberreadsm_wt>=20)].slippery.unique(),rotation=45, horizontalalignment='right')
        plt.title(i, fontsize=14)
        plt.axhline(y=wtvalues.gfpm[i], color='grey', linewidth=2)
        plt.axhspan(0,1.3,color=sns.xkcd_rgb['medium blue'], alpha=0.2)
        plt.ylim(ymin=0)
        f.savefig('./figures/slippery/slipperyvar_subgroup16_slipperymutantsperevent_'+'_'.join(i.split(' '))+'_U.png', \
          dpi = 300, format='png', bbox_inches='tight', frameon=True)

    except:
        pass

for i in nofslist:
    try:
        f=plt.figure(figsize=(len(df6[(df6.fsevent==i)&(df6.gfpm_wt<25)&(df6.numberreadsm_wt>=20)].slippery.unique())/2,3))
        ax=sns.swarmplot(data=df6[(df6.fsevent==i)&(df6.gfpm_wt<25)&(df6.numberreadsm_wt>=20)], x='slippery',y='gfpm_wt', color=sns.xkcd_rgb['medium blue'])
        ax.set_xticklabels(df6[(df6.fsevent==i)&(df6.gfpm_wt<25)&(df6.numberreadsm_wt>=20)].slippery.unique(),rotation=45, horizontalalignment='right')
        plt.title(i, fontsize=14)
        plt.axhline(y=wtvalues.gfpm[i], color='grey', linewidth=2)
        plt.axhspan(0,1.3,color=sns.xkcd_rgb['medium blue'], alpha=0.2)
        plt.ylim(ymin=0)
        f.savefig('./figures/slippery/slipperyvar_subgroup16_slipperymutantsperevent_'+'_'.join(i.split(' '))+'_U.png', \
          dpi = 300, format='png', bbox_inches='tight', frameon=True)

    except:
        pass



############# for figure 1

nofslistnoHSV=[x for x in nofslist if x!='herpes simplex - thymidine kinase']

# for minus 1

indexes=[]
samevarseq=[]
for var, group in fsdf[fsdf.subset=='fsmain'].groupby(by='varseq162'):
    samevarseq.append(var)
    indexes.append(list(fsdf[fsdf.varseq162==var].index))


df=fsdf[(fsdf.numberreadsm_wt>=20)&(fsdf.peaksm_wt==1)&(fsdf.laststopinframem1<12)]
bcpeakcorr = pd.DataFrame(columns=['gfpm_wt','fsevent'])
for i in samevarseq:
    indexes=list(df[(df.varseq162==i)&(df.gfpm_wt<25)].index)
    fsindexes=fsdf.loc[indexes,['gfpm_wt','fsevent']]
    event=fsdf[(fsdf.subset=='fsmain')&(fsdf.varseq162==i)].fsevent.values[0]
    fsindexes['fsevent']=event
    bcpeakcorr=pd.concat([bcpeakcorr,fsindexes])

bcpeakcorr['controls']='wild-type'


df=fsdf[(fsdf.numberreadsm_del>=20)&(fsdf.peaksm_del==1)]
bcpeakcorr2 = pd.DataFrame(columns=['gfpm_del','fsevent'])
for i in samevarseq:
    indexes=list(df[(df.varseq162==i)&(df.gfpm_del<15)].index)
    fsindexes=fsdf.loc[indexes,['gfpm_del','fsevent']]
    event=fsdf[(fsdf.subset=='fsmain')&(fsdf.varseq162==i)].fsevent.values[0]
    fsindexes['fsevent']=event
    bcpeakcorr2=pd.concat([bcpeakcorr2,fsindexes])
           
bcpeakcorr2['controls']='single nucleotide deletion'

bcpeakcorr.rename(columns={'gfpm_wt':'gfpm'}, inplace=True)
bcpeakcorr2.rename(columns={'gfpm_del':'gfpm'}, inplace=True)



df=fsdf[(fsdf.numberreadsstop_wt>=20)&(fsdf.peaksstop_wt==1)]
bcpeakcorr4 = pd.DataFrame(columns=['gfpstop_wt','fsevent'])
for i in samevarseq:
    indexes=list(df[(df.varseq162==i)&(df.gfpstop_wt<25)].index)
    fsindexes=fsdf.loc[indexes,['gfpstop_wt','fsevent']]
    event=fsdf[(fsdf.subset=='fsmain')&(fsdf.varseq162==i)].fsevent.values[0]
    fsindexes['fsevent']=event
    bcpeakcorr4=pd.concat([bcpeakcorr4,fsindexes])


bcpeakcorr4['controls']='stop control'
bcpeakcorr4.rename(columns={'gfpstop_wt':'gfpm'}, inplace=True)


          
controlscomb=pd.concat([bcpeakcorr, bcpeakcorr2, bcpeakcorr4], ignore_index=True)
#nofslist.remove('herpes simplex - thymidine kinase')

f=plt.figure(figsize=(6,3))
ax=sns.pointplot(data=controlscomb[[x in eventlist+nofslistnoHSV for x in controlscomb.fsevent]], x='fsevent',y='gfpm', hue='controls',order=bcpeakcorr[[x in eventlist+nofslistnoHSV for x in bcpeakcorr.fsevent]].pivot_table(index='fsevent',values='gfpm').sort_values(by='gfpm').index, \
                 scale=0.6, join=False, dodge=0.5, errwidth=3)
plt.ylim(0,12)
ax.set_xticklabels(bcpeakcorr[[x in eventlist+nofslistnoHSV for x in bcpeakcorr.fsevent]].pivot_table(index='fsevent',values='gfpm').sort_values(by='gfpm').index, rotation=45, horizontalalignment='right')
ax.legend_.remove()
plt.xlabel('')
plt.ylabel('% GFP fluorescence')
plt.axhspan(0,1.3, color=sns.xkcd_rgb['medium blue'], alpha=0.2)
#f.savefig('./figures/fsmain_multiple_barcodes_GFP_snd_stop.png', \
#          dpi = 300, format='png', bbox_inches='tight', frameon=True)
f.savefig('./figures/fsmain_multiple_barcodes_GFP_snd_stop_withoutHSV_nostop.png', \
          dpi = 300, format='png', bbox_inches='tight', frameon=True)
f.show()


controlscomb[[x in eventlist+nofslistnoHSV for x in controlscomb.fsevent]].pivot_table(index='fsevent', columns='controls', values='gfpm', aggfunc=np.count_nonzero).to_csv('./n_numbers/Fig1c.csv')

for event in controlscomb.fsevent.unique():
    print(event)
    print(scipy.stats.mannwhitneyu(controlscomb[(controlscomb.fsevent==event)&(controlscomb.controls=='wild-type')].gfpm,\
        controlscomb[(controlscomb.fsevent==event)&(controlscomb.controls!='wild-type')].gfpm))
'''
human - PEG10
MannwhitneyuResult(statistic=99.0, pvalue=0.089706932778028547)
simian srv1 retrovirus
MannwhitneyuResult(statistic=13.0, pvalue=8.4693094590744522e-05)
PLRV luteovirus
MannwhitneyuResult(statistic=164.0, pvalue=9.2250003493727327e-07)
Rous sarcoma virus
MannwhitneyuResult(statistic=349.0, pvalue=0.48244774831262111)
HIV HXB2
MannwhitneyuResult(statistic=21.0, pvalue=2.014283110299411e-08)
human T-lymphotropic virus
MannwhitneyuResult(statistic=27.0, pvalue=0.0030086435028803186)
SIVmac239
MannwhitneyuResult(statistic=37.0, pvalue=8.1932263743166577e-07)
human - HERV-K10
MannwhitneyuResult(statistic=1.0, pvalue=2.5827962097199238e-09)
SARS coronavirus
MannwhitneyuResult(statistic=8.0, pvalue=4.7530998483694518e-06)
human - ccr5
MannwhitneyuResult(statistic=5.0, pvalue=0.33027460260083646)
west nile virus
MannwhitneyuResult(statistic=44.0, pvalue=0.026964735186385206)

'''
# Test how many variants passed filtering for mCherry

df=fsdf[fsdf.library=='fs_designed']

df['passedfiltering']=df.gfpm_wt.apply(lambda x: 1 if x>0 else 0)


fractionpassed=pd.Series()
for var, group in fsdf[fsdf.subset=='fsmain'].groupby(by='fsevent'):
    if var!='herpes simplex - thymidine kinase':
        fractionpassed.loc[var]=float(df.loc[list(fsdf[fsdf.varseq162==fsdf[(fsdf.subset=='fsmain')&(fsdf.fsevent==var)].varseq162.values[0]].index),'passedfiltering'].sum())/\
                          len(list(fsdf[fsdf.varseq162==fsdf[(fsdf.subset=='fsmain')&(fsdf.fsevent==var)].varseq162.values[0]].index))
    

f=plt.figure(figsize=(6,3))
ax=fractionpassed.sort_values(ascending=False).plot(kind='bar', color=sns.xkcd_rgb['medium blue'])    
ax.set_xticklabels(ax.get_xticklabels(), rotation=45, horizontalalignment='right')
plt.ylabel('fraction of wild-type sequences\npassing filtering') 
f.savefig('./figures/fraction_of_sequences_passing_filtering.png', \
          dpi = 300, format='png', bbox_inches='tight', frameon=True)


npassed=pd.DataFrame()
for var, group in fsdf[fsdf.subset=='fsmain'].groupby(by='fsevent'):
    if var!='herpes simplex - thymidine kinase':
        npassed.loc[var,'passed'],npassed.loc[var,'total']=float(df.loc[list(fsdf[fsdf.varseq162==fsdf[(fsdf.subset=='fsmain')&(fsdf.fsevent==var)].varseq162.values[0]].index),'passedfiltering'].sum()),\
                          len(list(fsdf[fsdf.varseq162==fsdf[(fsdf.subset=='fsmain')&(fsdf.fsevent==var)].varseq162.values[0]].index))


npassed.to_csv('./n_numbers/FigS4d.csv')

# for plus 1

df=fsdf[(fsdf.numberreadsp_wt>=20)&(fsdf.peaksp_wt==1)&(fsdf.laststopinframep1<12)]
bcpeakcorr = pd.DataFrame(columns=['gfpp_wt','fsevent'])
for i in samevarseq:
    indexes=list(df[(df.varseq162==i)&(df.gfpp_wt<25)].index)
    fsindexes=fsdf.loc[indexes,['gfpp_wt','fsevent']]
    event=fsdf[(fsdf.subset=='fsmain')&(fsdf.varseq162==i)].fsevent.values[0]
    fsindexes['fsevent']=event
    bcpeakcorr=pd.concat([bcpeakcorr,fsindexes])


bcpeakcorr['controls']='wild-type'


df=fsdf[(fsdf.numberreadsp_del>=20)&(fsdf.peaksp_del==1)]
bcpeakcorr2 = pd.DataFrame(columns=['gfpm_del','fsevent'])
for i in samevarseq:
    indexes=list(df[(df.varseq162==i)&(df.gfpp_del<25)].index)
    fsindexes=fsdf.loc[indexes,['gfpp_del','fsevent']]
    event=fsdf[(fsdf.subset=='fsmain')&(fsdf.varseq162==i)].fsevent.values[0]
    fsindexes['fsevent']=event
    bcpeakcorr2=pd.concat([bcpeakcorr2,fsindexes])
           
bcpeakcorr2['controls']='single nucleotide deletion'

bcpeakcorr.rename(columns={'gfpp_wt':'gfpp'}, inplace=True)
bcpeakcorr2.rename(columns={'gfpp_del':'gfpp'}, inplace=True)



df=fsdf[(fsdf.numberreadsstop_wt>=20)&(fsdf.peaksstop_wt==1)]
bcpeakcorr4 = pd.DataFrame(columns=['gfpstop_wt','fsevent'])
for i in samevarseq:
    indexes=list(df[(df.varseq162==i)&(df.gfpstop_wt<25)].index)
    fsindexes=fsdf.loc[indexes,['gfpstop_wt','fsevent']]
    event=fsdf[(fsdf.subset=='fsmain')&(fsdf.varseq162==i)].fsevent.values[0]
    fsindexes['fsevent']=event
    bcpeakcorr4=pd.concat([bcpeakcorr4,fsindexes])


bcpeakcorr4['controls']='stop control'
bcpeakcorr4.rename(columns={'gfpstop_wt':'gfpp'}, inplace=True)


          
controlscomb=pd.concat([bcpeakcorr, bcpeakcorr2, bcpeakcorr4], ignore_index=True)


f=plt.figure(figsize=(3,3))
ax=sns.pointplot(data=controlscomb[[x in eventlist+nofslistnoHSV+plist for x in controlscomb.fsevent]], 
                 x='fsevent',y='gfpp', hue='controls',
                 order=bcpeakcorr[[x in eventlist+nofslistnoHSV+plist for x in bcpeakcorr.fsevent]].pivot_table(index='fsevent',values='gfpp').sort_values(by='gfpp').index, \
                 scale=0.6, join=False, dodge=0.5, errwidth=3)
plt.ylim(0,20)
ax.set_xticklabels(ax.get_xticklabels(), rotation=45, horizontalalignment='right')
ax.legend_.remove()
plt.xlabel('')
plt.ylabel('% GFP fluorescence')
plt.axhspan(0,1.3, color=sns.xkcd_rgb['medium blue'], alpha=0.2)
#f.savefig('./figures/fsmain_multiple_barcodes_GFPp1_snd_stop_allevents.png', \
#          dpi = 300, format='png', bbox_inches='tight', frameon=True)
f.savefig('./figures/fsmain_multiple_barcodes_GFPp1_snd_stop_allevents_withoutHSV.png', \
          dpi = 300, format='png', bbox_inches='tight', frameon=True)
f.show()

controlscomb[[x in eventlist+nofslistnoHSV for x in controlscomb.fsevent]].pivot_table(index='fsevent', columns='controls', values='gfpm', aggfunc=np.count_nonzero).to_csv('./n_numbers/Fig1d.csv')

for event in controlscomb.fsevent.unique():
    print(event)
    print(scipy.stats.ttest_ind(controlscomb[(controlscomb.fsevent==event)&(controlscomb.controls=='wild-type')].gfpp,\
        controlscomb[(controlscomb.fsevent==event)&(controlscomb.controls!='wild-type')].gfpp))
'''
influenza a virus
Ttest_indResult(statistic=4.786143711691305, pvalue=0.0019980252621817677)
PRRSV - nsp2TF
Ttest_indResult(statistic=0.83670141022975963, pvalue=0.42056588455850441)
HCV - F protein
Ttest_indResult(statistic=nan, pvalue=nan)
human - OAZ1
Ttest_indResult(statistic=14.546633411806818, pvalue=1.4023442383259715e-14)

'''



wtvaluesall=pd.DataFrame()
for var in fsdf[(fsdf.subset=='fsmain')|(fsdf.changes=='endogenous sequence')].fsevent.unique():
    vrsq=fsdf[((fsdf.subset=='fsmain')|(fsdf.changes=='endogenous sequence'))&(fsdf.fsevent==var)].varseq162.values[0]
    indexes=list(fsdf[(fsdf.gfpm_wt<25)&(fsdf.varseq162==vrsq)&(fsdf.peaksm_wt==1)&(fsdf.numberreadsm_wt>=20)].index)  
    wtvaluesall.loc[var,'gfpm']=fsdf.loc[indexes,'gfpm_wt'].median()
    indexes=list(fsdf[(fsdf.gfpp_wt<25)&(fsdf.varseq162==vrsq)&(fsdf.peaksp_wt==1)&(fsdf.numberreadsp_wt>=20)].index)  
    wtvaluesall.loc[var,'gfpp']=fsdf.loc[indexes,'gfpp_wt'].median()
    indexes=list(fsdf[(fsdf.gfpstop_wt<25)&(fsdf.varseq162==vrsq)&(fsdf.peaksstop_wt==1)&(fsdf.numberreadsstop_wt>=20)].index)  
    wtvaluesall.loc[var,'gfpstop']=fsdf.loc[indexes,'gfpstop_wt'].median()

    
f=plt.figure(figsize=(3,3))
plt.scatter(wtvaluesall.gfpm, wtvaluesall.gfpp)
plt.xlim(-0.2,10.2)
plt.ylim(-0.2,20.2)
plt.axhspan(-0.2,1.3, color=sns.xkcd_rgb['medium blue'], alpha=0.2)
plt.axvspan(-0.2,1.3, color=sns.xkcd_rgb['medium blue'], alpha=0.2)
f.savefig('./figures/readout_gfpm_vs_gfpp.png', \
          dpi = 300, format='png', bbox_inches='tight', frameon=True)
f.show()



### For rnaperdna

indexes=[]
samevarseq=[]
for var, group in fsdf[fsdf.subset=='fsmain'].groupby(by='varseq162'):
    samevarseq.append(var)
    indexes.append(list(fsdf[fsdf.varseq162==var].index))


df=fsdf[(fsdf.numberreadsm_wt>=20)&(fsdf.peaksm_wt==1)&(fsdf.fsevent!='herpes simplex - thymidine kinase')]
bcpeakcorr = pd.DataFrame(columns=['rnaperdna','fsevent'])
for i in samevarseq:
    indexes=list(df[(df.varseq162==i)&(df.gfpm_wt<25)].index)
    fsindexes=fsdf.loc[indexes,['rnaperdna','fsevent']]
    event=fsdf[(fsdf.subset=='fsmain')&(fsdf.varseq162==i)].fsevent.values[0]
    fsindexes['fsevent']=event
    bcpeakcorr=pd.concat([bcpeakcorr,fsindexes])

eventlist.reverse()
f=plt.figure(figsize=(8,3))
ax=sns.boxplot(data=bcpeakcorr, x='fsevent', y='rnaperdna', fliersize=0, color=sns.xkcd_rgb['light blue'],
               order=eventlist+nofslistnoHSV+plist)
ax=sns.swarmplot(data=bcpeakcorr, x='fsevent', y='rnaperdna', color=sns.xkcd_rgb['medium blue'],
                 order=eventlist+nofslistnoHSV+plist)
ax.set_xticklabels(eventlist+nofslistnoHSV+plist, rotation=45, horizontalalignment='right')
plt.ylim(-9,5)
plt.xlabel('')
plt.ylabel('log2(RNA/DNA reads)')
f.savefig('./figures/rnaperdna_per_PRF_event_wildtype.png', \
          dpi = 300, format='png', bbox_inches='tight', frameon=True)
eventlist.reverse()

bcpeakcorr.pivot_table(index='fsevent',values='rnaperdna',aggfunc=np.count_nonzero).to_csv('./n_numbers/FigS4e.csv')

#%%

###################
# SPACER
##################

df=fsdf[((fsdf.subset=="fs_spacervar")|(fsdf.subset=='fsmain'))&(fsdf.numberreadsm_wt>=10)&(fsdf.peaksm_wt==1)&\
        (fsdf.laststopinframem1<12)]

df['subgroup']=df.changes.apply(lambda x: int(x[9]) if type(x)==str else 4)

def parsechanges(x):
    if type(x)==str:
        if 'deleted' in x:
            ntchanges=-int(x.split(': ')[1][0])
            template=''
            stopmutant='stopmutant' in x
        elif 'inserted' in x:
            ntchanges=int(x.split(': ')[1][0])
            template=x.split('sequence ')[1].split(' ')[0]
            stopmutant='stopmutant' in x
        else:
            ntchanges=np.nan
            template=''       
            stopmutant='stopmutant' in x
    else:
        ntchanges=np.nan
        template=''       
        stopmutant=False
    return pd.Series([ntchanges, template, stopmutant], index=['ntchanges','template','stopmutant'])
            
df[['ntchanges','template','stopmutant']]=df.changes.apply(lambda x: parsechanges(x))

df['ntchanges']=df.index.map(lambda x: 0 if df.subset[x]=='fsmain' else df.ntchanges[x])



f=plt.figure(figsize=(8,3))
sns.boxplot(data=df[([x in eventlist for x in df.fsevent])&(df.gfpm_wt<25)], x='ntchanges', y='percgfpm_wt',\
                      color=sns.xkcd_rgb['light blue'])
ax=sns.swarmplot(data=df[([x in eventlist for x in df.fsevent])&(df.gfpm_wt<25)], x='ntchanges', y='percgfpm_wt',\
                      color=sns.xkcd_rgb['medium blue'])
ax.set_xticklabels(range(-6,10,1))
plt.ylim(0,220)
plt.legend(bbox_to_anchor=(1.55,1))
plt.axhline(y=100, color='gray', linewidth=2)
f.savefig('./figures/spacer/spacer_insertionsanddeletions_eventlist_percgfpm.png', \
  dpi = 300, format='png', bbox_inches='tight', frameon=True)

df[([x in eventlist for x in df.fsevent])&(df.gfpm_wt<25)].pivot_table(
        index='ntchanges',values='percgfpm_wt',aggfunc=np.count_nonzero).to_csv('./n_numbers/Fig4b.csv')

f=plt.figure(figsize=(8,3))
sns.boxplot(data=df[([x in eventlist for x in df.fsevent])&(df.gfpm_wt<25)], x='ntchanges', y='gfpm_wt',\
                      color=sns.xkcd_rgb['light blue'])
ax=sns.swarmplot(data=df[([x in eventlist for x in df.fsevent])&(df.gfpm_wt<25)], x='ntchanges', y='gfpm_wt',\
                      color=sns.xkcd_rgb['medium blue'])
ax.set_xticklabels(range(-6,10,1))
plt.ylim(0,6)
plt.legend(bbox_to_anchor=(1.55,1))
plt.axhline(y=100, color='gray', linewidth=2)
f.savefig('./figures/spacer/spacer_insertionsanddeletions_eventlist_gfpm.png', \
  dpi = 300, format='png', bbox_inches='tight', frameon=True)



#############

df3=df[df.subgroup==3]

wtseqs=pd.Series()

for event in fsdf[fsdf.subset=='fsmain'].fsevent.unique():
    wtseqs.loc[event]=fsdf[(fsdf.subset=='fsmain')&(fsdf.fsevent==event)].varseq162.values[0]

df3['codon']=df3.changes.apply(lambda x: x.split('Codon ')[-1][0:3])
df3['aa']=df3.codon.apply(lambda x: str(Seq(x).translate()))

df3['wtcodon']=df3.fsevent.apply(lambda x: wtseqs.loc[x][41:44])
df3['wtaa']=df3.wtcodon.apply(lambda x: str(Seq(x).translate()))

codonorder=np.sort(df3.codon.unique())
aaorder=np.sort(df3.aa.unique())


###

f=plt.figure(figsize=(8,2))
ax=sns.heatmap(data=df3[([x in eventlist for x in df3.fsevent])&
                        (df3.gfpm_wt<25)].pivot_table(index='wtcodon',columns='codon', values='gfpm_wt'))
f.savefig('./figures/spacer/spacervar_subgroup3_heatmap_eventlist_gfp_codon_by_wtcodon.png', \
  dpi = 300, format='png', bbox_inches='tight', frameon=True)


f=plt.figure(figsize=(8,2))
ax=sns.heatmap(data=df3[([x in eventlist for x in df3.fsevent])&
                        (df3.gfpm_wt<25)].pivot_table(index='wtcodon',columns='codon', values='percgfpm_wt'), center=100)
f.savefig('./figures/spacer/spacervar_subgroup3_heatmap_eventlist_gfpperc_codon_by_wtcodon.png', \
  dpi = 300, format='png', bbox_inches='tight', frameon=True)

f=plt.figure(figsize=(6,1.5))
ax=sns.heatmap(data=df3[((df3.smoothedpeaks==1)|(df3.rawpeaks==1))&([x in eventlist for x in df3.fsevent])&
                        (df3.gfpm_wt<25)].pivot_table(index='wtaa',columns='aa', values='percgfpm_wt'), center=100)
f.savefig('./figures/spacer/spacervar_subgroup3_heatmap_eventlist_gfpperc_aa_by_wtaa.png', \
  dpi = 300, format='png', bbox_inches='tight', frameon=True)


### tRNA adaptation index vs frameshifting efficiency


tai0=pd.DataFrame(index=df3.index)

def gettaimap(vrsq):
    tai=[]
    for i in range(len(vrsq)/3-1):
        tai.append(tAIdict[vrsq[0+3*i:3+3*i]])
    return pd.Series(tai)


def gettaidiffmap(vrsq):
    tai=[]
    for i in range(len(vrsq)/3-1):
        tai.append(tAIdict[vrsq[2+3*i:5+3*i]]-tAIdict[vrsq[3+3*i:6+3*i]])
    return pd.Series(tai)

tai0=df3.varseq162.apply(lambda x: gettaimap(x.upper()))
taim1=df3.varseq162.apply(lambda x: gettaimap(x[2:].upper()))

tai0all=fsdf[(fsdf.numberreadsm_wt>=20)&(fsdf.peaksm_wt==1)&(fsdf.gfpm_wt<25)&(fsdf.library=='fs_designed')&
             (['Y' not in vs for vs in fsdf.varseq162])].varseq162.apply(lambda x: gettaimap(x.upper()))
taim1all=fsdf[(fsdf.numberreadsm_wt>=20)&(fsdf.peaksm_wt==1)&(fsdf.gfpm_wt<25)&(fsdf.library=='fs_designed')&
              (['Y' not in vs for vs in fsdf.varseq162])].varseq162.apply(lambda x: gettaimap(x[2:].upper()))

tai0slip=fsdf[(fsdf.numberreadsm_wt>=20)&(fsdf.peaksm_wt==1)&(fsdf.gfpm_wt<25)&(fsdf.subset=='fs_slipperyvar')&
              (['Y' not in vs for vs in fsdf.varseq162])].varseq162.apply(lambda x: gettaimap(x.upper()))
taim1slip=fsdf[(fsdf.numberreadsm_wt>=20)&(fsdf.peaksm_wt==1)&(fsdf.gfpm_wt<25)&(fsdf.subset=='fs_slipperyvar')&
               (['Y' not in vs for vs in fsdf.varseq162])].varseq162.apply(lambda x: gettaimap(x[2:].upper()))



############### only for relevant regions

with sns.axes_style(style='dark'):
    f, ax=plt.subplots(figsize=(4,3))
    ax.plot(tai0slip.transpose().iloc[0:15].index.map(lambda x: pearsonr(tai0slip.transpose().loc[x], 
            fsdf[(fsdf.numberreadsm_wt>=20)&(fsdf.peaksm_wt==1)&(fsdf.gfpm_wt<25)&(fsdf.subset=='fs_slipperyvar')&(['Y' not in vs for vs in fsdf.varseq162])].gfpm_wt)[0]))
    ax.tick_params('y',colors=sns.xkcd_rgb['medium blue'])
    plt.axhline(y=0, linewidth=2, color=sns.xkcd_rgb['medium blue'])
    ax2=ax.twinx()
    ax2.plot(np.log10(tai0slip.transpose().iloc[0:15].index.map(lambda x: pearsonr(tai0slip.transpose().loc[x], 
            fsdf[(fsdf.numberreadsm_wt>=20)&(fsdf.peaksm_wt==1)&(fsdf.gfpm_wt<25)&(fsdf.subset=='fs_slipperyvar')&(['Y' not in vs for vs in fsdf.varseq162])].gfpm_wt)[1])), 
        color=sns.xkcd_rgb['green'], alpha=0.5)
    ax2.tick_params('y',colors=sns.xkcd_rgb['green'])
    plt.axvspan(12,13, alpha=0.2)
    f.savefig('./figures/correlation_tai0slip_gfp_aroundslippery_designedlibrary_first15.png', \
      dpi = 300, format='png', bbox_inches='tight', frameon=True)


with sns.axes_style(style='dark'):
    f, ax=plt.subplots(figsize=(4,3))
    ax.plot(tai0slip[[y in eventlist[:-1] for y in fsdf.fsevent[tai0slip.index]]].transpose().iloc[0:15].index.map(lambda x: 
        pearsonr(tai0slip[[y in eventlist[:-1] for y in fsdf.fsevent[tai0slip.index]]].transpose().loc[x], 
            fsdf[[y in eventlist[:-1] for y in fsdf.fsevent]&(fsdf.numberreadsm_wt>=20)&
                 (fsdf.peaksm_wt==1)&(fsdf.gfpm_wt<25)&(fsdf.subset=='fs_slipperyvar')&
                 (['Y' not in vs for vs in fsdf.varseq162])].gfpm_wt)[0]))
    ax.tick_params('y',colors=sns.xkcd_rgb['medium blue'])
    plt.axhline(y=0, linewidth=2, color=sns.xkcd_rgb['medium blue'])
    ax2=ax.twinx()
    ax2.plot(np.log10(tai0slip[[y in eventlist[:-1] for y in fsdf.fsevent[tai0slip.index]]].transpose().iloc[0:15].index.map(
            lambda x: pearsonr(tai0slip[[y in eventlist[:-1] for y in fsdf.fsevent[tai0slip.index]]].transpose().loc[x], 
            fsdf[[y in eventlist[:-1] for y in fsdf.fsevent]&(fsdf.numberreadsm_wt>=20)&
                 (fsdf.peaksm_wt==1)&(fsdf.gfpm_wt<25)&(fsdf.subset=='fs_slipperyvar')&
                 (['Y' not in vs for vs in fsdf.varseq162])].gfpm_wt)[1])), 
        color=sns.xkcd_rgb['green'], alpha=0.5)
    ax2.tick_params('y',colors=sns.xkcd_rgb['green'])
    plt.axvspan(12,13, alpha=0.2)
#    plt.axvspan(x=12, linewidth=4, alpha=0.3)
    f.savefig('./figures/correlation_tai0slip_gfp_aroundslippery_designedlibrary_eventlist_first15.png', \
      dpi = 300, format='png', bbox_inches='tight', frameon=True)




with sns.axes_style(style='dark'):
    f, ax=plt.subplots(figsize=(4,3))
    ax.plot(tai0slip[[y in nofslist for y in fsdf.fsevent[tai0slip.index]]].transpose().iloc[0:15].index.map(lambda x: 
        pearsonr(tai0slip[[y in nofslist for y in fsdf.fsevent[tai0slip.index]]].transpose().loc[x], 
            fsdf[[y in nofslist for y in fsdf.fsevent]&(fsdf.numberreadsm_wt>=20)&
                 (fsdf.peaksm_wt==1)&(fsdf.gfpm_wt<25)&(fsdf.subset=='fs_slipperyvar')&
                 (['Y' not in vs for vs in fsdf.varseq162])].gfpm_wt)[0]))
    ax.tick_params('y',colors=sns.xkcd_rgb['medium blue'])
    plt.axhline(y=0, linewidth=2, color=sns.xkcd_rgb['medium blue'])
    ax2=ax.twinx()
    ax2.plot(np.log10(tai0slip[[y in nofslist for y in fsdf.fsevent[tai0slip.index]]].transpose().iloc[0:15].index.map(
            lambda x: pearsonr(tai0slip[[y in nofslist for y in fsdf.fsevent[tai0slip.index]]].transpose().loc[x], 
            fsdf[[y in nofslist for y in fsdf.fsevent]&(fsdf.numberreadsm_wt>=20)&
                 (fsdf.peaksm_wt==1)&(fsdf.gfpm_wt<25)&(fsdf.subset=='fs_slipperyvar')&
                 (['Y' not in vs for vs in fsdf.varseq162])].gfpm_wt)[1])), 
        color=sns.xkcd_rgb['green'], alpha=0.5)
    ax2.tick_params('y',colors=sns.xkcd_rgb['green'])
    plt.axvspan(12,13, alpha=0.2)
    f.savefig('./figures/correlation_tai0slip_gfp_aroundslippery_designedlibrary_nofslist_first15.png', \
      dpi = 300, format='png', bbox_inches='tight', frameon=True)


with sns.axes_style(style='dark'):
    f, ax=plt.subplots(figsize=(4,3))
    ax.plot(tai0all[[y in eventlist[:-1] for y in fsdf.fsevent[tai0all.index]]].transpose().iloc[0:15].index.map(lambda x: 
        pearsonr(tai0all[[y in eventlist[:-1] for y in fsdf.fsevent[tai0all.index]]].transpose().loc[x], 
            fsdf[[y in eventlist[:-1] for y in fsdf.fsevent]&(fsdf.numberreadsm_wt>=20)&
                 (fsdf.peaksm_wt==1)&(fsdf.gfpm_wt<25)&
                 (['Y' not in vs for vs in fsdf.varseq162])].gfpm_wt)[0]))
    ax.tick_params('y',colors=sns.xkcd_rgb['medium blue'])
    plt.axhline(y=0, linewidth=2, color=sns.xkcd_rgb['medium blue'])
    ax2=ax.twinx()
    ax2.plot(np.log10(tai0all[[y in eventlist[:-1] for y in fsdf.fsevent[tai0all.index]]].transpose().iloc[0:15].index.map(
            lambda x: pearsonr(tai0all[[y in eventlist[:-1] for y in fsdf.fsevent[tai0all.index]]].transpose().loc[x], 
            fsdf[[y in eventlist[:-1] for y in fsdf.fsevent]&(fsdf.numberreadsm_wt>=20)&
                 (fsdf.peaksm_wt==1)&(fsdf.gfpm_wt<25)&
                 (['Y' not in vs for vs in fsdf.varseq162])].gfpm_wt)[1])), 
        color=sns.xkcd_rgb['green'], alpha=0.5)
    ax2.tick_params('y',colors=sns.xkcd_rgb['green'])
    plt.axvspan(12,13, alpha=0.2)
    f.savefig('./figures/correlation_tai0all_gfp_aroundslippery_designedlibrary_eventlist_first15.png', \
      dpi = 300, format='png', bbox_inches='tight', frameon=True)



with sns.axes_style(style='dark'):
    f, ax=plt.subplots(figsize=(4,3))
    ax.plot(tai0all.transpose().iloc[0:15].index.map(lambda x: pearsonr(tai0all.transpose().loc[x], 
            fsdf[(fsdf.numberreadsm_wt>=20)&(fsdf.peaksm_wt==1)&(fsdf.gfpm_wt<25)&(fsdf.library=='fs_designed')&(['Y' not in vs for vs in fsdf.varseq162])].gfpm_wt)[0]))
    ax.tick_params('y',colors=sns.xkcd_rgb['medium blue'])
    plt.axhline(y=0, linewidth=2, color=sns.xkcd_rgb['medium blue'])
    ax2=ax.twinx()
    ax2.plot(np.log10(tai0all.transpose().iloc[0:15].index.map(lambda x: pearsonr(tai0all.transpose().loc[x], 
            fsdf[(fsdf.numberreadsm_wt>=20)&(fsdf.peaksm_wt==1)&(fsdf.gfpm_wt<25)&(fsdf.library=='fs_designed')&(['Y' not in vs for vs in fsdf.varseq162])].gfpm_wt)[1])), 
        color=sns.xkcd_rgb['green'], alpha=0.5)
    ax2.tick_params('y',colors=sns.xkcd_rgb['green'])
    plt.axvspan(12,13, alpha=0.2)
    f.savefig('./figures/correlation_tai0all_gfp_aroundslippery_designedlibrary_full_first15.png', \
      dpi = 300, format='png', bbox_inches='tight', frameon=True)



#%% UPSTREAM

df=fsdf[(fsdf.subset=="fs_upstreamvar")]
df=df.dropna(axis=1, how='all')

df['subgroup']=df.changes.apply(lambda x: int(x[9]) if type(x)==str else 4)

df1=df[(df.subgroup==1)&(df.peaksm_wt==1)&(df.numberreadsm_wt>=20)]

df1['slippery']=df1.slippery.apply(lambda x: x.upper())

df1['upstrpos']=df1.changes.apply(lambda x: int(x.split('Number ')[1][0]))
df1['codon']=df1.changes.apply(lambda x: x[-3:])
df1['aa']=df1.codon.apply(lambda x: str(Seq(x).translate()))

df1['codonusage']=df1.codon.apply(lambda x: codonusagedict[x])
df1['tai']=df1.codon.apply(lambda x: tAIdict[x])

aacategory={
    'R':'charged',
    'H':'charged',
    'K':'charged',
    'D':'charged',
    'E':'charged',
    'S':'polar',
    'T':'polar',
    'N':'polar',
    'Q':'polar',
    'C':'polar',
    'Y':'polar',
    'G':'unpolar',
    'P':'unpolar',
    'A':'unpolar',
    'V':'unpolar',
    'I':'unpolar',
    'L':'unpolar',
    'M':'unpolar',
    'F':'unpolar',
    'W':'unpolar',
    }

aacategory={
    'R':'charged',
    'H':'polar',
    'K':'charged',
    'D':'charged',
    'E':'charged',
    'S':'polar',
    'T':'polar',
    'N':'polar',
    'Q':'polar',
    'C':'polar',
    'Y':'aromatic',
    'G':'unpolar',
    'P':'unpolar',
    'A':'unpolar',
    'V':'unpolar',
    'I':'unpolar',
    'L':'unpolar',
    'M':'polar',
    'F':'aromatic',
    'W':'aromatic',
    }

df1['lastcodonpos']=df1.index.map(lambda x: bool(df1.codon[x][-1]==df1.slippery[x][0]) if df1.upstrpos[x]==1 else True)
df1['iswt']=df1.index.map(lambda x: bool(df1.varseq162[x]==wtseqs[df1.fsevent[x]]))
df1['aaclass']=df1.aa.apply(lambda x: aacategory[x])


# plot the effect of codons preceding the slippery site on frameshifting

for event in eventlist:
    f=plt.figure(figsize=(18,8))
    ax=f.add_subplot(121)
    plt.plot(df1[(df1.fsevent==event)&(df1.gfpm_wt<25)&(df1.upstrpos==1)&(df1.aaclass=='unpolar')].tai, df1[(df1.fsevent==event)&(df1.gfpm_wt<25)&(df1.upstrpos==1)&(df1.aaclass=='unpolar')].gfpm_wt,marker='o',linewidth=0, color=sns.xkcd_rgb['light blue'])
    plt.plot(df1[(df1.fsevent==event)&(df1.gfpm_wt<25)&(df1.upstrpos==1)&(df1.aaclass=='polar')].tai, df1[(df1.fsevent==event)&(df1.gfpm_wt<25)&(df1.upstrpos==1)&(df1.aaclass=='polar')].gfpm_wt,marker='o',linewidth=0, color=sns.xkcd_rgb['green'])
    plt.plot(df1[(df1.fsevent==event)&(df1.gfpm_wt<25)&(df1.upstrpos==1)&(df1.aaclass=='charged')].tai, df1[(df1.fsevent==event)&(df1.gfpm_wt<25)&(df1.upstrpos==1)&(df1.aaclass=='charged')].gfpm_wt,marker='o',linewidth=0, color=sns.xkcd_rgb['red'])
#    plt.plot(df1[(df1.fsevent==event)&(df1.gfpm_wt<25)&(df1.upstrpos==1)&(df1.aaclass=='neg')].tai, df1[(df1.fsevent==event)&(df1.gfpm_wt<25)&(df1.upstrpos==1)&(df1.aaclass=='neg')].gfpm_wt,marker='o',linewidth=0, color=sns.xkcd_rgb['blue'])
    plt.plot(df1[(df1.fsevent==event)&(df1.gfpm_wt<25)&(df1.upstrpos==1)&(df1.lastcodonpos==True)].tai, df1[(df1.fsevent==event)&(df1.gfpm_wt<25)&(df1.upstrpos==1)&(df1.lastcodonpos==True)].gfpm_wt,marker='o',linewidth=0, color=sns.xkcd_rgb['black'], markersize=5)
    for label, x, y in zip(df1[(df1.fsevent==event)&(df1.gfpm_wt<25)&(df1.upstrpos==1)].codon, df1[(df.fsevent==event)&(df1.gfpm_wt<25)&(df1.upstrpos==1)].tai, df1[(df1.fsevent==event)&(df1.gfpm_wt<25)&(df1.upstrpos==1)].gfpm_wt):
        ax.annotate(label, xy=(x+0.02,y),fontsize=12, color='g')
    for label, x, y in zip(df1[(df1.fsevent==event)&(df1.gfpm_wt<25)&(df1.upstrpos==1)].aa, df1[(df.fsevent==event)&(df1.gfpm_wt<25)&(df1.upstrpos==1)].tai, df1[(df1.fsevent==event)&(df1.gfpm_wt<25)&(df1.upstrpos==1)].gfpm_wt):
        ax.annotate(label, xy=(x-0.03,y), fontsize=12, color='r')
    plt.xlim(-0.08, 1.08)
    plt.ylim(-0.2,wtvalues.at[event,'gfpm']+5)
    plt.axhline(y=wtvalues.at[event,'gfpm'], linewidth=2, color='gray', alpha=0.6)
    ax=f.add_subplot(122)
    plt.plot(df1[(df1.fsevent==event)&(df1.gfpm_wt<25)&(df1.upstrpos==2)&(df1.aaclass=='unpolar')].tai, df1[(df1.fsevent==event)&(df1.gfpm_wt<25)&(df1.upstrpos==2)&(df1.aaclass=='unpolar')].gfpm_wt,marker='o',linewidth=0, color=sns.xkcd_rgb['light blue'])
    plt.plot(df1[(df1.fsevent==event)&(df1.gfpm_wt<25)&(df1.upstrpos==2)&(df1.aaclass=='polar')].tai, df1[(df1.fsevent==event)&(df1.gfpm_wt<25)&(df1.upstrpos==2)&(df1.aaclass=='polar')].gfpm_wt,marker='o',linewidth=0, color=sns.xkcd_rgb['green'])
    plt.plot(df1[(df1.fsevent==event)&(df1.gfpm_wt<25)&(df1.upstrpos==2)&(df1.aaclass=='charged')].tai, df1[(df1.fsevent==event)&(df1.gfpm_wt<25)&(df1.upstrpos==2)&(df1.aaclass=='charged')].gfpm_wt,marker='o',linewidth=0, color=sns.xkcd_rgb['red'])
#    plt.plot(df1[(df1.fsevent==event)&(df1.gfpm_wt<25)&(df1.upstrpos==2)&(df1.aaclass=='neg')].tai, df1[(df1.fsevent==event)&(df1.gfpm_wt<25)&(df1.upstrpos==2)&(df1.aaclass=='neg')].gfpm_wt,marker='o',linewidth=0, color=sns.xkcd_rgb['blue'])
#    plt.plot(df1[(df1.fsevent==event)&(df1.gfpm_wt<25)&(df1.upstrpos==2)&(df1.lastcodonpos==True)].tai, df1[(df1.fsevent==event)&(df1.gfpm_wt<25)&(df1.upstrpos==2)&(df1.lastcodonpos==True)].gfpm_wt,marker='o',linewidth=0, color=sns.xkcd_rgb['black'], markersize=5)
    for label, x, y in zip(df1[(df.fsevent==event)&(df1.gfpm_wt<25)&(df1.upstrpos==2)].codon, df1[(df.fsevent==event)&(df1.gfpm_wt<25)&(df1.upstrpos==2)].tai, df1[(df1.fsevent==event)&(df1.gfpm_wt<25)&(df1.upstrpos==2)].gfpm_wt):
        ax.annotate(label, xy=(x+0.02,y),fontsize=12, color='g')
    for label, x, y in zip(df1[(df.fsevent==event)&(df1.gfpm_wt<25)&(df1.upstrpos==2)].aa, df1[(df.fsevent==event)&(df1.gfpm_wt<25)&(df1.upstrpos==2)].tai, df1[(df1.fsevent==event)&(df1.gfpm_wt<25)&(df1.upstrpos==2)].gfpm_wt):
        ax.annotate(label, xy=(x-0.03,y), fontsize=12, color='r')
    plt.ylim(-0.2,wtvalues.at[event,'gfpm']+5)
    plt.xlim(-0.08, 1.08)
    plt.axhline(y=wtvalues.at[event,'gfpm'], linewidth=2, color='gray', alpha=0.6)
    f.savefig('./figures/immediateupstream/upstreamvar_subgroup1_gfp_' + '_'.join(event.split(' ')) + '_percodon.png', \
      dpi = 300, format='png', bbox_inches='tight', frameon=True)


# all combined, by aa class

f=plt.figure(figsize=(3.5,3))
ax=sns.boxplot(data=df1[([x in eventlist for x in df1.fsevent])&(df1.gfpm_wt<25)&(df1.lastcodonpos==True)], x='aaclass',y='percgfpm_wt',hue='upstrpos', order=['unpolar','polar','charged'], hue_order=[2,1], palette=[sns.xkcd_rgb['light blue'], sns.xkcd_rgb['light green']])
ax=sns.swarmplot(data=df1[([x in eventlist for x in df1.fsevent])&(df1.gfpm_wt<25)&(df1.lastcodonpos==True)], x='aaclass',y='percgfpm_wt',hue='upstrpos', order=['unpolar','polar','charged'], hue_order=[2,1], split=True)
plt.axhline(y=100, color='grey', linewidth=2)
plt.ylim(-5,200)
ax.legend_.remove()
f.savefig('./figures/immediateupstream/uupstreamvar_subgroup1_by_aaclass_gfpperc_eventlist_lastcodonpositionthesame.png', \
  dpi = 300, format='png', bbox_inches='tight', frameon=True)


f=plt.figure(figsize=(3.5,3))
ax=sns.boxplot(data=df1[([x in eventlist for x in df1.fsevent])&(df1.gfpm_wt<25)&(df1.lastcodonpos==False)], x='aaclass',y='percgfpm_wt',hue='upstrpos', order=['unpolar','polar','charged'], hue_order=[2,1], palette=[sns.xkcd_rgb['light blue'], sns.xkcd_rgb['light green']])
ax=sns.swarmplot(data=df1[([x in eventlist for x in df1.fsevent])&(df1.gfpm_wt<25)&(df1.lastcodonpos==False)], x='aaclass',y='percgfpm_wt',hue='upstrpos', order=['unpolar','polar','charged'], hue_order=[2,1], split=True)
plt.axhline(y=100, color='grey', linewidth=2)
plt.ylim(-5,200)
ax.legend_.remove()
f.savefig('./figures/immediateupstream/upstreamvar_subgroup1_by_aaclass_gfpperc_eventlist_lastcodonpositiondifferent.png', \
  dpi = 300, format='png', bbox_inches='tight', frameon=True)




#######################

df2=df[df.subgroup==2]
#


f=plt.figure(figsize=(6,3))
ax=sns.swarmplot(data=df2[df2.gfpm_wt<25], x='fsevent',y='gfpm_wt', color=sns.xkcd_rgb['medium blue'])
wtvalues.loc[df2[df2.gfpm_wt<25].fsevent.unique(),'gfpm'].plot(marker='o',linewidth=0, color=sns.xkcd_rgb['light blue'])
ax.set_xticklabels(df2[df2.gfpm_wt<25].fsevent.unique(), rotation=45, horizontalalignment='right')
plt.xlim(-0.5,len(df2[df2.gfpm_wt<25].fsevent.unique())-0.5)
plt.ylim(0,10)
plt.axhspan(0, 1.3, color=sns.xkcd_rgb['medium blue'], alpha=0.2)
f.savefig('./figures/upstream/upstreamvar_subgroup2_gfpm_wt_overview_withwtvalues.png', \
  dpi = 300, format='png', bbox_inches='tight', frameon=True)


f=plt.figure(figsize=(6,3))
ax=sns.swarmplot(data=df2[df2.gfpp_wt<25], x='fsevent',y='gfpp_wt', color=sns.xkcd_rgb['medium blue'])
wtvalues.loc[df2[df2.gfpp_wt<25].fsevent.unique(),'gfpp'].plot(marker='o',linewidth=0, color=sns.xkcd_rgb['light blue'])
ax.set_xticklabels(df2[df2.gfpp_wt<25].fsevent.unique(), rotation=45, horizontalalignment='right')
plt.xlim(-0.5,len(df2[df2.gfpp_wt<25].fsevent.unique())-0.5)
plt.ylim(0,22)
plt.axhspan(0, 1.3, color=sns.xkcd_rgb['medium blue'], alpha=0.2)
f.savefig('./figures/upstream/upstreamvar_subgroup2_gfpp_wt_overview_withwtvalues.png', \
  dpi = 300, format='png', bbox_inches='tight', frameon=True)





#######################
# upstream recoding/structure variants

df234=df[(df.subgroup==2)|(df.subgroup==3)|(df.subgroup==4)]


plotorder=['subgroup 2: upstream sequence recoded randomly',
       'subgroup 3: upstream antarna inverse fold',
       'subgroup 4: upstream sequence no secondary structure',
       'subgroup 4: upstream sequence hairpin at 5prime end, from pos 0',
       'subgroup 4: upstream sequence hairpin at 3prime end, until pos 30']


f=plt.figure(figsize=(4,3))
ax=sns.boxplot(data=df234[[x in eventlist for x in df234.fsevent]&(df234.gfpm_wt<25)&(df234.peaksm_wt==1)], x='changes',y='percgfpm_wt', fliersize=0, color=sns.xkcd_rgb['light blue'], order=plotorder)
ax=sns.swarmplot(data=df234[[x in eventlist for x in df234.fsevent]&(df234.gfpm_wt<25)&(df234.peaksm_wt==1)], x='changes',y='percgfpm_wt', color=sns.xkcd_rgb['medium blue'], order=plotorder)
plt.axhline(y=100, linewidth=2, color='gray')
ax.set_xticklabels(['recoded randomly','structure-mimicking','no secondary structure','hairpin 5prime','hairpin 3prime'], rotation=45, horizontalalignment='right')
plt.ylim(0,250)
f.savefig('./figures/upstream/upstreamvar_subgroup234_gfpperc_eventlist_combined_boxplot.png', \
  dpi = 300, format='png', bbox_inches='tight', frameon=True)

f=plt.figure(figsize=(4,3))
ax=sns.boxplot(data=df234[[x in nofslist for x in df234.fsevent]&(df234.gfpm_wt<25)&(df234.peaksm_wt==1)], x='changes',y='percgfpm_wt', fliersize=0, color=sns.xkcd_rgb['light blue'], order=plotorder)
ax=sns.swarmplot(data=df234[[x in nofslist for x in df234.fsevent]&(df234.gfpm_wt<25)&(df234.peaksm_wt==1)], x='changes',y='percgfpm_wt', color=sns.xkcd_rgb['medium blue'], order=plotorder)
plt.axhline(y=100, linewidth=2, color='gray')
ax.set_xticklabels(['recoded randomly','structure-mimicking','no secondary structure','hairpin 5prime','hairpin 3prime'], rotation=45, horizontalalignment='right')
plt.ylim(0,550)
f.savefig('./figures/upstream/upstreamvar_subgroup234_gfpperc_nofslist_combined_boxplot.png', \
  dpi = 300, format='png', bbox_inches='tight', frameon=True)


f=plt.figure(figsize=(4,3))
ax=sns.boxplot(data=df234[[x in plist for x in df234.fsevent]&(df234.gfpp_wt<25)&(df234.peaksp_wt==1)], x='changes',y='percgfpp_wt', fliersize=0, color=sns.xkcd_rgb['light blue'], order=plotorder)
ax=sns.swarmplot(data=df234[[x in plist for x in df234.fsevent]&(df234.gfpp_wt<25)&(df234.peaksp_wt==1)], x='changes',y='percgfpp_wt', color=sns.xkcd_rgb['medium blue'], order=plotorder)
plt.axhline(y=100, linewidth=2, color='gray')
ax.set_xticklabels(['recoded randomly','structure-mimicking','no secondary structure','hairpin 5prime','hairpin 3prime'], rotation=45, horizontalalignment='right')
plt.ylim(0,110)
f.savefig('./figures/upstream/upstreamvar_subgroup234_gfpperc_plist_combined_boxplot.png', \
  dpi = 300, format='png', bbox_inches='tight', frameon=True)


f=plt.figure(figsize=(4,3))
ax=sns.boxplot(data=df234[[x in nofslist for x in df234.fsevent]&(df234.gfpm_wt<25)&(df234.peaksm_wt==1)], 
                          x='changes',y='gfpm_wt', fliersize=0, color=sns.xkcd_rgb['light blue'], order=plotorder)
ax=sns.swarmplot(data=df234[[x in nofslist for x in df234.fsevent]&(df234.gfpm_wt<25)&(df234.peaksm_wt==1)], 
                            x='changes',y='gfpm_wt', color=sns.xkcd_rgb['medium blue'], order=plotorder)
ax.set_xticklabels(['recoded randomly','structure-mimicking','no secondary structure','hairpin 5prime','hairpin 3prime'], rotation=45, horizontalalignment='right')
plt.ylim(0,9)
plt.axhspan(-0.2, 1.3, color=sns.xkcd_rgb['medium blue'], alpha=0.2)
f.savefig('./figures/upstream/upstreamvar_subgroup234_gfp_nofslist_combined_boxplot.png', \
  dpi = 300, format='png', bbox_inches='tight', frameon=True)


#nofslist.remove('herpes simplex - thymidine kinase')
f=plt.figure(figsize=(5,3))
#ax=sns.boxplot(data=df234[[x in nofslist for x in df234.fsevent]&(df234.gfpm_wt<25)&(df234.peaksm_wt==1)], 
#                          x='changes',y='percgfpm_wt', hue='fsevent', fliersize=0, order=plotorder)
ax=sns.swarmplot(data=df234[[x in nofslist for x in df234.fsevent]&(df234.gfpm_wt<25)&(df234.peaksm_wt==1)], 
                         x='changes',y='percgfpm_wt', hue='fsevent', split=True, order=plotorder)
plt.axhline(y=100, linewidth=2, color='gray')
ax.set_xticklabels(['recoded randomly','structure-mimicking','no secondary structure','hairpin 5prime','hairpin 3prime'], rotation=45, horizontalalignment='right')
plt.ylim(-5,805)
plt.legend(bbox_to_anchor=(2,1))
f.savefig('./figures/upstream/upstreamvar_subgroup234_gfpperc_nofslist_swarmplot.png', \
  dpi = 300, format='png', bbox_inches='tight', frameon=True)

df234[[x in nofslist for x in df234.fsevent]&(df234.gfpm_wt<25)&(df234.peaksm_wt==1)].pivot_table(
        index='changes',values='percgfpm_wt', columns='fsevent',aggfunc=np.count_nonzero)

f=plt.figure(figsize=(5,3))
#ax=sns.boxplot(data=df234[[x in nofslist for x in df234.fsevent]&(df234.gfpm_wt<25)&(df234.peaksm_wt==1)], 
#                          x='changes',y='gfpm_wt', hue='fsevent', fliersize=0, order=plotorder)
ax=sns.swarmplot(data=df234[[x in nofslist for x in df234.fsevent]&(df234.gfpm_wt<25)&(df234.peaksm_wt==1)], 
                         x='changes',y='gfpm_wt', hue='fsevent', split=True, order=plotorder)
ax.set_xticklabels(['recoded randomly','structure-mimicking','no secondary structure','hairpin 5prime','hairpin 3prime'], rotation=45, horizontalalignment='right')
plt.ylim(-0.2,9)
plt.legend(bbox_to_anchor=(2,1))
plt.axhspan(-0.2, 1.3, color=sns.xkcd_rgb['medium blue'], alpha=0.2)
f.savefig('./figures/upstream/upstreamvar_subgroup234_gfp_nofslist_swarmplot.png', \
  dpi = 300, format='png', bbox_inches='tight', frameon=True)


#######################
# replace upstream region from the 5' end

df5=df[(df.subgroup==5)]

altstarts=['gaagctgctgcaagagaagctgcagctagggaggctgcagctagggaggctgctgcaaga',
 'CTACAACACAAGCGTTCGTTACCCTGCCAC',
 'TGCCACCGTGCCATTGAGTCCATTTGTCCG',
 'AGATGCAGCACCGACCCCTCCCACCACTCG']


def findendostartpos(x):
    endo=fsdf[(fsdf.fsevent==df5.fsevent[x])&(fsdf.subset=='fsmain')].varseq162.unique()[0]
    return difflib.SequenceMatcher(a=df.varseq162[x], b=endo).get_matching_blocks()[-2][1]-38


df5['endostartpos']=df5.index.map(lambda x: findendostartpos(x))



f=plt.figure(figsize=(3,3))
ax=sns.regplot(data=df5[(df5.gfpp_wt<25)&(df5.fsevent=='human - OAZ1')], x='endostartpos',y='percgfpp_wt', fit_reg=False, color=sns.xkcd_rgb['medium blue'])
ax=sns.regplot(data=df5[(df5.gfpm_wt<25)&(df5.fsevent=='HIV HXB2')], x='endostartpos',y='percgfpm_wt', fit_reg=False, color=sns.xkcd_rgb['medium green'])
plt.axhline(y=100, color='gray', linewidth=2)
ax.set_xticks([-40,-30,-20,-10,0])
ax.set_yticks([0,50,100,150,200])
plt.ylim(-5,205)
f.savefig('./figures/upstream/upstream_df5_oaz1vshiv.png', \
  dpi = 300, format='png', bbox_inches='tight', frameon=True)



###############

#%%
df=fsdf[(fsdf.subset=="fs_downstrvar")]
df=df.dropna(axis=1, how='all')
df['subgroup']=df.changes.apply(lambda x: int(x[9]) if type(x)==str else 5)


for i in df[df.subgroup==1].index:
    df.loc[i,'ntreplaced']=int(df.loc[i,'changes'].split(' ')[-1][:-3])

for i in df[df.subgroup==1].index:
    df.loc[i,'ntkept']=120-df.loc[i,'ntreplaced']


pkiss=pd.read_table('../additional/pKiss_structures_MFE_downstrvar.txt', header=None)

pkissdf=pd.DataFrame()
for item in pkiss.values:
    if item[0][0]=='>':
        var=int(item[0][1:])
    elif item[0][0]=='-':
        items=item[0].split(' ')
        pkissdf.loc[var, 'mfe']=float(items[0])
        pkissdf.loc[var, 'structure']=items[-1]

df=df.join(pkissdf)


f=plt.figure(figsize=(3,3))
plt.scatter(df.dg_downstream, df.mfe, alpha=0.1)
plt.xlim(-71,1)
plt.ylim(-81,1)
plt.xlabel('MFE (Vienna RNA)')
plt.ylabel('MFE (pKiss)')
plt.title('Pearson r=' + '{:.2g}'.format(scipy.stats.pearsonr(df.dg_downstream, df.mfe)[0]) + 
          ', p=' + '{:.2g}'.format(scipy.stats.pearsonr(df.dg_downstream, df.mfe)[1]), fontsize=12)
f.savefig('./figures/downstream/MFE_Vienna_vs_pKiss.png', \
          dpi = 300, format='png', bbox_inches='tight', frameon=True)


# replace downstream region from the 3' end

f=plt.figure(figsize=(6,3))
sns.pointplot(df[(df.peaksm_wt==1)&(df.subgroup==1)&(df.percgfpm_wt<200)&(df.fsevent=='HIV HXB2')].ntkept, df[(df.peaksm_wt==1)
&(df.subgroup==1)&(df.gfpm_wt<25)&(df.fsevent=='HIV HXB2')].percgfpm_wt, 
                 color=eventcolor['HIV HXB2'], scale=0.5, linewidth=2)
sns.pointplot(df[(df.peaksm_wt==1)&(df.subgroup==1)&(df.percgfpm_wt<200)&(df.fsevent=='SARS coronavirus')].ntkept, df[(df.peaksm_wt==1)&
                 (df.subgroup==1)&(df.gfpm_wt<25)&(df.fsevent=='SARS coronavirus')].percgfpm_wt, 
color=eventcolor['SARS coronavirus'], scale=0.5, linewidth=2)
sns.pointplot(df[(df.peaksm_wt==1)&(df.subgroup==1)&(df.percgfpm_wt<200)&(df.fsevent=='simian srv1 retrovirus')].ntkept, df[(df.peaksm_wt==1)&
                 (df.subgroup==1)&(df.gfpm_wt<25)&(df.fsevent=='simian srv1 retrovirus')].percgfpm_wt, 
color=eventcolor['simian srv1 retrovirus'], scale=0.5, linewidth=2)
sns.pointplot(df[(df.peaksm_wt==1)&(df.subgroup==1)&(df.percgfpm_wt<200)&(df.fsevent=='SIVmac239')].ntkept, df[(df.peaksm_wt==1)&
                 (df.subgroup==1)&(df.gfpm_wt<25)&(df.fsevent=='SIVmac239')].percgfpm_wt, 
color=eventcolor['SIVmac239'], scale=0.5, linewidth=2)

plt.axhline(y=100, color='grey', linewidth=2)
plt.ylim(-5,200)
plt.ylabel('GFP fluorescence')
plt.xlabel('nucleotides after slippery site unchanged')
f.savefig('./figures/downstream/downstreamvar_gfp_subgroup1_hivsarsherv.png', \
          dpi = 300, format='png', bbox_inches='tight', frameon=True)

f=plt.figure(figsize=(6,3))
for i in eventcolor.keys():
    if i!='PLRV luteovirus':
        sns.pointplot(df[(df.peaksm_wt==1)&(df.subgroup==1)&(df.percgfpm_wt<200)&(df.fsevent==i)].ntkept, df[(df.peaksm_wt==1)
        &(df.subgroup==1)&(df.gfpm_wt<25)&(df.fsevent==i)].percgfpm_wt, 
                     color=eventcolor[i], scale=0.5, linewidth=2)
plt.ylim(-5,200)
plt.axhline(y=100, color='grey', linewidth=2)
plt.ylabel('GFP fluorescence')
plt.xlabel('nucleotides after slippery site unchanged')
f.savefig('./figures/downstream/downstreamvar_gfp_subgroup1_allexceptplrv.png', \
          dpi = 300, format='png', bbox_inches='tight', frameon=True)

df[(df.peaksm_wt==1)&(df.subgroup==1)&(df.percgfpm_wt<200)&(df.fsevent==i)].pivot_table(
    index='fsevent',columns='ntkept',values='percgfpm_wt', aggfunc=np.count_nonzero).to_csv('./n_numbers/Fig3b.csv')
# Replace native with sequence that is predicted to fold into the same secondary structure

df[df.subgroup==3]['event']=df[df.subgroup==3].changes.apply(lambda x: x.split(' ')[1].split('_')[0])
df['event']=df[df.subgroup==3].changes.apply(lambda x: x.split(' ')[-1].split('_')[0])
df['method']=df[df.subgroup==3].changes.apply(lambda x: x.split(' ')[-2])


f=plt.figure(figsize=(4,3))
ax=sns.swarmplot(data=df[[x in eventlist for x in df.fsevent]&(df.subgroup==3)&(df.gfpm_wt<25)&(df.laststopinframem1<12)], x='fsevent',y='gfpm_wt', \
                         order=eventlist, color=sns.xkcd_rgb['medium blue'])
wtvalues.loc[eventlist,'gfpm'].plot(marker='v',linewidth=0, color=sns.xkcd_rgb['dark red'], alpha=0.9)
ax.set_xticklabels(eventlist, rotation=45,horizontalalignment='right')
plt.xlim(-0.5,len(eventlist)-0.5)
plt.xlabel('')
plt.ylim(-0.2,12)
plt.axhspan(-0.2, 1.3, color=sns.xkcd_rgb['medium blue'], alpha=0.2)
plt.ylabel('% GFP fluorescence')
f.savefig('./figures/downstream/downstreamvar_subgroup3_eventlist.png', \
          dpi = 300, format='png', bbox_inches='tight', frameon=True)



f=plt.figure(figsize=(3,3))
ax=sns.swarmplot(data=df[[x in nofslist for x in df.fsevent]&(df.subgroup==3)&(df.gfpm_wt<25)&(df.laststopinframem1<12)], x='fsevent',y='gfpm_wt', \
                         order=nofslist, color=sns.xkcd_rgb['medium blue'])
wtvalues.loc[nofslist,'gfpm'].plot(marker='v',linewidth=0, color=sns.xkcd_rgb['dark red'], alpha=0.6)
ax.set_xticklabels(nofslist, rotation=45,horizontalalignment='right')
plt.xlim(-0.5,len(nofslist)-0.5)
plt.xlabel('')
plt.ylim(-0.2,12)
plt.axhspan(-0.2, 1.3, color=sns.xkcd_rgb['medium blue'], alpha=0.2)
plt.ylabel('% GFP fluorescence')
f.savefig('./figures/downstream/downstreamvar_subgroup3_nofslist.png', \
          dpi = 300, format='png', bbox_inches='tight', frameon=True)


f=plt.figure(figsize=(3,3))
ax=sns.swarmplot(data=df[[x in plist for x in df.fsevent]&(df.subgroup==3)&(df.gfpp_wt<35)], x='fsevent',y='gfpp_wt', \
                         order=plist, color=sns.xkcd_rgb['medium blue'])
wtvalues.loc[plist,'gfpp'].plot(marker='v',linewidth=0, color=sns.xkcd_rgb['dark red'], alpha=0.6)
ax.set_xticklabels(plist, rotation=45,horizontalalignment='right')
plt.xlim(-0.5,len(plist)-0.5)
plt.xlabel('')
plt.ylim(-0.2,25)
plt.axhspan(-0.2, 1.3, color=sns.xkcd_rgb['medium blue'], alpha=0.2)
plt.ylabel('% GFP fluorescence')
f.savefig('./figures/downstream/downstreamvar_subgroup3_plist_p1.png', \
          dpi = 300, format='png', bbox_inches='tight', frameon=True)


f=plt.figure(figsize=(6,3))
ax=sns.swarmplot(data=df[[x in list(eventlist)+list(nofslist) for x in df.fsevent]&(df.subgroup==3)&(df.gfpm_wt<25)], x='fsevent',y='gfpm_wt', \
                          hue='method',
                         order=wtvalues.loc[list(eventlist)+list(nofslist)].sort_values(by='gfpm', ascending=False).index, palette=[sns.xkcd_rgb['medium blue'],sns.xkcd_rgb['dark blue']])
wtvalues.loc[list(eventlist)+list(nofslist)].sort_values(by='gfpm', ascending=False).gfpm.plot(marker='v',linewidth=0, color=sns.xkcd_rgb['dark red'], alpha=0.6)
ax.set_xticklabels(wtvalues.loc[list(eventlist)+list(nofslist)].sort_values(by='gfpm', ascending=False).index, rotation=45,horizontalalignment='right')
plt.xlim(-0.5,len(list(eventlist)+list(nofslist))-0.5)
plt.xlabel('')
plt.ylim(-0.2,10.2)
plt.axhspan(-0.2, 1.3, color=sns.xkcd_rgb['medium blue'], alpha=0.2)
plt.ylabel('% GFP fluorescence')
ax.legend_.remove()
f.savefig('./figures/downstream/downstreamvar_subgroup3_eventlist_fslist.png', \
          dpi = 300, format='png', bbox_inches='tight', frameon=True)



df[[x in list(eventlist)+list(nofslist) for x in df.fsevent]&(df.subgroup==3)&(df.gfpm_wt<25)].pivot_table(
    index='fsevent',columns='method',values='gfpm_wt', aggfunc=np.count_nonzero).to_csv('./n_numbers/Fig3f.csv')



############# scanning mutagenesis of downstream region

df2=df[(df.subgroup==2)]

df2['nt']=df2.changes.apply(lambda x: x.split(' at position ')[0][-1])
df2['pos']=df2.changes.apply(lambda x: int(x.split(' at position ')[1].split(' ')[0]))

def fraction_background(event):
    try:
        return float(len(df2[(df2.fsevent==event)&(df2.gfpm_wt<1.3)&(df2.laststopinframem1<12)]))/\
        len(df2[(df2.fsevent==event)&(df2.gfpm_wt<25)&(df2.laststopinframem1<12)])
    except:
        return np.nan
    
def fraction_background_p1(event):
    try:
        return float(len(df2[(df2.fsevent==event)&(df2.gfpp_wt<1.3)&(df2.laststopinframep1<12)]))/\
        len(df2[(df2.fsevent==event)&(df2.gfpp_wt<25)&(df2.laststopinframep1<12)])
    except:
        return np.nan

def fraction_half(event):
    try:
        return float(len(df2[(df2.fsevent==event)&(df2.gfpm_wt<25)&(df2.percgfpm_wt<50)&(df2.laststopinframem1<12)]))/\
        len(df2[(df2.fsevent==event)&(df2.gfpm_wt<25)&(df2.laststopinframem1<12)])
    except:
        return np.nan

def fraction_half_p1(event):
    try:
        return float(len(df2[(df2.fsevent==event)&(df2.gfpp_wt<25)&(df2.percgfpp_wt<50)&(df2.laststopinframep1<12)]))/\
        len(df2[(df2.fsevent==event)&(df2.gfpp_wt<25)&(df2.laststopinframep1<12)])
    except:
        return np.nan

positionsaffecting=pd.DataFrame()
for i in df2.fsevent.unique():
    positionsaffecting.loc[i,'background']=fraction_background(i)
    positionsaffecting.loc[i,'background_p1']=fraction_background_p1(i)
    positionsaffecting.loc[i,'half']=fraction_half(i)
    positionsaffecting.loc[i,'half_p1']=fraction_half_p1(i)
    positionsaffecting.loc[i,'wtvalues']=wtvalues.loc[i,'gfpm']
    positionsaffecting.loc[i,'wtvalues_p1']=wtvalues.loc[i,'gfpp']


with sns.axes_style('dark'):
    ax=positionsaffecting.loc[eventlist,['background','half']].dropna().sort_values(by='half').plot(kind='barh', figsize=(3,4))
    ax.xaxis.grid()
    plt.xlim(0,0.8)
    plt.xticks([0,0.2,0.4,0.6,0.8,1])
    plt.xlabel('fraction of point mutations')
    ax.legend_.remove()
    f=ax.get_figure()
    f.savefig('./figures/downstream/downstreamvar_subset2_eventlist_fractionofpointmutationsaffectingPRF.png', \
          dpi = 300, format='png', bbox_inches='tight', frameon=True)


def n_background(event):
    try:
        return (len(df2[(df2.fsevent==event)&(df2.gfpm_wt<1.3)&(df2.laststopinframem1<12)])),\
        len(df2[(df2.fsevent==event)&(df2.gfpm_wt<25)&(df2.laststopinframem1<12)])
    except:
        return np.nan
    
def n_background_p1(event):
    try:
        return (len(df2[(df2.fsevent==event)&(df2.gfpp_wt<1.3)&(df2.laststopinframep1<12)])),\
        len(df2[(df2.fsevent==event)&(df2.gfpp_wt<25)&(df2.laststopinframep1<12)])
    except:
        return np.nan

def n_half(event):
    try:
        return (len(df2[(df2.fsevent==event)&(df2.gfpm_wt<25)&(df2.percgfpm_wt<50)&(df2.laststopinframem1<12)])),\
        len(df2[(df2.fsevent==event)&(df2.gfpm_wt<25)&(df2.laststopinframem1<12)])
    except:
        return np.nan

def n_half_p1(event):
    try:
        return (len(df2[(df2.fsevent==event)&(df2.gfpp_wt<25)&(df2.percgfpp_wt<50)&(df2.laststopinframep1<12)])),\
        len(df2[(df2.fsevent==event)&(df2.gfpp_wt<25)&(df2.laststopinframep1<12)])
    except:
        return np.nan

npositionsaffecting=pd.DataFrame()
for i in df2.fsevent.unique():
    npositionsaffecting.loc[i,'background'], npositionsaffecting.loc[i,'total']=n_background(i)
    npositionsaffecting.loc[i,'background_p1'], npositionsaffecting.loc[i,'total_p1']=n_background_p1(i)
    npositionsaffecting.loc[i,'half'], npositionsaffecting.loc[i,'total_half']=n_half(i)
    npositionsaffecting.loc[i,'half_p1'], npositionsaffecting.loc[i,'total_half_p1']=n_half_p1(i)
    
npositionsaffecting.to_csv('./n_numbers/Fig3c.csv')

'''
                            background      half
SARS coronavirus              0.333333  0.529412
SIVmac239                     0.040816  0.061224
simian srv1 retrovirus        0.382979  0.510638
PLRV luteovirus               0.574468  0.319149
human - HERV-K10              0.277778  0.370370
human T-lymphotropic virus    0.344828  0.034483
HIV HXB2                      0.117647  0.137255
west nile virus               0.260870  0.347826

'''

with sns.axes_style('dark'):
    ax=positionsaffecting.loc[plist,['background_p1','half_p1']].sort_values(by='half_p1').plot(kind='barh', figsize=(3,1.7))
    ax.xaxis.grid()
    plt.xlim(0,1)
    plt.xticks([0,0.2,0.4,0.6,0.8,1])
    plt.xlabel('fraction of point mutations')
    ax.legend_.remove()
    f=ax.get_figure()
    f.savefig('./figures/downstream/downstreamvar_subset2_plist_fractionofpointmutationsaffectingPRF.png', \
          dpi = 300, format='png', bbox_inches='tight', frameon=True)
'''
                   background_p1   half_p1
human - OAZ1            0.019231  0.461538
influenza a virus       0.400000  0.500000
HCV - F protein         0.666667  0.777778

'''
# Compare downstream regions with single point mutations

eventstotest=['human - HERV-K10',
 'simian srv1 retrovirus',
 'SARS coronavirus',
 'HIV HXB2']

f=plt.figure(figsize=(9,3))
ax=f.add_subplot(111)
sns.pointplot(data=df2[[x in eventstotest for x in df2.fsevent.astype(str)]&
            (df2.gfpm_wt<25)&(df2.percgfpm_wt<200)], x='pos',y='percgfpm_wt', hue='fsevent',
    hue_order=eventstotest,
    palette=[eventcolor[x] for x in eventstotest], 
            scale=0.3, linewidth=2, errwidth=2, dodge=0.3)
ax.legend_.remove()
plt.ylim(0,200)
plt.axhline(y=100, linewidth=2, color='grey')
f.savefig('./figures/downstream/downstreamvar_subset2_mostevents.png', \
          dpi = 300, format='png', bbox_inches='tight', frameon=True)

eventstotest=['human - HERV-K10','PLRV luteovirus',
 'simian srv1 retrovirus']

f=plt.figure(figsize=(9,3))
ax=f.add_subplot(111)
sns.pointplot(data=df2[[x in eventstotest for x in df2.fsevent.astype(str)]&
            (df2.gfpm_wt<25)&(df2.percgfpm_wt<200)], x='pos',y='percgfpm_wt', hue='fsevent',
    hue_order=eventstotest,
    palette=[eventcolor[x] for x in eventstotest], 
            scale=0.3, linewidth=2, errwidth=2, dodge=0.3)
ax.legend_.remove()
plt.ylim(0,200)
plt.axhline(y=100, linewidth=2, color='grey')
f.savefig('./figures/downstream/downstreamvar_subset2_hervplrvsrv1.png', \
          dpi = 300, format='png', bbox_inches='tight', frameon=True)

plistcolors={'human - OAZ1':'green', 'HCV - F protein':'magenta','influenza a virus':'orange'}

f=plt.figure(figsize=(9,3))
ax=f.add_subplot(111)
sns.pointplot(data=df2[[x in plist for x in df2.fsevent.astype(str)]&
            (df2.gfpp_wt<25)&(df2.percgfpp_wt<200)], x='pos',y='percgfpp_wt', hue='fsevent',
    hue_order=plist,
    palette=[plistcolors[x] for x in plist], 
            scale=0.3, linewidth=2, errwidth=2, dodge=0.3)
ax.legend_.remove()
plt.ylim(0,200)
plt.axhline(y=100, linewidth=2, color='grey')
f.savefig('./figures/downstream/downstreamvar_subset2_plist.png', \
          dpi = 300, format='png', bbox_inches='tight', frameon=True)


# For PRF events separately

for event in df2.fsevent.unique():
    f=plt.figure(figsize=(len(df2[(df2.gfpm_wt<25)&(df2.percgfpm_wt<200)&(df.fsevent==event)].pos.unique())/2,3))
    sns.swarmplot(data=df2[(df2.gfpm_wt<25)&(df2.percgfpm_wt<200)&(df.fsevent==event)], 
                x='pos', y='gfpm_wt',hue='nt', hue_order=['A','C'], 
                order=np.sort(df2[(df2.gfpm_wt<25)&(df2.percgfpm_wt<200)&(df.fsevent==event)].pos.unique()),
                    palette=[sns.xkcd_rgb['green'],sns.xkcd_rgb['dark green']])
    plt.axhline(y=wtvalues.loc[event,'gfpm'], color='grey', linewidth=2)
    plt.ylabel('% GFP intensity')
    plt.xlabel('position relative to FS site')
#    plt.ylim(-3, df2[(df2.gfp<25)&(df.fsevent==event)&(df2.gfpperc<200)].gfpperc.max()/10*10+10)
    plt.legend(bbox_to_anchor=(1.1,1.5))
    plt.title(event, fontsize=14)
    plt.ylim(0)
    plt.axhspan(0,1.3,color=sns.xkcd_rgb['medium blue'], alpha=0.2)
    f.savefig('./figures/downstream/downstreamvar_subgroup2_gfp_'+'_'.join(event) +'.png', \
              dpi = 300, format='png', bbox_inches='tight', frameon=True)

for event in plist:
    f=plt.figure(figsize=(len(df2[(df2.gfpp_wt<25)&(df2.percgfpp_wt<200)&(df.fsevent==event)].pos.unique())/2,3))
    sns.swarmplot(data=df2[(df2.gfpp_wt<25)&(df2.percgfpp_wt<200)&(df.fsevent==event)], x='pos', y='gfpp_wt',hue='nt', palette=[sns.xkcd_rgb['green'],sns.xkcd_rgb['dark green']])
    ax.set_xticklabels(np.sort(df[(df.subgroup==1)&(df.gfpp_wt<25)&(df2.percgfpp_wt<200)&(df.fsevent==event)].ntreplaced.unique().astype(int)))
    plt.axhline(y=wtvalues.loc[event,'gfpp'], color='grey', linewidth=2)
    plt.ylabel('% wild-type GFP intensity')
    plt.xlabel('position relative to FS site')
    plt.ylim(0)
    plt.axhspan(0,1.3,color=sns.xkcd_rgb['medium blue'], alpha=0.2)
    plt.legend(bbox_to_anchor=(1.1,1.5))
    plt.title(event, fontsize=14)
    f.savefig('./figures/downstream/downstreamvar_subgroup2_plist_gfpp_'+'_'.join(event) +'.png', \
              dpi = 300, format='png', bbox_inches='tight', frameon=True)


#### Secondary structure variants

df4=df[df.subgroup==4]

df4['origin']=df4.changes.apply(lambda x: x.split(' ')[2])
df4['structure']=df4.changes.apply(lambda x: x.split(' ')[-1])
df4['GUpairing']=df4.changes.apply(lambda x: 'no GU pairing' not in x)

f=plt.figure(figsize=(6,3))
ax=sns.swarmplot(data=df4[[x in eventlist[1:] for x in df4.fsevent]&(df4.gfpm_wt<25)&(df4.origin=='hiv')], 
                          x='fsevent', y='gfpm_wt', hue='dg_downstream',\
                          order=eventlist[1:], palette='Blues_d')
wtvalues.loc[eventlist[1:],'gfpm'].plot(marker='v',linewidth=0, color=sns.xkcd_rgb['dark red'], alpha=0.9)
ax.set_xticklabels(eventlist[1:], rotation=45,horizontalalignment='right')
plt.xlim(-0.5,len(eventlist[1:])-0.5)
plt.ylim(-0.3,12.3)
plt.ylabel('% GFP fluorescence')
plt.axhspan(-0.3, 1.3, color=sns.xkcd_rgb['medium blue'], alpha=0.2)
ax.legend_.remove()
f.savefig('./figures/downstream/downstreamvar_subgroup4_gfp_eventlist_hiv.png', \
          dpi = 300, format='png', bbox_inches='tight', frameon=True)

f=plt.figure(figsize=(9,3))
ax=sns.swarmplot(data=df4[(df4.gfpm_wt<25)&(df4.origin=='hiv')], 
                          x='fsevent', y='gfpm_wt', hue='dg_downstream',\
                          order=df4[(df4.gfpm_wt<25)&(df4.origin=='hiv')].fsevent.unique(),
                           palette='Blues_d')
wtvalues.loc[df4[(df4.gfpm_wt<25)&(df4.origin=='hiv')].fsevent.unique(),'gfpm'].plot(marker='v',linewidth=0, color=sns.xkcd_rgb['dark red'], alpha=0.9)
ax.set_xticklabels(df4[(df4.gfpm_wt<25)&(df4.origin=='hiv')].fsevent.unique(), rotation=45,horizontalalignment='right')
plt.xlim(-0.5,len(df4[(df4.gfpm_wt<25)&(df4.origin=='hiv')].fsevent.unique())-0.5)
plt.ylim(-0.3,12.3)
plt.ylabel('% GFP fluorescence')
plt.axhspan(-0.3, 1.3, color=sns.xkcd_rgb['medium blue'], alpha=0.2)
ax.legend_.remove()
f.savefig('./figures/downstream/downstreamvar_subgroup4_gfp_allm1_hiv_swarmplot_withdeltag.png', \
          dpi = 300, format='png', bbox_inches='tight', frameon=True)

f=plt.figure(figsize=(9,3))
ax=sns.swarmplot(data=df4[(df4.gfpm_wt<25)&(df4.origin=='srvl')], 
                          x='fsevent', y='gfpm_wt', hue='dg_downstream',\
                          order=df4[(df4.gfpm_wt<25)&(df4.origin=='srvl')].fsevent.unique(),
                           palette='Blues_d')
wtvalues.loc[df4[(df4.gfpm_wt<25)&(df4.origin=='srvl')].fsevent.unique(),'gfpm'].plot(marker='v',linewidth=0, color=sns.xkcd_rgb['dark red'], alpha=0.9)
ax.set_xticklabels(df4[(df4.gfpm_wt<25)&(df4.origin=='srvl')].fsevent.unique(), rotation=45,horizontalalignment='right')
plt.xlim(-0.5,len(df4[(df4.gfpm_wt<25)&(df4.origin=='srvl')].fsevent.unique())-0.5)
plt.ylim(-0.3,12.3)
plt.ylabel('% GFP fluorescence')
plt.axhspan(-0.3, 1.3, color=sns.xkcd_rgb['medium blue'], alpha=0.2)
ax.legend_.remove()
f.savefig('./figures/downstream/downstreamvar_subgroup4_gfp_allm1_srvl_swarmplot_withdeltag.png', \
          dpi = 300, format='png', bbox_inches='tight', frameon=True)


f=plt.figure(figsize=(4,3))
ax=sns.swarmplot(data=df4[(df4.gfpm_wt<25)&(df4.origin=='nsp2f')], 
                          x='fsevent', y='gfpp_wt', hue='dg_downstream',\
                          order=df4[(df4.gfpp_wt<25)&(df4.origin=='nsp2f')].fsevent.unique(),
                           palette='Blues_d')
wtvalues.loc[df4[(df4.gfpp_wt<25)&(df4.origin=='nsp2f')].fsevent.unique(),'gfpp'].plot(marker='v',linewidth=0, color=sns.xkcd_rgb['dark red'], alpha=0.9)
ax.set_xticklabels(df4[(df4.gfpp_wt<25)&(df4.origin=='nsp2f')].fsevent.unique(), rotation=45,horizontalalignment='right')
plt.xlim(-0.5,len(df4[(df4.gfpp_wt<25)&(df4.origin=='nsp2f')].fsevent.unique())-0.5)
plt.ylim(-0.3,25)
plt.ylabel('% GFP fluorescence')
plt.axhspan(-0.3, 1.3, color=sns.xkcd_rgb['medium blue'], alpha=0.2)
ax.legend_.remove()
f.savefig('./figures/downstream/downstreamvar_subgroup4_gfpp_allp1_nsp2f_swarmplot_withdeltag.png', \
          dpi = 300, format='png', bbox_inches='tight', frameon=True)





f=plt.figure(figsize=(9,3))
ax=sns.boxplot(data=df4[(df4.origin!='nsp2f')&
                        [x in list(eventlist[1:])+list(nofslist) for x in df4.fsevent]&(df4.gfpm_wt<25)], 
                         x='fsevent', y='gfpm_wt',hue='origin',\
                          order=eventlist[1:]+nofslist, palette=[sns.xkcd_rgb['light blue'], sns.xkcd_rgb['light green']], 
                                                            fliersize=0)
ax=sns.swarmplot(data=df4[(df4.origin!='nsp2f')&
                        [x in list(eventlist[1:])+list(nofslist) for x in df4.fsevent]&(df4.gfpm_wt<25)], 
                         x='fsevent', y='gfpm_wt',hue='origin',\
                          order=eventlist[1:]+nofslist, palette=[sns.xkcd_rgb['medium blue'], sns.xkcd_rgb['dark green']], 
                                    size=3, split=True)
#wtvalues.loc[eventlist+nofslist,'gfp_nrpeaks1'].plot(marker='o',linewidth=0, color=sns.xkcd_rgb['dark blue'])
ax.set_xticklabels(eventlist[1:]+nofslist, rotation=45,horizontalalignment='right')
plt.xlim(-0.5,len(eventlist[1:]+nofslist)-0.5)
plt.ylim(-0.3,10)
plt.ylabel('% GFP fluorescence')
plt.axhspan(-0.3, 1.3, color=sns.xkcd_rgb['medium blue'], alpha=0.2)
ax.legend_.remove()
f.savefig('./figures/downstream/downstreamvar_subgroup4_gfp_eventlist_nofslist_boxswarmplot_byorigin.png', 
          dpi = 300, format='png', bbox_inches='tight', frameon=True)

    
df4[(df4.origin!='nsp2f')&
                        [x in list(eventlist[1:])+list(nofslist) for x in df4.fsevent]&(df4.gfpm_wt<25)].pivot_table(
    index='fsevent',columns='origin',values='gfpm_wt', aggfunc=np.count_nonzero).to_csv('./n_numbers/Fig3g.csv')


f=plt.figure(figsize=(9,3))
ax=sns.boxplot(data=df4[(df4.origin!='nsp2f')&
                        [x in list(eventlist[1:])+list(nofslist) for x in df4.fsevent]&(df4.gfpm_wt<25)], 
                         x='fsevent', y='percgfpm_wt',hue='origin',\
                          order=eventlist[1:]+nofslist, palette=[sns.xkcd_rgb['light blue'], sns.xkcd_rgb['light green']], 
                                                            fliersize=0)
ax=sns.swarmplot(data=df4[(df4.origin!='nsp2f')&
                        [x in list(eventlist[1:])+list(nofslist) for x in df4.fsevent]&(df4.gfpm_wt<25)], 
                         x='fsevent', y='percgfpm_wt',hue='origin',\
                          order=eventlist[1:]+nofslist, palette=[sns.xkcd_rgb['medium blue'], sns.xkcd_rgb['dark green']], 
                                    size=3, split=True)
#wtvalues.loc[eventlist+nofslist,'gfp_nrpeaks1'].plot(marker='o',linewidth=0, color=sns.xkcd_rgb['dark blue'])
ax.set_xticklabels(eventlist[1:]+nofslist, rotation=45,horizontalalignment='right')
plt.xlim(-0.5,len(eventlist[1:]+nofslist)-0.5)
plt.axhline(y=100, linewidth=2, color='grey')
plt.ylim(0,500)
plt.ylabel('% GFP fluorescence')
ax.legend_.remove()
f.savefig('./figures/downstream/downstreamvar_subgroup4_percgfp_eventlist_nofslist_boxswarmplot_byorigin.png', 
          dpi = 300, format='png', bbox_inches='tight', frameon=True)







invfold=pd.DataFrame(index=eventlist[1:]+nofslist)

invfold['wtval']=wtvalues.loc[eventlist+nofslist,'gfpm']
invfold['hivmean']=invfold.index.map(lambda x: df4[(df4.fsevent==x)&(df4.gfpm_wt<25)&(df4.origin=='hiv')].gfpm_wt.mean())
invfold['hivmedian']=invfold.index.map(lambda x: df4[(df4.fsevent==x)&(df4.gfpm_wt<25)&(df4.origin=='hiv')].gfpm_wt.median())
invfold['srvlmean']=invfold.index.map(lambda x: df4[(df4.fsevent==x)&(df4.gfpm_wt<25)&(df4.origin=='srvl')].gfpm_wt.mean())
invfold['srvlmedian']=invfold.index.map(lambda x: df4[(df4.fsevent==x)&(df4.gfpm_wt<25)&(df4.origin=='srvl')].gfpm_wt.median())
invfold['hivstd']=invfold.index.map(lambda x: df4[(df4.fsevent==x)&(df4.gfpm_wt<25)&(df4.origin=='hiv')].gfpm_wt.std())
invfold['srvlstd']=invfold.index.map(lambda x: df4[(df4.fsevent==x)&(df4.gfpm_wt<25)&(df4.origin=='srvl')].gfpm_wt.std())

f=plt.figure(figsize=(3,3))
plt.scatter(invfold.wtval, invfold.hivmean, color=sns.xkcd_rgb['medium blue'])
plt.scatter(invfold.wtval, invfold.srvlmean, color=sns.xkcd_rgb['green'])
plt.xlim(0,10)
plt.ylim(0.5,3)
plt.axhspan(0, 1.3, color=sns.xkcd_rgb['medium blue'], alpha=0.15)
plt.axvspan(0, 1.3, color=sns.xkcd_rgb['medium blue'], alpha=0.15)
f.savefig('./figures/downstream/downstreamvar_subgroup4_meangfp_vs_wtvalues_eventlist_nofslist.png', hue='dg_downstream',\
          dpi = 300, format='png', bbox_inches='tight', frameon=True)

f=plt.figure(figsize=(3,3))
plt.scatter(invfold.wtval, invfold.hivmedian, color=sns.xkcd_rgb['medium blue'])
plt.scatter(invfold.wtval, invfold.srvlmedian, color=sns.xkcd_rgb['green'])
plt.xlim(0,10)
plt.ylim(0.5,3)
plt.axhspan(0, 1.3, color=sns.xkcd_rgb['medium blue'], alpha=0.15)
plt.axvspan(0, 1.3, color=sns.xkcd_rgb['medium blue'], alpha=0.15)
f.savefig('./figures/downstream/downstreamvar_subgroup4_mediangfp_vs_wtvalues_eventlist_nofslist.png', hue='dg_downstream',\
          dpi = 300, format='png', bbox_inches='tight', frameon=True)

f=plt.figure(figsize=(12,12))
ax=f.add_subplot(111)
plt.scatter(invfold.wtval, invfold.hivmean, color=sns.xkcd_rgb['medium blue'])
plt.scatter(invfold.wtval, invfold.srvlmean, color=sns.xkcd_rgb['green'])
plt.xlim(0,15)
plt.ylim(0.5,3)
plt.axhspan(0, 1.3, color=sns.xkcd_rgb['medium blue'], alpha=0.15)
plt.axvspan(0, 1.3, color=sns.xkcd_rgb['medium blue'], alpha=0.15)
for label, x, y in zip(invfold.index, invfold.wtval, invfold.hivmean):
    ax.annotate(label, xy=(x+0.02,y),fontsize=10, color='b')
for label, x, y in zip(invfold.index, invfold.wtval, invfold.srvlmean):
    ax.annotate(label, xy=(x+0.02,y),fontsize=10, color='g')
f.savefig('./figures/downstream/downstreamvar_subgroup4_meangfp_vs_wtvalues_eventlist_nofslist_withannotation.png', hue='dg_downstream',\
          dpi = 300, format='png', bbox_inches='tight', frameon=True)



df['wtvalue']=df.fsevent.apply(lambda x: wtvalues.loc[x,'gfpm'])

f=plt.figure(figsize=(4,3))
sns.regplot(data=df[(df.gfpm_wt<25)&(df.subgroup==2)], 
                    x='wtvalue', y='gfpm_wt', x_jitter=0.1, 
                    scatter_kws={'alpha':0.2}, fit_reg=False)
plt.xlim(0,10)
plt.ylim(0,15)
f.savefig('./figures/downstream/downstreamvar_subgroup2_regplot_gfp_vs_wtvalues_eventlist_nofslist.png', hue='dg_downstream',\
          dpi = 300, format='png', bbox_inches='tight', frameon=True)


f=plt.figure(figsize=(4,3))
sns.regplot(data=df[(df.gfpm_wt<25)&((df.subgroup==3)|(df.subgroup==4))], 
                    x='wtvalue', y='gfpm_wt', x_jitter=0.1, 
                    scatter_kws={'alpha':0.2}, fit_reg=False)
plt.xlim(0,10)
plt.ylim(0,15)
f.savefig('./figures/downstream/downstreamvar_subgroup34_regplot_gfp_vs_wtvalues_eventlist_nofslist.png', hue='dg_downstream',\
          dpi = 300, format='png', bbox_inches='tight', frameon=True)


# plot correlations with secondary structure together


dgcorrs40=pd.DataFrame()
dgcorrs40p=pd.DataFrame()

for i in eventlist:
    dgcorrs40.loc[i,'point mutations']=pearsonr(df[(df.laststopinframem1<12)&(df.subgroup==2)&(df.fsevent==i)&(df.percgfpm_wt<300)&(df.gfpm_wt<25)][['gfpm_wt','dg40']].dropna().gfpm_wt, \
                   df[(df.laststopinframem1<12)&(df.subgroup==2)&(df.fsevent==i)&(df.gfpm_wt<25)&(df.percgfpm_wt<300)][['gfpm_wt','dg40']].dropna().dg40)[0]
    dgcorrs40p.loc[i,'point mutations']=pearsonr(df[(df.laststopinframem1<12)&(df.subgroup==2)&(df.fsevent==i)&(df.percgfpm_wt<300)&(df.gfpm_wt<25)][['gfpm_wt','dg40']].dropna().gfpm_wt, \
                   df[(df.laststopinframem1<12)&(df.subgroup==2)&(df.fsevent==i)&(df.gfpm_wt<25)&(df.percgfpm_wt<300)][['gfpm_wt','dg40']].dropna().dg40)[1]

for i in eventlist:
    dgcorrs40.loc[i,'structural variants']=pearsonr(df[(df.laststopinframem1<12)&((df.subgroup==3)|(df.subgroup==4))&(df.fsevent==i)&(df.percgfpm_wt<300)&(df.gfpm_wt<25)][['gfpm_wt','dg40']].dropna().gfpm_wt, \
                   df[(df.laststopinframem1<12)&((df.subgroup==3)|(df.subgroup==4))&(df.fsevent==i)&(df.percgfpm_wt<300)&(df.gfpm_wt<25)][['gfpm_wt','dg40']].dropna().dg40)[0]
    dgcorrs40p.loc[i,'structural variants']=pearsonr(df[(df.laststopinframem1<12)&((df.subgroup==3)|(df.subgroup==4))&(df.fsevent==i)&(df.percgfpm_wt<300)&(df.gfpm_wt<25)][['gfpm_wt','dg40']].dropna().gfpm_wt, \
                   df[(df.laststopinframem1<12)&((df.subgroup==3)|(df.subgroup==4))&(df.fsevent==i)&(df.percgfpm_wt<300)&(df.gfpm_wt<25)][['gfpm_wt','dg40']].dropna().dg40)[1]

for i in eventlist:
    dgcorrs40.loc[i,'all variants']=pearsonr(df[(df.laststopinframem1<12)&(df.endoslippery==True)&(df.fsevent==i)&(df.percgfpm_wt<300)&(df.gfpm_wt<25)][['gfpm_wt','dg40']].dropna().gfpm_wt, \
                   fsdf[(df.laststopinframem1<12)&(fsdf.endoslippery==True)&(fsdf.fsevent==i)&(df.percgfpm_wt<300)&(df.gfpm_wt<25)][['gfpm_wt','dg40']].dropna().dg40)[0]
    dgcorrs40p.loc[i,'all variants']=pearsonr(df[(df.laststopinframem1<12)&(df.endoslippery==True)&(df.fsevent==i)&(df.percgfpm_wt<300)&(df.gfpm_wt<25)][['gfpm_wt','dg40']].dropna().gfpm_wt, \
                   df[(df.laststopinframem1<12)&(df.endoslippery==True)&(df.fsevent==i)&(df.percgfpm_wt<300)&(df.gfpm_wt<25)][['gfpm_wt','dg40']].dropna().dg40)[1]

cg=sns.clustermap(data=dgcorrs40.dropna(), figsize=(4,4), 
                  annot=dgcorrs40p.dropna(), 
                  cbar_kws={'ticks':[-0.4,0,0.4]},
               annot_kws={'fontsize':12}, fmt='.2g', col_cluster=False)

cg=sns.clustermap(data=dgcorrs40.dropna(), figsize=(4,4), 
                  annot=dgcorrs40p.dropna().iloc[cg.dendrogram_row.reordered_ind], 
                  cbar_kws={'ticks':[-0.4,0,0.4]},
               annot_kws={'fontsize':12}, fmt='.2g', col_cluster=False)
plt.setp(cg.ax_heatmap.yaxis.get_majorticklabels(), rotation=0)
cg.savefig('./figures/clustermap_structure_dg40_correlations_bysubgroup_pearson.png',
              dpi = 300, format='png', bbox_inches='tight', frameon=True)



# plot correlations with secondary structure together - dg_downstream


dgcorrs40=pd.DataFrame()
dgcorrs40p=pd.DataFrame()

for i in eventlist:
    dgcorrs40.loc[i,'point mutations']=pearsonr(df[(df.laststopinframem1<12)&(df.subgroup==2)&(df.fsevent==i)&(df.percgfpm_wt<300)&(df.gfpm_wt<25)][['gfpm_wt','dg_downstream']].dropna().gfpm_wt, \
                   df[(df.laststopinframem1<12)&(df.subgroup==2)&(df.fsevent==i)&(df.gfpm_wt<25)&(df.percgfpm_wt<300)][['gfpm_wt','dg_downstream']].dropna().dg_downstream)[0]
    dgcorrs40p.loc[i,'point mutations']=pearsonr(df[(df.laststopinframem1<12)&(df.subgroup==2)&(df.fsevent==i)&(df.percgfpm_wt<300)&(df.gfpm_wt<25)][['gfpm_wt','dg_downstream']].dropna().gfpm_wt, \
                   df[(df.laststopinframem1<12)&(df.subgroup==2)&(df.fsevent==i)&(df.gfpm_wt<25)&(df.percgfpm_wt<300)][['gfpm_wt','dg_downstream']].dropna().dg_downstream)[1]

for i in eventlist:
    dgcorrs40.loc[i,'structural variants']=pearsonr(df[(df.laststopinframem1<12)&((df.subgroup==3)|(df.subgroup==4))&(df.fsevent==i)&(df.percgfpm_wt<300)&(df.gfpm_wt<25)][['gfpm_wt','dg_downstream']].dropna().gfpm_wt, \
                   df[(df.laststopinframem1<12)&((df.subgroup==3)|(df.subgroup==4))&(df.fsevent==i)&(df.percgfpm_wt<300)&(df.gfpm_wt<25)][['gfpm_wt','dg_downstream']].dropna().dg_downstream)[0]
    dgcorrs40p.loc[i,'structural variants']=pearsonr(df[(df.laststopinframem1<12)&((df.subgroup==3)|(df.subgroup==4))&(df.fsevent==i)&(df.percgfpm_wt<300)&(df.gfpm_wt<25)][['gfpm_wt','dg_downstream']].dropna().gfpm_wt, \
                   df[(df.laststopinframem1<12)&((df.subgroup==3)|(df.subgroup==4))&(df.fsevent==i)&(df.percgfpm_wt<300)&(df.gfpm_wt<25)][['gfpm_wt','dg_downstream']].dropna().dg_downstream)[1]

for i in eventlist:
    dgcorrs40.loc[i,'all variants']=pearsonr(df[(df.laststopinframem1<12)&(df.endoslippery==True)&(df.fsevent==i)&(df.percgfpm_wt<300)&(df.gfpm_wt<25)][['gfpm_wt','dg_downstream']].dropna().gfpm_wt, \
                   fsdf[(df.laststopinframem1<12)&(fsdf.endoslippery==True)&(fsdf.fsevent==i)&(df.percgfpm_wt<300)&(df.gfpm_wt<25)][['gfpm_wt','dg_downstream']].dropna().dg_downstream)[0]
    dgcorrs40p.loc[i,'all variants']=pearsonr(df[(df.laststopinframem1<12)&(df.endoslippery==True)&(df.fsevent==i)&(df.percgfpm_wt<300)&(df.gfpm_wt<25)][['gfpm_wt','dg_downstream']].dropna().gfpm_wt, \
                   df[(df.laststopinframem1<12)&(df.endoslippery==True)&(df.fsevent==i)&(df.percgfpm_wt<300)&(df.gfpm_wt<25)][['gfpm_wt','dg_downstream']].dropna().dg_downstream)[1]

cg=sns.clustermap(data=dgcorrs40.dropna(), figsize=(4,4), 
                  annot=dgcorrs40p.dropna(), 
                  cbar_kws={'ticks':[-0.4,0,0.4]},
               annot_kws={'fontsize':12}, fmt='.2g', col_cluster=False)

cg=sns.clustermap(data=dgcorrs40.dropna(), figsize=(4,4), 
                  annot=dgcorrs40p.dropna().iloc[cg.dendrogram_row.reordered_ind], 
                  cbar_kws={'ticks':[-0.4,0,0.4]},
               annot_kws={'fontsize':12}, fmt='.2g', col_cluster=False)
plt.setp(cg.ax_heatmap.yaxis.get_majorticklabels(), rotation=0)
cg.savefig('./figures/clustermap_structure_dgdownstream_correlations_bysubgroup_pearson.png',
              dpi = 300, format='png', bbox_inches='tight', frameon=True)



# plot correlations with secondary structure together - mfe [pKiss]


dgcorrs40=pd.DataFrame()
dgcorrs40p=pd.DataFrame()

for i in eventlist:
    dgcorrs40.loc[i,'point mutations']=pearsonr(df[(df.laststopinframem1<12)&(df.subgroup==2)&(df.fsevent==i)&(df.percgfpm_wt<300)&(df.gfpm_wt<25)][['gfpm_wt','mfe']].dropna().gfpm_wt, \
                   df[(df.laststopinframem1<12)&(df.subgroup==2)&(df.fsevent==i)&(df.gfpm_wt<25)&(df.percgfpm_wt<300)][['gfpm_wt','mfe']].dropna().mfe)[0]
    dgcorrs40p.loc[i,'point mutations']=pearsonr(df[(df.laststopinframem1<12)&(df.subgroup==2)&(df.fsevent==i)&(df.percgfpm_wt<300)&(df.gfpm_wt<25)][['gfpm_wt','mfe']].dropna().gfpm_wt, \
                   df[(df.laststopinframem1<12)&(df.subgroup==2)&(df.fsevent==i)&(df.gfpm_wt<25)&(df.percgfpm_wt<300)][['gfpm_wt','mfe']].dropna().mfe)[1]

for i in eventlist:
    dgcorrs40.loc[i,'structural variants']=pearsonr(df[(df.laststopinframem1<12)&((df.subgroup==3)|(df.subgroup==4))&(df.fsevent==i)&(df.percgfpm_wt<300)&(df.gfpm_wt<25)][['gfpm_wt','mfe']].dropna().gfpm_wt, \
                   df[(df.laststopinframem1<12)&((df.subgroup==3)|(df.subgroup==4))&(df.fsevent==i)&(df.percgfpm_wt<300)&(df.gfpm_wt<25)][['gfpm_wt','mfe']].dropna().mfe)[0]
    dgcorrs40p.loc[i,'structural variants']=pearsonr(df[(df.laststopinframem1<12)&((df.subgroup==3)|(df.subgroup==4))&(df.fsevent==i)&(df.percgfpm_wt<300)&(df.gfpm_wt<25)][['gfpm_wt','mfe']].dropna().gfpm_wt, \
                   df[(df.laststopinframem1<12)&((df.subgroup==3)|(df.subgroup==4))&(df.fsevent==i)&(df.percgfpm_wt<300)&(df.gfpm_wt<25)][['gfpm_wt','mfe']].dropna().mfe)[1]

for i in eventlist:
    dgcorrs40.loc[i,'all variants']=pearsonr(df[(df.laststopinframem1<12)&(df.endoslippery==True)&(df.fsevent==i)&(df.percgfpm_wt<300)&(df.gfpm_wt<25)][['gfpm_wt','mfe']].dropna().gfpm_wt, \
                   df[(df.laststopinframem1<12)&(df.endoslippery==True)&(df.fsevent==i)&(df.percgfpm_wt<300)&(df.gfpm_wt<25)][['gfpm_wt','mfe']].dropna().mfe)[0]
    dgcorrs40p.loc[i,'all variants']=pearsonr(df[(df.laststopinframem1<12)&(df.endoslippery==True)&(df.fsevent==i)&(df.percgfpm_wt<300)&(df.gfpm_wt<25)][['gfpm_wt','mfe']].dropna().gfpm_wt, \
                   df[(df.laststopinframem1<12)&(df.endoslippery==True)&(df.fsevent==i)&(df.percgfpm_wt<300)&(df.gfpm_wt<25)][['gfpm_wt','mfe']].dropna().mfe)[1]

cg=sns.clustermap(data=dgcorrs40.dropna(), figsize=(4,4), 
                  annot=dgcorrs40p.dropna(), 
                  cbar_kws={'ticks':[-0.4,0,0.4]},
               annot_kws={'fontsize':12}, fmt='.2g', col_cluster=False)

cg=sns.clustermap(data=dgcorrs40.dropna(), figsize=(4,4), 
                  annot=dgcorrs40p.dropna().iloc[cg.dendrogram_row.reordered_ind], 
                  cbar_kws={'ticks':[-0.4,0,0.4]},
               annot_kws={'fontsize':12}, fmt='.2g', col_cluster=False)
plt.setp(cg.ax_heatmap.yaxis.get_majorticklabels(), rotation=0)
cg.savefig('./figures/clustermap_structure_mfe_correlations_bysubgroup_pearson.png',
              dpi = 300, format='png', bbox_inches='tight', frameon=True)


# Plot correlations for individual events

for event in eventlist:
    f=plt.figure(figsize=(4,3))
    ax=plt.scatter(df[(df.subgroup==2)&(df.fsevent==event)&(df.gfpm_wt<25)&(df.percgfpm_wt<300)][['gfpm_wt','dg40','rnaperdna']].dropna().dg40, \
              df[(df.subgroup==2)&(df.fsevent==event)&(df.gfpm_wt<25)&(df.percgfpm_wt<300)][['gfpm_wt','dg40','rnaperdna']].dropna().gfpm_wt,
                 c=df[(df.subgroup==2)&(df.fsevent==event)&(df.gfpm_wt<25)&(df.percgfpm_wt<300)][['gfpm_wt','dg40','rnaperdna']].dropna().rnaperdna, cmap='Blues')
    plt.axhline(y=wtvalues.loc[event,'gfpm'], linewidth=2, color='grey')
    plt.axvline(x=fsdf[(fsdf.subset=='fsmain')&(fsdf.fsevent==event)].dg40.unique(), linewidth=2, color='grey')
    plt.colorbar(ax, label='log2(RNA/DNA reads)')
    plt.xlabel('MFE (Vienna RNA)')
    plt.ylabel('% GFP fluorescence')
    plt.ylim(0)
    plt.title(event+'\nr={:.2f}'.format(pearsonr(df[(df.subgroup==2)&(df.fsevent==event)&(df.gfpm_wt<25)&(df.percgfpm_wt<300)][['gfpm_wt','dg40']].dropna().gfpm_wt, \
                   df[(df.subgroup==2)&(df.fsevent==event)&(df.gfpm_wt<25)&(df.percgfpm_wt<300)][['gfpm_wt','dg40']].dropna().dg40)[0]) + \
    ', p={:.1e}'.format(pearsonr(df[(df.subgroup==2)&(df.fsevent==event)&(df.gfpm_wt<25)&(df.percgfpm_wt<300)][['gfpm_wt','dg40']].dropna().gfpm_wt, \
                   df[(df.subgroup==2)&(df.fsevent==event)&(df.gfpm_wt<25)&(df.percgfpm_wt<300)][['gfpm_wt','dg40']].dropna().dg40)[1]), \
                   fontsize=14)
    f.savefig('./figures/downstream/downstreamvar_subgroup2_correlationsgfpwithdeltag40_'+'_'.join(event) +'.png', \
              dpi = 300, format='png', bbox_inches='tight', frameon=True)

for event in eventlist:
    f=plt.figure(figsize=(4,3))
    ax=plt.scatter(df[((df.subgroup==3)|(df.subgroup==4))&(df.fsevent==event)&(df.gfpm_wt<25)&(df.percgfpm_wt<300)][['gfpm_wt','dg40','rnaperdna']].dropna().dg40, \
              df[((df.subgroup==3)|(df.subgroup==4))&(df.fsevent==event)&(df.gfpm_wt<25)&(df.percgfpm_wt<300)][['gfpm_wt','dg40','rnaperdna']].dropna().gfpm_wt, 
                 c=df[((df.subgroup==3)|(df.subgroup==4))&(df.fsevent==event)&(df.gfpm_wt<25)&(df.percgfpm_wt<300)][['gfpm_wt','dg40','rnaperdna']].dropna().rnaperdna, cmap='Blues')
    plt.axhline(y=wtvalues.loc[event,'gfpm'], linewidth=2, color='grey')
    plt.axvline(x=fsdf[(fsdf.subset=='fsmain')&(fsdf.fsevent==event)].dg40.unique(), linewidth=2, color='grey')
    plt.colorbar(ax, label='log2(RNA/DNA reads)')
    plt.xlabel('MFE (Vienna RNA)')
    plt.ylabel('% GFP fluorescence')
    plt.ylim(0)
    plt.title(event+'\nr={:.2f}'.format(pearsonr(df[((df.subgroup==3)|(df.subgroup==4))&(df.fsevent==event)&(df.gfpm_wt<25)&(df.percgfpm_wt<300)][['gfpm_wt','dg40']].dropna().gfpm_wt, \
                   df[((df.subgroup==3)|(df.subgroup==4))&(df.fsevent==event)&(df.gfpm_wt<25)&(df.percgfpm_wt<300)][['gfpm_wt','dg40']].dropna().dg40)[0]) + \
    ', p={:.1e}'.format(pearsonr(df[((df.subgroup==3)|(df.subgroup==4))&(df.fsevent==event)&(df.gfpm_wt<25)&(df.percgfpm_wt<300)][['gfpm_wt','dg40']].dropna().gfpm_wt, \
                   df[((df.subgroup==3)|(df.subgroup==4))&(df.fsevent==event)&(df.gfpm_wt<25)&(df.percgfpm_wt<300)][['gfpm_wt','dg40']].dropna().dg40)[1]), \
                   fontsize=14)
    f.savefig('./figures/downstream/downstreamvar_subgroup34_correlationsgfpwithdeltag40_'+'_'.join(event) +'.png', \
              dpi = 300, format='png', bbox_inches='tight', frameon=True)

for event in eventlist:
    f=plt.figure(figsize=(4,3))
    ax=plt.scatter(df[(df.endoslippery==True)&(df.fsevent==event)&(df.gfpm_wt<25)&(df.percgfpm_wt<300)][['gfpm_wt','dg40','rnaperdna']].dropna().dg40, \
              df[(df.endoslippery==True)&(df.fsevent==event)&(df.gfpm_wt<25)&(df.percgfpm_wt<300)][['gfpm_wt','dg40','rnaperdna']].dropna().gfpm_wt,
                 c=df[(df.endoslippery==True)&(df.fsevent==event)&(df.gfpm_wt<25)&(df.percgfpm_wt<300)][['gfpm_wt','dg40','rnaperdna']].dropna().rnaperdna, cmap='Blues')
    plt.axhline(y=wtvalues.loc[event,'gfpm'], linewidth=2, color='grey')
    plt.axvline(x=fsdf[(fsdf.subset=='fsmain')&(fsdf.fsevent==event)].dg40.unique(), linewidth=2, color='grey')
    plt.colorbar(ax, label='log2(RNA/DNA reads)')
    plt.xlabel('MFE (Vienna RNA)')
    plt.ylabel('% GFP fluorescence')
    plt.ylim(0)
    plt.title(event+'\nr={:.2f}'.format(pearsonr(df[(df.endoslippery==True)&(df.fsevent==event)&(df.gfpm_wt<25)&(df.percgfpm_wt<300)][['gfpm_wt','dg40','rnaperdna']].dropna().gfpm_wt, \
                   df[(df.endoslippery==True)&(df.fsevent==event)&(df.gfpm_wt<25)&(df.percgfpm_wt<300)][['gfpm_wt','dg40','rnaperdna']].dropna().dg40)[0]) + \
    ', p={:.1e}'.format(pearsonr(df[(df.endoslippery==True)&(df.fsevent==event)&(df.gfpm_wt<25)&(df.percgfpm_wt<300)][['gfpm_wt','dg40','rnaperdna']].dropna().gfpm_wt, \
                   df[(df.endoslippery==True)&(df.fsevent==event)&(df.gfpm_wt<25)&(df.percgfpm_wt<300)][['gfpm_wt','dg40','rnaperdna']].dropna().dg40)[1]), \
                   fontsize=14)
    f.savefig('./figures/downstream/downstreamvar_allvariants_correlationsgfpwithdeltag40_'+'_'.join(event) +'.png', \
              dpi = 300, format='png', bbox_inches='tight', frameon=True)
    
    
for event in eventlist:
    f=plt.figure(figsize=(3,3))
    ax=plt.scatter(df[(df.endoslippery==True)&(df.fsevent==event)&(df.gfpm_wt<25)&(df.percgfpm_wt<300)][['dg40','rnaperdna']].dropna().dg40, \
              df[(df.endoslippery==True)&(df.fsevent==event)&(df.gfpm_wt<25)&(df.percgfpm_wt<300)][['dg40','rnaperdna']].dropna().rnaperdna)
    plt.axvline(x=fsdf[(fsdf.subset=='fsmain')&(fsdf.fsevent==event)].dg40.unique(), linewidth=2, color='grey')
    plt.axhline(y=fsdf[(fsdf.subset=='fsmain')&(fsdf.fsevent==event)&(fsdf.gfpm_wt<25)].rnaperdna.dropna().median(), linewidth=2, color='grey')
    plt.xlabel('MFE (Vienna RNA)')
    plt.ylabel('log2(RNA/DNA reads)')
    plt.title(event+'\nr={:.2f}'.format(pearsonr(df[(df.endoslippery==True)&(df.fsevent==event)&(df.gfpm_wt<25)&(df.percgfpm_wt<300)][['dg40','rnaperdna']].replace(to_replace=-np.inf, value=np.nan).dropna().rnaperdna, \
                   df[(df.endoslippery==True)&(df.fsevent==event)&(df.gfpm_wt<25)&(df.percgfpm_wt<300)][['dg40','rnaperdna']].replace(to_replace=-np.inf, value=np.nan).dropna().dg40)[0]) + \
    ', p={:.1e}'.format(pearsonr(df[(df.endoslippery==True)&(df.fsevent==event)&(df.gfpm_wt<25)&(df.percgfpm_wt<300)][['dg40','rnaperdna']].replace(to_replace=-np.inf, value=np.nan).dropna().rnaperdna, \
                   df[(df.endoslippery==True)&(df.fsevent==event)&(df.gfpm_wt<25)&(df.percgfpm_wt<300)][['dg40','rnaperdna']].replace(to_replace=-np.inf, value=np.nan).dropna().dg40)[1]), \
                   fontsize=14)
    f.savefig('./figures/downstream/downstreamvar_allvariants_correlationsrnaperdnawithdeltag40_'+'_'.join(event) +'.png', \
              dpi = 300, format='png', bbox_inches='tight', frameon=True)
    

# with mfe from pKiss

for event in eventlist:
    f=plt.figure(figsize=(4,3))
    ax=plt.scatter(df[(df.subgroup==2)&(df.fsevent==event)&(df.gfpm_wt<25)&(df.percgfpm_wt<300)][['gfpm_wt','mfe','rnaperdna']].dropna().mfe, \
              df[(df.subgroup==2)&(df.fsevent==event)&(df.gfpm_wt<25)&(df.percgfpm_wt<300)][['gfpm_wt','mfe','rnaperdna']].dropna().gfpm_wt,
                 c=df[(df.subgroup==2)&(df.fsevent==event)&(df.gfpm_wt<25)&(df.percgfpm_wt<300)][['gfpm_wt','mfe','rnaperdna']].dropna().rnaperdna, cmap='Blues')
    plt.axvline(x=df[(df.varseq162==wtseqs[event])].mfe.dropna().unique(), linewidth=2, color='grey')
    plt.axhline(y=wtvalues.loc[event,'gfpm'], linewidth=2, color='grey')
    plt.colorbar(ax, label='log2(RNA/DNA reads)')
    plt.xlabel('MFE (pKiss)')
    plt.ylabel('% GFP fluorescence')
    plt.ylim(0)
    plt.title(event+'\nr={:.2f}'.format(pearsonr(df[(df.subgroup==2)&(df.fsevent==event)&(df.gfpm_wt<25)&(df.percgfpm_wt<300)][['gfpm_wt','mfe','rnaperdna']].dropna().gfpm_wt, \
                   df[(df.subgroup==2)&(df.fsevent==event)&(df.gfpm_wt<25)&(df.percgfpm_wt<300)][['gfpm_wt','mfe','rnaperdna']].dropna().mfe)[0]) + \
    ', p={:.1e}'.format(pearsonr(df[(df.subgroup==2)&(df.fsevent==event)&(df.gfpm_wt<25)&(df.percgfpm_wt<300)][['gfpm_wt','mfe','rnaperdna']].dropna().gfpm_wt, \
                   df[(df.subgroup==2)&(df.fsevent==event)&(df.gfpm_wt<25)&(df.percgfpm_wt<300)][['gfpm_wt','mfe','rnaperdna']].dropna().mfe)[1]), \
                   fontsize=14)
    f.savefig('./figures/downstream/downstreamvar_subgroup2_correlationsgfpwithMFEpKiss_'+'_'.join(event) +'.png', \
              dpi = 300, format='png', bbox_inches='tight', frameon=True)
    
    
    

for event in eventlist:
    f=plt.figure(figsize=(4,3))
    ax=plt.scatter(df[((df.subgroup==3)|(df.subgroup==4))&(df.fsevent==event)&(df.gfpm_wt<25)&(df.percgfpm_wt<300)][['gfpm_wt','mfe','rnaperdna']].dropna().mfe, \
              df[((df.subgroup==3)|(df.subgroup==4))&(df.fsevent==event)&(df.gfpm_wt<25)&(df.percgfpm_wt<300)][['gfpm_wt','mfe','rnaperdna']].dropna().gfpm_wt,
                 c=df[((df.subgroup==3)|(df.subgroup==4))&(df.fsevent==event)&(df.gfpm_wt<25)&(df.percgfpm_wt<300)][['gfpm_wt','mfe','rnaperdna']].dropna().rnaperdna, cmap='Blues')
    plt.axvline(x=df[(df.varseq162==wtseqs[event])].mfe.dropna().unique(), linewidth=2, color='grey')
    plt.axhline(y=wtvalues.loc[event,'gfpm'], linewidth=2, color='grey')
    plt.colorbar(ax, label='log2(RNA/DNA reads)')
    plt.xlabel('MFE (pKiss)')
    plt.ylabel('% GFP fluorescence')
    plt.ylim(0)
    plt.title(event+'\nr={:.2f}'.format(pearsonr(df[((df.subgroup==3)|(df.subgroup==4))&(df.fsevent==event)&(df.gfpm_wt<25)&(df.percgfpm_wt<300)][['gfpm_wt','mfe','rnaperdna']].dropna().gfpm_wt, \
                   df[((df.subgroup==3)|(df.subgroup==4))&(df.fsevent==event)&(df.gfpm_wt<25)&(df.percgfpm_wt<300)][['gfpm_wt','mfe','rnaperdna']].dropna().mfe)[0]) + \
    ', p={:.1e}'.format(pearsonr(df[((df.subgroup==3)|(df.subgroup==4))&(df.fsevent==event)&(df.gfpm_wt<25)&(df.percgfpm_wt<300)][['gfpm_wt','mfe','rnaperdna']].dropna().gfpm_wt, \
                   df[((df.subgroup==3)|(df.subgroup==4))&(df.fsevent==event)&(df.gfpm_wt<25)&(df.percgfpm_wt<300)][['gfpm_wt','mfe','rnaperdna']].dropna().mfe)[1]), \
                   fontsize=14)
    f.savefig('./figures/downstream/downstreamvar_subgroup34_correlationsgfpwithMFEpKiss_'+'_'.join(event) +'.png', \
              dpi = 300, format='png', bbox_inches='tight', frameon=True)
    
    
    
    
for event in eventlist:
    f=plt.figure(figsize=(4,3))
    ax=plt.scatter(df[(df.endoslippery==True)&(df.fsevent==event)&(df.gfpm_wt<25)&(df.percgfpm_wt<300)][['gfpm_wt','mfe','rnaperdna']].dropna().mfe, \
              df[(df.endoslippery==True)&(df.fsevent==event)&(df.gfpm_wt<25)&(df.percgfpm_wt<300)][['gfpm_wt','mfe','rnaperdna']].dropna().gfpm_wt,
                 c=df[(df.endoslippery==True)&(df.fsevent==event)&(df.gfpm_wt<25)&(df.percgfpm_wt<300)][['gfpm_wt','mfe','rnaperdna']].dropna().rnaperdna, cmap='Blues')
    plt.axvline(x=df[(df.varseq162==wtseqs[event])].mfe.dropna().unique(), linewidth=2, color='grey')
    plt.axhline(y=wtvalues.loc[event,'gfpm'], linewidth=2, color='grey')
    plt.colorbar(ax, label='log2(RNA/DNA reads)')
    plt.xlabel('MFE (pKiss)')
    plt.ylabel('% GFP fluorescence')
    plt.ylim(0)
    plt.title(event+'\nr={:.2f}'.format(pearsonr(df[(df.endoslippery==True)&(df.fsevent==event)&(df.gfpm_wt<25)&(df.percgfpm_wt<300)][['gfpm_wt','mfe','rnaperdna']].dropna().gfpm_wt, \
                   df[(df.endoslippery==True)&(df.fsevent==event)&(df.gfpm_wt<25)&(df.percgfpm_wt<300)][['gfpm_wt','mfe','rnaperdna']].dropna().mfe)[0]) + \
    ', p={:.1e}'.format(pearsonr(df[(df.endoslippery==True)&(df.fsevent==event)&(df.gfpm_wt<25)&(df.percgfpm_wt<300)][['gfpm_wt','mfe','rnaperdna']].dropna().gfpm_wt, \
                   df[(df.endoslippery==True)&(df.fsevent==event)&(df.gfpm_wt<25)&(df.percgfpm_wt<300)][['gfpm_wt','mfe','rnaperdna']].dropna().mfe)[1]), \
                   fontsize=14)
    f.savefig('./figures/downstream/downstreamvar_allvariants_correlationsgfpwithMFEpKiss_'+'_'.join(event) +'.png', \
              dpi = 300, format='png', bbox_inches='tight', frameon=True)
    
    
for event in eventlist:
    f=plt.figure(figsize=(3,3))
    ax=plt.scatter(df[(df.endoslippery==True)&(df.fsevent==event)&(df.gfpm_wt<25)&(df.percgfpm_wt<300)][['mfe','rnaperdna']].dropna().mfe, \
              df[(df.endoslippery==True)&(df.fsevent==event)&(df.gfpm_wt<25)&(df.percgfpm_wt<300)][['mfe','rnaperdna']].dropna().rnaperdna)
    plt.axvline(x=df[(df.varseq162==wtseqs[event])].mfe.dropna().unique(), linewidth=2, color='grey')
    plt.axhline(y=fsdf[(fsdf.subset=='fsmain')&(fsdf.fsevent==event)&(fsdf.gfpm_wt<25)].rnaperdna.dropna().median(), linewidth=2, color='grey')
    plt.xlabel('MFE (pKiss)')
    plt.ylabel('log2(RNA/DNA reads)')
    plt.title(event+'\nr={:.2f}'.format(pearsonr(df[(df.endoslippery==True)&(df.fsevent==event)&(df.gfpm_wt<25)&(df.percgfpm_wt<300)][['mfe','rnaperdna']].replace(to_replace=-np.inf, value=np.nan).dropna().rnaperdna, \
                   df[(df.endoslippery==True)&(df.fsevent==event)&(df.gfpm_wt<25)&(df.percgfpm_wt<300)][['mfe','rnaperdna']].replace(to_replace=-np.inf, value=np.nan).dropna().mfe)[0]) + \
    ', p={:.1e}'.format(pearsonr(df[(df.endoslippery==True)&(df.fsevent==event)&(df.gfpm_wt<25)&(df.percgfpm_wt<300)][['mfe','rnaperdna']].replace(to_replace=-np.inf, value=np.nan).dropna().rnaperdna, \
                   df[(df.endoslippery==True)&(df.fsevent==event)&(df.gfpm_wt<25)&(df.percgfpm_wt<300)][['mfe','rnaperdna']].replace(to_replace=-np.inf, value=np.nan).dropna().mfe)[1]), \
                   fontsize=14)
    f.savefig('./figures/downstream/downstreamvar_allvariants_correlationsrnaperdnawithMFEpKiss_'+'_'.join(event) +'.png', \
              dpi = 300, format='png', bbox_inches='tight', frameon=True)
    

    
#%%
df=fsdf[(fsdf.subset=="fs_synthvar")]
df=df.dropna(axis=1, how='all')
df['subgroup']=df.changes.apply(lambda x: x.split(' ')[0])

df['upstream']=df.varseq162.apply(lambda x: x[:30])
df['slippery']=df.varseq162.apply(lambda x: x[30:42])
df['downstream']=df.varseq162.apply(lambda x: x[42:])
  

f=plt.figure(figsize=(10,3))
#ax=sns.boxplot(data=df[(df.gfp<10)&(df.smoothedpeaks==1)], x='fsevent', y='gfp', color=sns.xkcd_rgb['light blue'], whis=1)
ax=sns.swarmplot(data=df[(df.gfpm_wt<10)&(df.numberreadsm_wt>=10)&(df.peaksm_wt==1)], x='fsevent', y='gfpm_wt', hue='dg_downstream', palette='Blues_d')
ax.set_xticklabels(df[(df.gfpm_wt<10)&(df.numberreadsm_wt>=10)&(df.peaksm_wt==1)].fsevent.unique(), rotation=45,horizontalalignment='right')
plt.ylim(-0.2,10.2)
plt.axhspan(-0.2,1.7, color=sns.xkcd_rgb['medium blue'], alpha=0.2)
plt.ylabel('% GFP fluorescence')
ax.legend_.remove()
plt.xlabel('upstream and slippery site')
f.savefig('./figures/synthetic/synthvar_swarmplot_gfp_by_fsevent_max10_withdg.png', \
          dpi = 300, format='png', bbox_inches='tight', frameon=True)

f=plt.figure(figsize=(10,3))
#ax=sns.boxplot(data=df[(df.gfp<10)&(df.smoothedpeaks==1)], x='fsevent', y='gfp', color=sns.xkcd_rgb['light blue'], whis=1)
ax=sns.swarmplot(data=df[(df.gfpp_wt<25)&(df.numberreadsp_wt>=20)&(df.peaksp_wt==1)], x='fsevent', y='gfpp_wt', hue='dg_downstream', palette='Blues_d')
ax.set_xticklabels(df[(df.gfpp_wt<25)&(df.numberreadsp_wt>=20)&(df.peaksp_wt==1)].fsevent.unique(), rotation=45,horizontalalignment='right')
plt.ylim(-0.2,25.2)
plt.axhspan(-0.2,1.7, color=sns.xkcd_rgb['medium blue'], alpha=0.2)
plt.ylabel('% GFP fluorescence')
ax.legend_.remove()
plt.xlabel('upstream and slippery site')
f.savefig('./figures/synthetic/synthvar_swarmplot_gfp_by_fsevent_max10_withdg.png', \
          dpi = 300, format='png', bbox_inches='tight', frameon=True)




f=sns.clustermap(df[(df.gfpm_wt<25)&(df.numberreadsm_wt>=20)&(df.subgroup=='completely')].pivot_table(index='downstream',
                    columns='slippery',values='gfpm_wt').dropna(), 
            vmin=1.3, figsize=(3,6), cbar_kws={'ticks':[2,5,8,11]})
ax=f.ax_heatmap
ax.set_yticklabels('')
f.savefig('./figures/synthetic/synthvar_clustermap_completelysynth_slipperyvsdownstream_abovebackground_small.png', \
          dpi = 300, format='png', bbox_inches='tight', frameon=True)


f=sns.clustermap(df[(df.gfpp_wt<25)&(df.numberreadsp_wt>=20)&(df.subgroup=='completely')].pivot_table(index='downstream',
                    columns='slippery',values='gfpp_wt').dropna(), 
            vmin=1.3, figsize=(3,6), cbar_kws={'ticks':[2,5,8,11]})
ax=f.ax_heatmap
ax.set_yticklabels('')
f.savefig('./figures/synthetic/synthvar_clustermap_completelysynth_gfpp1_slipperyvsdownstream_abovebackground_small.png', \
          dpi = 300, format='png', bbox_inches='tight', frameon=True)


    
#%%
df=fsdf[(fsdf.subset=='fs_frameshifts')]

df['typechanges']=df.changes.apply(lambda x: x.split(' ')[0])
dfframeshifts=df[(df.gfpm_wt<25)&(df.peaksm_wt==1)&(df.fsevent!='')&(df.numberreadsm_wt>=20)&(df.type==-1)].pivot_table(index='fsevent', columns='changes', values='gfpm_wt', aggfunc=np.median)
stopinwt=df[(df.gfpm_wt<25)&(df.peaksm_wt==1)&(df.fsevent!='')&(df.numberreadsm_wt>=20)&(df.typechanges=='endogenous')&(df.laststopinframem1>11)].fsevent.unique()

dfframeshifts.drop('frameshift mimicking', axis=1, inplace=True)

dfframeshifts.to_csv('./tables/frameshifts_table.csv')

ax=dfframeshifts.dropna().T.plot(kind='bar', figsize=(3,3))
f=ax.get_figure()
ax.legend_.remove()
plt.axhspan(0,1.3,alpha=0.2)
f.savefig('./figures/frameshifts/hervk10_all.png', \
          dpi = 300, format='png', bbox_inches='tight', frameon=True)

df[(df.gfpm_wt<25)&(df.peaksm_wt==1)&(df.fsevent!='')&(df.numberreadsm_wt>=20)&(df.type==-1)].pivot_table(index='fsevent', columns='changes', values='gfpm_wt', aggfunc=np.count_nonzero).dropna().to_csv('./n_numbers/FigS4b.csv')

f=plt.figure(figsize=(8,8))               
plt.scatter(dfframeshifts['endogenous sequence'],
            dfframeshifts['stop codons mutated'])
for i in dfframeshifts.index:
    plt.annotate(str(i), xy=(dfframeshifts.loc[i, 'endogenous sequence']+0.03,
            dfframeshifts.loc[i, 'stop codons mutated']+0.03), fontsize=10)
plt.scatter(dfframeshifts.loc[stopinwt,'endogenous sequence'],
            dfframeshifts.loc[stopinwt,'stop codons mutated'], color='r')
plt.plot([0,8],[0,8], '--', color='grey')
plt.xlim(0,12)
plt.ylim(0,8)
f.savefig('./figures/frameshifts/endogneous_vs_stop_minus1_withannotation.png', \
          dpi = 300, format='png', bbox_inches='tight', frameon=True)



f=plt.figure(figsize=(3,3))               
plt.scatter(dfframeshifts['endogenous sequence'],
            dfframeshifts['stop codons mutated'])
plt.scatter(dfframeshifts.loc[stopinwt,'endogenous sequence'],
            dfframeshifts.loc[stopinwt,'stop codons mutated'], color='r')
plt.plot([0,8],[0,8], '--', color='grey')
plt.xlim(0,12)
plt.ylim(0,8)
f.savefig('./figures/frameshifts/endogneous_vs_stop_minus1.png', \
          dpi = 300, format='png', bbox_inches='tight', frameon=True)





for ch in df.changes.unique():
    try:
        with sns.axes_style('dark'):
            ax=df[(df.changes==ch)&(df.gfpm_wt<25)&(df.peaksm_wt==1)&([x not in eventlist+nofslist+plist+[''] for x in df.fsevent])&(df.numberreadsm_wt>=10)&(df.laststopinframem1<12)].pivot_table(index='fsevent',values='gfpm_wt', 
                  aggfunc=np.median).sort_values('gfpm_wt', ascending=False).plot(kind='barh', 
                figsize=(4,len(df[(df.changes==ch)&(df.gfpm_wt<25)&([x not in eventlist+nofslist+plist+[''] for x in df.fsevent])&(df.peaksm_wt==1)&(df.numberreadsm_wt>=10)&(df.laststopinframem1<12)].pivot_table(index='fsevent',
                                  values='gfpm_wt'))*0.3))
            sns.swarmplot(data=df[(df.changes==ch)&(df.gfpm_wt<25)&(df.peaksm_wt==1)&
                    ([x not in eventlist+nofslist+plist+[''] for x in df.fsevent])&
                    (df.numberreadsm_wt>=10)&(df.laststopinframem1<12)], y='fsevent',x='gfpm_wt',
                     order=df[(df.changes==ch)&(df.gfpm_wt<25)&(df.peaksm_wt==1)&
                    ([x not in eventlist+nofslist+plist+[''] for x in df.fsevent])&
                    (df.numberreadsm_wt>=10)&(df.laststopinframem1<12)].pivot_table(
                    index='fsevent',values='gfpm_wt', 
                  aggfunc=np.median).sort_values('gfpm_wt', ascending=False).index, orient='h',
                                  color=sns.xkcd_rgb['dark blue'])
            plt.axvspan(0,1.3, color=sns.xkcd_rgb['medium blue'], alpha=0.2)
            ax.xaxis.grid()
            plt.title(ch)
            f=ax.get_figure()
            plt.xlim(0)
            ax.legend_.remove()
            f.savefig('./figures/frameshifts/frameshifts_nrpeaks1_'+'_'.join(ch.split(' ')) +'.png', \
                      dpi = 300, format='png', bbox_inches='tight', frameon=True)
    except:
        pass
  
for ch in df.changes.unique():
    try:
        with sns.axes_style('dark'):
            ax=df[(df.changes==ch)&(df.gfpp_wt<35)&(df.peaksp_wt==1)&(df.numberreadsp_wt>=10)&
                  ([x not in eventlist+nofslist+plist+[''] for x in df.fsevent])&
                  (df.laststopinframep1<12)].pivot_table(index='fsevent',values='gfpp_wt', 
                  aggfunc=np.median).sort_values('gfpp_wt', ascending=False).plot(kind='barh', 
                figsize=(4,len(df[(df.changes==ch)&(df.gfpp_wt<35)&(df.peaksp_wt==1)&
                (df.numberreadsp_wt>=10)&([x not in eventlist+nofslist+plist+[''] for x in df.fsevent])&(df.laststopinframep1<12)].pivot_table(index='fsevent',values='gfpp_wt'))*0.3))
            sns.swarmplot(data=df[(df.changes==ch)&(df.gfpp_wt<35)&(df.peaksp_wt==1)&
                    ([x not in eventlist+nofslist+plist+[''] for x in df.fsevent])&
                    (df.numberreadsp_wt>=10)&(df.laststopinframep1<12)], y='fsevent',x='gfpp_wt',
                     order=df[(df.changes==ch)&(df.gfpp_wt<35)&(df.peaksp_wt==1)&
                    ([x not in eventlist+nofslist+plist+[''] for x in df.fsevent])&
                    (df.numberreadsp_wt>=10)&(df.laststopinframep1<12)].pivot_table(
                    index='fsevent',values='gfpp_wt', 
                  aggfunc=np.median).sort_values('gfpp_wt', ascending=False).index, orient='h',
                                  color=sns.xkcd_rgb['dark blue'])
            plt.axvspan(0,1.3, color=sns.xkcd_rgb['medium blue'], alpha=0.2)
            ax.xaxis.grid()
            plt.xlim(0)
            plt.title(ch)
            f=ax.get_figure()
            ax.legend_.remove()
            f.savefig('./figures/frameshifts/frameshifts_nrpeaks1_p1_'+'_'.join(ch.split(' ')) +'.png', \
                  dpi = 300, format='png', bbox_inches='tight', frameon=True)
    except:
        pass


        
    
