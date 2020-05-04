#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 20 13:54:28 2019

@author: martinm
"""


import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import scipy.stats

sns.set_context('talk')


ratiocalcbins=pd.read_pickle('./mapping/minus1/ratiocalcbins1.pkl') 

xval=ratiocalcbins['median'].apply(lambda x: np.log2(x))

xvalstop=pd.read_pickle('./mapping/stopvector/xvals_corrected.pkl') 

# combined dataframes

fsdf=pd.read_pickle('./dataframes/fsdf.pkl')

wtvalues=pd.read_pickle('./dataframes/wtvalues.pkl')

# Plot Isolated clones vs FACSseq measurements

# Load data from isolated library clones tested in isolation

expisol=pd.read_pickle('./dataframes/isolatedHIVclones.pkl')


expisol['gfpm_wt']=expisol['first'].apply(lambda x: fsdf.gfpm_wt[x] if x in fsdf.index else np.nan)
expisol['numberreadsm_wt']=expisol['first'].apply(lambda x: fsdf.numberreadsm_wt[x] if x in fsdf.index else np.nan)
expisol['peaksm_wt']=expisol['first'].apply(lambda x: fsdf.peaksm_wt[x] if x in fsdf.index else np.nan)


f=plt.figure(figsize=(3,3))
plt.scatter(expisol[(expisol.gfpmeanperc>1.3)&(expisol.subset=='hivdbvar')&(expisol.gfpmeanperc<25)][['gfpm_wt','gfpmeanperc']].dropna().gfpm_wt,\
                    expisol[(expisol.gfpmeanperc>1.3)&(expisol.subset=='hivdbvar')&(expisol.gfpmeanperc<25)][['gfpm_wt','gfpmeanperc']].dropna().gfpmeanperc, 
                            color=sns.xkcd_rgb['medium blue'], alpha=0.7) 
plt.xlim(1,5)
plt.ylim(2,6)
f.savefig('./figures/hiv/gfpm_wt_vs_isolatedclones_abovebackground_smallerthan25.png', \
   dpi = 300, format='png', bbox_inches='tight', frameon=True)

len(expisol[(expisol.gfpmeanperc>1.3)&(expisol.subset=='hivdbvar')&(expisol.gfpmeanperc<25)][['gfpm_wt','gfpmeanperc']].dropna())
# n=69


scipy.stats.pearsonr(expisol[(expisol.gfpmeanperc>1.3)&(expisol.subset=='hivdbvar')&(expisol.gfpmeanperc<25)][['gfpm_wt','gfpmeanperc','mchmean']].dropna().gfpm_wt,\
                    expisol[(expisol.gfpmeanperc>1.3)&(expisol.subset=='hivdbvar')&(expisol.gfpmeanperc<25)][['gfpm_wt','gfpmeanperc','mchmean']].dropna().gfpmeanperc)
# (r=0.51911183650540027, p=4.8761201601090479e-06)



#%% HIV

# Correlation with secondary structure 

df=fsdf[(fsdf.subset=="hivdbvar")&(fsdf.peaksm_wt==1)&(fsdf.gfpm_wt<25)]
df=df.dropna(axis=1, how='all')

corrsdg=pd.DataFrame()
for col in [x for x in df.columns if ('dg' in x)&('up' not in x)]:
    corrsdg.loc[col.split('dg')[1], 'pearsonr']=scipy.stats.pearsonr(df[col], df['gfpm_wt'])[0]
    corrsdg.loc[col.split('dg')[1], 'pearsonp']=scipy.stats.pearsonr(df[col], df['gfpm_wt'])[1]

with sns.axes_style(style='dark'):
    f, ax1=plt.subplots(figsize=(5,3))
    ax1.plot(corrsdg.iloc[1:].loc[:,'pearsonr'], color=sns.xkcd_rgb['medium blue'])
    ax1.tick_params('y',colors=sns.xkcd_rgb['medium blue'])
    ax2=ax1.twinx()
    ax2.plot(np.log10(corrsdg.iloc[1:].loc[:,'pearsonp']), color=sns.xkcd_rgb['light blue'])
    ax2.tick_params('y',colors=sns.xkcd_rgb['light blue'])
    f.savefig('./figures/hiv/hivdb_pearsoncorrelation_gfp_dgvariouslengths.png', \
              dpi = 300, format='png', bbox_inches='tight', frameon=True)


# Plot distribution of HIV clinical isolates

f=plt.figure(figsize=(6,3))
df.gfpm_wt[df.gfpm_wt<12].hist(bins=50, linewidth=0)
plt.axvline(x=wtvalues.gfpm['HIV HXB2'], color='grey', linewidth=2)
f.savefig('./figures/hiv/hivdb_overallhistogram_gfp.png', \
  dpi = 300, format='png', bbox_inches='tight', frameon=True)

# len(df.gfpm_wt[df.gfpm_wt<12]) - 368

# Compare subtypes and country of origin

df[['subtype','place','year','ncbiaccession']]=df.id.apply(lambda x: pd.Series(x.split('.')[0:3]+x.split('.')[-1:], index=['subtype','place','year','ncbiaccession']))

df[df.year!='-'].year=df[df.year!='-'].year.astype(int)


scipy.stats.f_oneway(df[df.gfpm_wt<8].pivot_table(index='id',columns='subtype', values='gfpm_wt').loc[:, 'A1'].dropna(),
         df[df.gfpm_wt<8].pivot_table(index='id',columns='subtype', values='gfpm_wt').loc[:, '01_AE'].dropna(),
         df[df.gfpm_wt<8].pivot_table(index='id',columns='subtype', values='gfpm_wt').loc[:, '02_AG'].dropna(),
         df[df.gfpm_wt<8].pivot_table(index='id',columns='subtype', values='gfpm_wt').loc[:, 'B'].dropna(),
         df[df.gfpm_wt<8].pivot_table(index='id',columns='subtype', values='gfpm_wt').loc[:, 'C'].dropna(),
         df[df.gfpm_wt<8].pivot_table(index='id',columns='subtype', values='gfpm_wt').loc[:, 'D'].dropna())
#F_onewayResult(statistic=12.263086591150492, pvalue=7.442248030369468e-11)

scipy.stats.f_oneway(df[df.gfp<8].pivot_table(index='id',columns='place', values='gfpm_wt').loc[:, 'US'].dropna(),
         df[df.gfpm_wt<8].pivot_table(index='id',columns='place', values='gfpm_wt').loc[:, 'BW'].dropna(),
         df[df.gfpm_wt<8].pivot_table(index='id',columns='place', values='gfpm_wt').loc[:, 'KE'].dropna(),
         df[df.gfpm_wt<8].pivot_table(index='id',columns='place', values='gfpm_wt').loc[:, 'CN'].dropna(),
         df[df.gfpm_wt<8].pivot_table(index='id',columns='place', values='gfpm_wt').loc[:, 'ZA'].dropna(),

         df[df.gfpm_wt<8].pivot_table(index='id',columns='place', values='gfpm_wt').loc[:, 'ZM'].dropna(),
         df[df.gfpm_wt<8].pivot_table(index='id',columns='place', values='gfpm_wt').loc[:, 'BR'].dropna(),
         df[df.gfpm_wt<8].pivot_table(index='id',columns='place', values='gfpm_wt').loc[:, 'JP'].dropna(),
         df[df.gfpm_wt<8].pivot_table(index='id',columns='place', values='gfpm_wt').loc[:, 'VN'].dropna(),

         df[df.gfpm_wt<8].pivot_table(index='id',columns='place', values='gfpm_wt').loc[:, 'ES'].dropna(),
         df[df.gfpm_wt<8].pivot_table(index='id',columns='place', values='gfpm_wt').loc[:, 'UG'].dropna(),
         df[df.gfpm_wt<8].pivot_table(index='id',columns='place', values='gfpm_wt').loc[:, 'IN'].dropna(),
         df[df.gfpm_wt<8].pivot_table(index='id',columns='place', values='gfpm_wt').loc[:, 'CA'].dropna(),
         df[df.gfpm_wt<8].pivot_table(index='id',columns='place', values='gfpm_wt').loc[:, 'DK'].dropna(),

         df[df.gfpm_wt<8].pivot_table(index='id',columns='place', values='gfpm_wt').loc[:, 'TH'].dropna())
#F_onewayResult(statistic=3.0710095184379296, pvalue=0.00018529506118266932)

scipy.stats.kruskal(df[df.gfp<8].pivot_table(index='id',columns='place', values='gfp').loc[:, 'US'].dropna(),
         df[df.gfp<8].pivot_table(index='id',columns='place', values='gfp').loc[:, 'BW'].dropna(),
         df[df.gfp<8].pivot_table(index='id',columns='place', values='gfp').loc[:, 'KE'].dropna(),
         df[df.gfp<8].pivot_table(index='id',columns='place', values='gfp').loc[:, 'CN'].dropna(),
         df[df.gfp<8].pivot_table(index='id',columns='place', values='gfp').loc[:, 'ZA'].dropna(),

         df[df.gfp<8].pivot_table(index='id',columns='place', values='gfp').loc[:, 'ZM'].dropna(),
         df[df.gfp<8].pivot_table(index='id',columns='place', values='gfp').loc[:, 'BR'].dropna(),
         df[df.gfp<8].pivot_table(index='id',columns='place', values='gfp').loc[:, 'JP'].dropna(),
         df[df.gfp<8].pivot_table(index='id',columns='place', values='gfp').loc[:, 'VN'].dropna(),

         df[df.gfp<8].pivot_table(index='id',columns='place', values='gfp').loc[:, 'ES'].dropna(),
         df[df.gfp<8].pivot_table(index='id',columns='place', values='gfp').loc[:, 'UG'].dropna(),
         df[df.gfp<8].pivot_table(index='id',columns='place', values='gfp').loc[:, 'IN'].dropna(),
         df[df.gfp<8].pivot_table(index='id',columns='place', values='gfp').loc[:, 'CA'].dropna(),
         df[df.gfp<8].pivot_table(index='id',columns='place', values='gfp').loc[:, 'DK'].dropna(),

         df[df.gfp<8].pivot_table(index='id',columns='place', values='gfp').loc[:, 'TH'].dropna())
#KruskalResult(statistic=48.73976772610456, pvalue=9.9098004949197159e-06)

scipy.stats.kruskal(df[df.gfp<8].pivot_table(index='id',columns='subtype', values='gfp').loc[:, 'A1'].dropna(),
         df[df.gfp<8].pivot_table(index='id',columns='subtype', values='gfp').loc[:, '01_AE'].dropna(),
         df[df.gfp<8].pivot_table(index='id',columns='subtype', values='gfp').loc[:, '02_AG'].dropna(),
         df[df.gfp<8].pivot_table(index='id',columns='subtype', values='gfp').loc[:, 'B'].dropna(),
         df[df.gfp<8].pivot_table(index='id',columns='subtype', values='gfp').loc[:, 'C'].dropna(),
         df[df.gfp<8].pivot_table(index='id',columns='subtype', values='gfp').loc[:, 'D'].dropna())
#KruskalResult(statistic=43.812016674374711, pvalue=2.5289221135398068e-08)

scipy.stats.ttest_ind(df[df.gfpm_wt<8].pivot_table(index='id',columns='subtype', values='gfpm_wt').loc[:, 'B'].dropna(),
         df[df.gfpm_wt<8].pivot_table(index='id',columns='subtype', values='gfpm_wt').loc[:, 'C'].dropna())
#Ttest_indResult(statistic=-7.3040113953524575, pvalue=5.5035295977024532e-12)

f=plt.figure(figsize=(4,3))
sns.boxplot(data=df[df.gfpm_wt<25], x='subtype',y='gfpm_wt', order=['A1','01_AE','02_AG','B','C','D'], color=sns.xkcd_rgb['light blue'])
ax=sns.swarmplot(data=df[df.gfpm_wt<25], x='subtype',y='gfpm_wt', order=['A1','01_AE','02_AG','B','C','D'], color=sns.xkcd_rgb['medium blue'], alpha=0.7, size=4)
plt.ylim(0,6)
ax.set_xticklabels(['A1','AE','AG','B','C','D'])
f.savefig('./figures/hiv/hivdb_by_subtype.png', \
  dpi = 300, format='png', bbox_inches='tight', frameon=True)

df[df.gfpm_wt<25].pivot_table(index='subtype',values='gfpm_wt',aggfunc=np.count_nonzero).loc[['A1','01_AE','02_AG','B','C','D']].to_csv('./n_numbers/Fig6c.csv')

f=plt.figure(figsize=(8,3))
sns.boxplot(data=df[df.gfpm_wt<25], x='place',y='gfpm_wt', order=df.place.value_counts().index[0:15], color=sns.xkcd_rgb['light blue'])
ax=sns.swarmplot(data=df[df.gfpm_wt<25], x='place',y='gfpm_wt', order=df.place.value_counts().index[0:15], color=sns.xkcd_rgb['medium blue'], alpha=0.7, size=4)
plt.ylim(0,6)
f.savefig('./figures/hiv/hivdb_by_place.png', \
  dpi = 300, format='png', bbox_inches='tight', frameon=True)

df[df.gfpm_wt<25].pivot_table(index='place',values='gfpm_wt',aggfunc=np.count_nonzero).to_csv('./n_numbers/FigS12c.csv')

df['decade']=df[df.year!='-'].year.apply(lambda x: str(x)[0:3])

scipy.stats.f_oneway(df[df.gfpm_wt<8].pivot_table(index='id',columns='decade', values='gfpm_wt').loc[:, '198'].dropna(),
         df[df.gfpm_wt<8].pivot_table(index='id',columns='decade', values='gfpm_wt').loc[:, '199'].dropna(),
         df[df.gfpm_wt<8].pivot_table(index='id',columns='decade', values='gfpm_wt').loc[:, '200'].dropna(),
         df[df.gfpm_wt<8].pivot_table(index='id',columns='decade', values='gfpm_wt').loc[:, '201'].dropna())
#F_onewayResult(statistic=0.99297547219173243, pvalue=0.39635063924458469)


f=plt.figure(figsize=(4,3))
sns.boxplot(data=df[(df.gfpm_wt<8)&(df.year!='-')], x='decade',y='gfpm_wt', order=df[(df.gfpm_wt<8)&(df.year!='-')].decade.sort_values().unique()[1:], color=sns.xkcd_rgb['light blue'])
ax=sns.swarmplot(data=df[(df.gfpm_wt<8)&(df.year!='-')], x='decade',y='gfpm_wt', order=df[(df.gfpm_wt<8)&(df.year!='-')].decade.sort_values().unique()[1:], color=sns.xkcd_rgb['medium blue'], alpha=0.7, size=4)
plt.ylim(0,6)
ax.set_xticklabels(df[(df.gfpm_wt<8)&(df.year!='-')].decade.sort_values().unique()[1:], rotation=45, horizontalalignment='right')
f.savefig('./figures/hiv//hivdb_by_decade.png', \
  dpi = 300, format='png', bbox_inches='tight', frameon=True)


df[df.gfpm_wt<25].pivot_table(index='decade',values='gfpm_wt',aggfunc=np.count_nonzero).to_csv('./n_numbers/FigS12d.csv')


##### compare with viral load (from Stanford HIV-DB)

titer=pd.read_table('../additional/viralload.txt')
titer=titer.rename(columns={'Accession':'ncbiaccession'})
titer['varseq162']=titer.Sequence.apply(lambda x: x[x.upper().find('TAATTTTTTA')-32:x.upper().find('TAATTTTTTA')+130].upper())

for i in titer.index:
    if len(titer.loc[i, 'varseq162'])<160:
        titer.drop(i, inplace=True)

titer['libindex']=titer.varseq162.apply(lambda x: df[[x[20:-60] in varsq for varsq in df.varseq162]].index.values)

titer['lenlibindex']=titer.libindex.apply(lambda x: len(x))

for i in titer[titer.lenlibindex>0].index:
    df.loc[titer.libindex[i][0], 'viralload']=int(titer.at[i, 'Viral load'])

meandiff=df[df.viralload>0].sort_values(by='gfpm_wt')
meandiff=df[df.viralload>0].sort_values(by='viralload')
#


f=plt.figure(figsize=(4,3))
#plt.scatter(df.gfpm_wt, df.viralload, color=sns.xkcd_rgb['medium blue'])
plt.scatter(meandiff.gfpm_wt, meandiff.viralload.rolling(10, center=True).mean(), color=sns.xkcd_rgb['medium blue'])
plt.axvline(x=wtvalues.gfpm['HIV HXB2'], color='grey', linewidth=2)
plt.xlim(0,8)
plt.ylim(-1,1000000)
f.savefig('./figures/hiv/hivdb_gfp_vs_virusload_10_to_minus60.png', \
  dpi = 300, format='png', bbox_inches='tight', frameon=True)

len(meandiff[(meandiff.gfpm_wt<8)&(meandiff.viralload.rolling(10, center=True).mean()<10000000)]) # 77
