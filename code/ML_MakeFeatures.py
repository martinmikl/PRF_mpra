#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Aug  9 12:45:55 2018

@author: martinm
"""

         
#%%
import pandas as pd
import seaborn as sns
from Bio.Seq import Seq

sns.set_context('talk')

#%%
import RNA

from sklearn.preprocessing import LabelEncoder
from sklearn.preprocessing import OneHotEncoder

#%%
#%% 
fsdf=pd.read_pickle('./dataframes/fsdf.pkl')

fsdesforml=pd.read_pickle('./dataframes/ml/fsdesforml.pkl')
hivdbvarforml=pd.read_pickle('./dataframes/ml/hivdbvarforml.pkl')

#####


fsdesmlx=pd.read_pickle('./dataframes/ml/fsdesmlx.pkl')
fsdestestx=pd.read_pickle('./dataframes/ml/fsdestestx.pkl')
fsdesmly=pd.read_pickle('./dataframes/ml/fsdesmly.pkl')
fsdestesty=pd.read_pickle('./dataframes/ml/fsdestesty.pkl')
  

#####

hivmlx=pd.read_pickle('./dataframes/ml/hivmlx.pkl')
hivtestx=pd.read_pickle('./dataframes/ml/hivtestx.pkl')
hivmly=pd.read_pickle('./dataframes/ml/hivmly.pkl')
hivtesty=pd.read_pickle('./dataframes/ml/hivtesty.pkl')
      
#####



#%% tai

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

tai0=pd.DataFrame(index=fsdesforml.index)

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

tai0=fsdesforml.varseq162.apply(lambda x: gettaimap(x.upper()))

taim1=fsdesforml.varseq162.apply(lambda x: gettaimap(x[2:].upper()))

taip1=fsdesforml.varseq162.apply(lambda x: gettaimap(x[1:].upper()))

taidiff=fsdesforml.varseq162.apply(lambda x: gettaidiffmap(x.upper()))

tai=tai0.join(taim1, lsuffix='_frame0', rsuffix = '_framem1')
tai=tai.join(taip1, rsuffix = '_framep1')
tai=tai.join(taidiff, rsuffix = '_diffframem1to0')


tai.to_pickle('./ml/features/tai_features_fsdes.pkl')


### hiv

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



taihiv0=hivdbvarforml.varseq162.apply(lambda x: gettaimap(x.upper()))

taihivm1=hivdbvarforml.varseq162.apply(lambda x: gettaimap(x[2:].upper()))

taihivp1=hivdbvarforml.varseq162.apply(lambda x: gettaimap(x[1:].upper()))

taihivdiff=hivdbvarforml.varseq162.apply(lambda x: gettaidiffmap(x.upper()))

taihiv=taihiv0.join(taihivm1, lsuffix='_frame0', rsuffix = '_framem1')
taihiv=taihiv.join(taihivp1, rsuffix = '_framep1')
taihiv=taihiv.join(taihivdiff, rsuffix = '_diffframem1to0')

taihiv.to_pickle('./ml/features/tai_features_hiv.pkl')


#%% aa class

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
    '*':'stop'
    }


def getaaclassmap(vrsq):
    aa=[]
    for i in range(len(vrsq)/3-1):
        aa.append(aacategory[Seq(vrsq[0+3*i:3+3*i]).translate()])
    return pd.Series(aa)

aa0=fsdesforml.varseq162.apply(lambda x: getaaclassmap(x.upper()))

aam1=fsdesforml.varseq162.apply(lambda x: getaaclassmap(x[2:].upper()))

aap1=fsdesforml.varseq162.apply(lambda x: getaaclassmap(x[1:].upper()))


aa0enc=pd.DataFrame(index=aa0.index)
for col in aa0.columns:
    labelenc=LabelEncoder()
    feature0=labelenc.fit_transform(aa0[col])
    onehotenc=OneHotEncoder(sparse=False)
    aa0enc=aa0enc.join(pd.DataFrame(onehotenc.fit_transform(feature0.reshape(feature0.shape[0],1)), index=aa0.index, columns=labelenc.classes_), rsuffix=col)

aam1enc=pd.DataFrame(index=aam1.index)
for col in aam1.columns:
    labelenc=LabelEncoder()
    feature0=labelenc.fit_transform(aam1[col])
    onehotenc=OneHotEncoder(sparse=False)
    aam1enc=aam1enc.join(pd.DataFrame(onehotenc.fit_transform(feature0.reshape(feature0.shape[0],1)), index=aam1.index, columns=labelenc.classes_), rsuffix=col)

aap1enc=pd.DataFrame(index=aap1.index)
for col in aap1.columns:
    labelenc=LabelEncoder()
    feature0=labelenc.fit_transform(aap1[col])
    onehotenc=OneHotEncoder(sparse=False)
    aap1enc=aap1enc.join(pd.DataFrame(onehotenc.fit_transform(feature0.reshape(feature0.shape[0],1)), index=aap1.index, columns=labelenc.classes_), rsuffix=col)


aa=aa0enc.join(aam1enc, lsuffix='_frame0', rsuffix = '_framem1')
aa=aa.join(aap1enc, rsuffix = '_framep1')

aa.to_pickle('./ml/features/aa_allfeaturesenc_fsdes.pkl')


###

aa0=hivdbvarforml.varseq162.apply(lambda x: getaaclassmap(x.upper()))

aam1=hivdbvarforml.varseq162.apply(lambda x: getaaclassmap(x[2:].upper()))

aap1=hivdbvarforml.varseq162.apply(lambda x: getaaclassmap(x[1:].upper()))


aa0enc=pd.DataFrame(index=aa0.index)
for col in aa0.columns:
    labelenc=LabelEncoder()
    feature0=labelenc.fit_transform(aa0[col])
    onehotenc=OneHotEncoder(sparse=False)
    aa0enc=aa0enc.join(pd.DataFrame(onehotenc.fit_transform(feature0.reshape(feature0.shape[0],1)), index=aa0.index, columns=labelenc.classes_), rsuffix=col)

aam1enc=pd.DataFrame(index=aam1.index)
for col in aam1.columns:
    labelenc=LabelEncoder()
    feature0=labelenc.fit_transform(aam1[col])
    onehotenc=OneHotEncoder(sparse=False)
    aam1enc=aam1enc.join(pd.DataFrame(onehotenc.fit_transform(feature0.reshape(feature0.shape[0],1)), index=aam1.index, columns=labelenc.classes_), rsuffix=col)

aap1enc=pd.DataFrame(index=aap1.index)
for col in aap1.columns:
    labelenc=LabelEncoder()
    feature0=labelenc.fit_transform(aap1[col])
    onehotenc=OneHotEncoder(sparse=False)
    aap1enc=aap1enc.join(pd.DataFrame(onehotenc.fit_transform(feature0.reshape(feature0.shape[0],1)), index=aap1.index, columns=labelenc.classes_), rsuffix=col)


aahiv=aa0enc.join(aam1enc, lsuffix='_frame0', rsuffix = '_framem1')
aahiv=aahiv.join(aap1enc, rsuffix = '_framep1')

aahiv.to_pickle('./ml/features/aa_allfeaturesenc_hiv.pkl')


#%% paired

sec=fsdesforml.varseq162.apply(lambda x: pd.Series(list(RNA.fold(x)[0])))
sec=sec.replace({'.':0,'(':1,')':1})
sec.to_pickle('./ml/features/sec_features_fsdesforml.pkl')


###

sechiv=hivdbvarforml.varseq162.apply(lambda x: pd.Series(list(RNA.fold(x)[0])))
sechiv=sechiv.replace({'.':0,'(':1,')':1})
sechiv.to_pickle('./ml/features/sec_features_hivdbvarforml.pkl')

#%% slippery

def classifyslippery(varsq):
    slippery=list(varsq[35:42])
    if (len(pd.Series(slippery[0:3]).value_counts())==1) & (len(pd.Series(slippery[3:6]).value_counts())==1):
        pattern=varsq[35:42]
    elif (len(pd.Series(slippery[0:2]).value_counts())==1) & (len(pd.Series(slippery[2:6]).value_counts())==1):
        pattern=varsq[35:42]
    else:
        pattern='NNNNNNN'
    return pattern

fsdf['slipclass']=fsdf.varseq162.apply(lambda x: classifyslippery(x))

slip=fsdf.slipclass.apply(lambda x: pd.Series([x[0],x[4],x[-1]], index=['first','second','last']))
slipenc=pd.DataFrame(index=slip.index)
for col in slip.columns:
    labelenc=LabelEncoder()
    feature0=labelenc.fit_transform(slip[col])
    onehotenc=OneHotEncoder(sparse=False)
    slipenc=slipenc.join(pd.DataFrame(onehotenc.fit_transform(feature0.reshape(feature0.shape[0],1)), index=slip.index, columns=labelenc.classes_), rsuffix=col)

slipencdes=slipenc.loc[fsdesforml.index]
slipencdes.to_pickle('./ml/features/slipenc_fsdesforml.pkl')

slipenchiv=slipenc.loc[hivdbvarforml.index]
slipenchiv.to_pickle('./ml/features/slipenc_hivdbvarforml.pkl')


