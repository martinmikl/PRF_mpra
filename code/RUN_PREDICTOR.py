#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 29 16:18:33 2020

@author: martinm
"""

'''
first argument:
    name of the csv file containing two columns: 
        an identifier, 
        the DNA sequence of the region containing a potential frameshifting site
                (210 nt in length, starting in the original frame, slippery site ends after 72 nt
                e.g.: the first 72 nt end with TTTTTTA in the case of HIV)
second argument:
    name of output file (optional)
third argument:
    name of feature set used (optional)
        all (default)
        aa: only amino acid class around the slippery site
        dg: only minimum free energy of regions downstream of the frameshifting site
        sec: only information of predicted pairedness of downstream positions
        slip: only slippery site position identity
        tai: only tAI of regions around the slippery site (amino acid positions 12-25 
                    in frame 0 and frame -1 and the difference between them)
        
    
    
packages:
    pandas 0.23.4
    numpy 1.11.3
    RNA 2.4.11
    sklearn 0.18.1
    Bio 1.66


'''

import pandas as pd
import numpy as np
import sys
import pickle
import RNA

from sklearn.preprocessing import LabelEncoder
from sklearn.preprocessing import OneHotEncoder

from Bio.Seq import Seq

filename=sys.argv[1]

try:
    outputfilename=sys.argv[2]
except:
    outputfilename='./predictions'

try:
    featureset=sys.argv[3]
except:
    featureset='all'

#%%

topredict=pd.read_csv(filename, header=None)
topredict.columns=['identifier','varseq']
topredict['varseq162']=topredict.varseq.apply(lambda x: x[30:-18])

#%% make features

# tai

tAIdict={}

with open ('./forpredictor/tAI_index_human_nar-02315.txt') as infile:
    table = infile.read()
    table_list = table.splitlines()
    
for row in range (0, len(table_list)):
    linearr=str.split(table_list[row], '\t')
    if len(str(linearr[2]))==3:
        if linearr[4]!=' ----' and linearr[4]!='----':
            tAIdict[linearr[2]]=float(linearr[4])
        else:
            tAIdict[linearr[2]]=float(0)


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

tai0=topredict.varseq162.apply(lambda x: gettaimap(x.upper()))

taim1=topredict.varseq162.apply(lambda x: gettaimap(x[2:].upper()))

taip1=topredict.varseq162.apply(lambda x: gettaimap(x[1:].upper()))

taidiff=topredict.varseq162.apply(lambda x: gettaidiffmap(x.upper()))

tai=tai0.join(taim1, lsuffix='_frame0', rsuffix = '_framem1')
tai=tai.join(taip1, rsuffix = '_framep1')
tai=tai.join(taidiff, rsuffix = '_diffframem1to0')

# aa class

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
        aa.append(aacategory[str(Seq(vrsq[0+3*i:3+3*i]).translate())])
    return pd.Series(aa)

aa0=topredict.varseq162.apply(lambda x: getaaclassmap(x.upper()))

aam1=topredict.varseq162.apply(lambda x: getaaclassmap(x[2:].upper()))

aap1=topredict.varseq162.apply(lambda x: getaaclassmap(x[1:].upper()))


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

# paired

sec=topredict.varseq162.apply(lambda x: pd.Series(list(RNA.fold(x)[0])))
sec=sec.replace({'.':0,'(':1,')':1})

# deltaG

dg=pd.DataFrame()

dg['dg_downstream']=topredict.varseq162.apply(lambda x: RNA.fold(x[42:])[1])

for length in range(0,70,10):
    dg['dgup'+str(length)]=topredict.varseq.apply(lambda x: RNA.fold(x[length:30+42])[1])
for length in range(20,130,10):
    dg['dg'+str(length)]=topredict.varseq.apply(lambda x: RNA.fold(x[30+42:30+42+length])[1])

# slippery

def classifyslippery(varsq):
    slippery=list(varsq[35:42])
    if (len(pd.Series(slippery[0:3]).value_counts())==1) & (len(pd.Series(slippery[3:6]).value_counts())==1):
        pattern=varsq[35:42]
    elif (len(pd.Series(slippery[0:2]).value_counts())==1) & (len(pd.Series(slippery[2:6]).value_counts())==1):
        pattern=varsq[35:42]
    else:
        pattern='NNNNNNN'
    return pattern

topredict['slipclass']=topredict.varseq162.apply(lambda x: classifyslippery(x))

slip=topredict.slipclass.apply(lambda x: pd.Series([x[0],x[4],x[-1]], index=['first','second','last']))
slipenc=pd.DataFrame(index=slip.index)
for col in slip.columns:
    labelenc=LabelEncoder()
    feature0=labelenc.fit_transform(slip[col])
    onehotenc=OneHotEncoder(sparse=False)
    slipenc=slipenc.join(pd.DataFrame(onehotenc.fit_transform(feature0.reshape(feature0.shape[0],1)), index=slip.index, columns=labelenc.classes_), rsuffix=col)


#%%

taicols=[str(x)+'_frame0' for x in np.arange(12,26)]+[str(x)+'_framem1' for x in np.arange(12,26)]+[str(x)+'_diffframem1to0' for x in np.arange(12,26)]
aacols=['charged'+str(x) for x in np.arange(10, 14)]+['polar'+str(x) for x in np.arange(10, 14)]+['unpolar'+str(x) for x in np.arange(10, 14)]

seccols=np.arange(42,162)
dgcols=['dg_downstream',
 'dg20',
 'dg30',
 'dg40',
 'dg50',
 'dg60',
 'dg70',
 'dg80',
 'dg90',
 'dg100',
 'dg110']
slipcols=['A','C','G','N','T','Asecond','Csecond','Gsecond','Nsecond','Tsecond','Alast','Clast','Glast','Nlast','Tlast']

for a in aacols:
    if a not in aa.columns:
        aa[a]=0


for a in slipcols:
    if a not in slipenc.columns:
        slipenc[a]=0
          
feat=tai.join(aa)
feat=feat.join(dg)
feat=feat.join(sec)
feat=feat.join(slipenc)

#%%

r2scoresoptcondpeaks1min20=pd.read_csv('./forpredictor/r2scoresoptcondpeaks1min20.csv', index_col='Unnamed: 0')

'''
cols={'all':taicols+aacols+list(seccols)+list(dgcols)+list(slipcols),'tai':taicols,'aa':aacols,'sec':seccols,'dg':dgcols, 'slip':slipcols}
ypred=pd.DataFrame()
for event in r2scoresoptcondpeaks1min20.index:
    if (event=='eventlist'):
        clf=pickle.load(open('./ml/models/eventlist_'+featureset+'.mdl', 'rb'))
        ypred['based on PRF events']=pd.Series(clf.predict(feat[cols[featureset]]), index=topredict.identifier) 
    elif (event=='full'):
        clf=pickle.load(open('./ml/models/full_'+featureset+'.mdl', 'rb'))
        ypred['based on full library']=pd.Series(clf.predict(feat[cols[featureset]]), index=topredict.identifier) 
    else:
        clf=pickle.load(open('./ml/models/'+event+'_'+featureset+'.mdl', 'rb'))
        ypred['based on '+event+' variants']=pd.Series(clf.predict(feat[cols[featureset]]), index=topredict.identifier) 
'''


cols={'all':taicols+aacols+list(seccols)+list(dgcols)+list(slipcols),'tai':taicols,'aa':aacols,'sec':seccols,'dg':dgcols, 'slip':slipcols}
ypred=pd.DataFrame()
clf=pickle.load(open('./forpredictor/models/classifier_eventlist_'+featureset+'_m1.mdl', 'rb'))
ypred['predicted to frameshift?']=pd.Series(clf.predict(feat[cols[featureset]]), 
     index=topredict.identifier).replace(to_replace=0, value='no').replace(to_replace=1, value='yes')
ypred['predicted probability']=pd.Series(clf.predict_proba(feat[cols[featureset]])[:,1], 
     index=topredict.identifier)
clf=pickle.load(open('./forpredictor/models/eventlist_'+featureset+'_m1.mdl', 'rb'))
ypred['predicted %GFP']=pd.Series(clf.predict(feat[cols[featureset]]), index=topredict.identifier)

ypred.to_csv(outputfilename+'.csv')

#%%
'''
with open('./forpredictor/testinput_randomsequences.csv', 'w') as f:
    for j in range(100):
        f.write('test random' + str(i+1)+','+''.join([random.choice(['A','C','G','T']) for i in range(210)])+'\n')
    
'''
    