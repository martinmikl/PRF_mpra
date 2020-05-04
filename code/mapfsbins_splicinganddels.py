#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 23 20:09:38 2017

@author: martinm
"""

'''
experiment is minus1, plus1, stopvector or RNA
The raw data (which can be obtained from GEO, accession GSE145684) should be split into portions (e.g. of 100,000 reads)
and mapping should be parallelized.
The sequencing data is paired end, and the splits from the two mates should be named 'Split1-' and a counter 
and 'Split2-' and a counter, respectively.
This file can be run from the code folder with
python2.7 mapfsbins_splicinganddels.py [filenumber/counter] [experiment: minus1, plus1, stopvector or RNA]
'''


import pandas as pd
from Bio import SeqIO
import sys

mapfolder='../rawdata/'

#%%
filenumber = sys.argv[1]
experiment = sys.argv[2]

fqread1 = SeqIO.to_dict(SeqIO.parse(mapfolder + experiment+ '/Split1-' + str(filenumber),'fastq'))
fqread2 = SeqIO.to_dict(SeqIO.parse(mapfolder + experiment+ '/Split2-' + str(filenumber),'fastq'))

deslift='GCCCCACGGAGGTGCCAC'
explift='CTCCCGGGCATGCGAATT'
upstream='GAGCTGTACAAGCCGGACCG'

sublib=pd.read_csv('./design/fslibraryall.csv', index_col='Unnamed: 0')
sublib['barcode']=sublib.varseq.apply(lambda x: x[18:30])

readsmap=pd.Series('',index=sublib.index)


for read in fqread1.keys():
    if (fqread1[read].seq.find(upstream)>-1):
        testbc=fqread1[read].seq[fqread1[read].seq.find(upstream)+21+18:fqread1[read].seq.find(upstream)+21+30]
        readsmap.loc[sublib[sublib.barcode==testbc].index]=readsmap.loc[sublib[sublib.barcode==testbc].index] +' '+ [str(fqread1[read].seq)+'_'+str(fqread2[read].seq.reverse_complement())]      
    elif (fqread2[read].seq.find(upstream)>-1):
        testbc=fqread2[read].seq[fqread2[read].seq.find(upstream)+21+18:fqread2[read].seq.find(upstream)+21+30]
        readsmap.loc[sublib[sublib.barcode==testbc].index]=readsmap.loc[sublib[sublib.barcode==testbc].index] +' '+ [str(fqread2[read].seq)+'_'+str(fqread1[read].seq.reverse_complement())]      

readsmap.to_pickle('./mapping/' + experiment + '/coveragePYTHON-' + str(filenumber) + '.pkl')

