#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sun Jun 16 19:52:10 2019

@author: martinm
"""

import pandas as pd
import numpy as np
from Bio import SeqIO
import difflib
from Bio.Seq import Seq


def classifyread_ir(intron):
    if (intron.loc['intronstart']=='0') & (intron.loc['intronend']=='0'):
        return 'wt'
    elif (intron.loc['intlen']==0):
        if len(intron.loc['mutation'])==2:
            return intron.loc['mutation'][0]+str(int(intron.loc['intronend'])-1)+intron.loc['mutation'][1]
        else:
            return 'X'+str(int(intron.loc['intronend'])-1)+'X'
    elif (intron.loc['intlen']==300):
        try:
            return intron.loc['intronstart'].split('_')[0] + 'mut' + intron.loc['intronend'].split('_')[0] + '|' + intron.loc['intronstart'].split('_')[1] + 'mut' + intron.loc['intronend'].split('_')[1]
        except:
            return 'undetermined'
    elif (intron.loc['intlen']<11):
        return str(intron.loc['intronstart'])+'del'+str(intron.loc['intronend'])
    elif (intron.loc['intlen']>10):
        return str(intron.loc['intronstart'])+'spliced'+str(intron.loc['intronend'])
    else:
        return 'undefined'

def getbinbarcodes(selectedreads, rawreads):
    binbarcodes=[]
    for sel in selectedreads:
        if len(sel.split('_'))==2:
            for raw in rawreads:
                if raw.find(sel.split('_')[1])>-1:
                    pos=raw.find(sel.split('_')[0])
                    if pos>7:
                        binbarcodes.append(raw[pos-8:pos])
    return ' '.join(binbarcodes)
            

def updateclassifier(ident):
    if '|wt|' in ident:
        if len(ident.split('_'))>1:
            identnew=[x for x in ident.split('_') if '|wt|' not in x]
            if len(identnew)==0:
                return 'wt'
            else:
                return '_'.join(identnew)
    else:
        return ident
                            
   
    
def map_intronsdels_fs_trial(rnareads, thresholdreads, startvariant, endvariant, kind='RNA', fstype='m1'):
    upstream='GAGCTGTACAAGCCGGACCG'
    if kind=='stop':
        upstream='GAGCTGTACAAGCCGGACTGATAGCTGACTAGTCGGACCGA'
    library=pd.read_pickle('/net/mraid08/export/genie/Runs/Martin/Library/alllibraries210.pkl')
    lib300=library[(library.library=='fs_designed')|(library.library=='fs_explorative')]
    lib300['barcode']=lib300.varseq.apply(lambda x: x[18:30])
    rna=pd.DataFrame()
    
    if kind=='bins':
        primers = SeqIO.to_dict(SeqIO.parse('/net/mraid08/export/genie/Runs/Martin/fsdesbins/primersmChRsr.txt',"tab"))
        binbc={}
        for x in primers.keys():
            binbc[str(primers[x].seq[6:14])]=x
        colnames=[]
        for i in range(1,17):
            for j in range(1,4):
                colnames.append('Bin#'+str(i)+'_'+str(j)+'_Fw')
        cov=pd.DataFrame(columns=colnames)

    if kind=='stop':       
        primers = {1:'HNHNHNCTGTGAAGgagctgtacaagccgGACtga',
        2	:'HNHNHNTGAGGTGGgagctgtacaagccgGACtga',
        3	:'HNHNHNTCCAAGGCgagctgtacaagccgGACtga',
        4	:'HNHNHNGTGAGTACgagctgtacaagccgGACtga',
        5	:'HNHNHNTCATCGGAgagctgtacaagccgGACtga',
        6	:'HNHNHNGCTTCTTCgagctgtacaagccgGACtga',
        7	:'HNHNHNTTAGTTCAgagctgtacaagccgGACtga',
        8	:'HNHNHNCCAGTATGgagctgtacaagccgGACtga',
        9	:'HNHNHNGCAGTCGAgagctgtacaagccgGACtga',
        10	:'HNHNHNAGGTATTGgagctgtacaagccgGACtga',
        11	:'HNHNHNACATGGAGgagctgtacaagccgGACtga',
        12	:'HNHNHNTTGAAGATgagctgtacaagccgGACtga'}
        
        binbc={}
        for x in primers.keys():
            binbc[primers[x][6:14]]='Bin'+str(x)
        
        colnames=[]
        for i in range(1,13):
            colnames.append('Bin'+str(i))
        cov=pd.DataFrame(columns=colnames)


    for var in range(startvariant, endvariant):
        varsq=lib300.varseq[var].upper()
        intronstarts=[]
        intronends=[]
        mutation=[]
        intronstarts1=[]
        intronends1=[]
        intronstarts2=[]
        intronends2=[]
        mutation1=[]
        mutation2=[]
        reads=[]
        readcount=[]
        rr_raw=pd.Series(rnareads.loc[var].split(' '))
        rr_raw.drop(0, inplace=True)
        rr_raw=rr_raw[[len(rr_raw.loc[x].split('_'))==2 for x in rr_raw.index]]
        if kind=='bins':
            if fstype=='m1':
                rr_trimmed=rr_raw.apply(lambda x: x[x.find(upstream):x.find('TTAGCTTAGGCGCGCC')] if x.find(upstream)>7 else 'too_short')
            elif fstype=='p1':
                rr_trimmed=rr_raw.apply(lambda x: x[x.find(upstream):x.find('TTCATTGGGCGCGCC')] if x.find(upstream)>7 else 'too_short')
        else:
            rr_trimmed=rr_raw.apply(lambda x: x[x.find(upstream):x.find('TTAGCTTAGGCGCGCC')] if 'TTAGCTTAGGCGCGCC' in x else x[x.find(upstream):x.find('TTCATTGGGCGCGCC')])
        rr_trimmed.drop(rr_trimmed[rr_trimmed=='too_short'].index, inplace=True)
        # do it with try because it always raises an error (index out of range) and can't figure out why
        try:
            if kind=='stop':
                rr_short=rr_trimmed.apply(lambda x: x.split('_')[0][:110] + '_' + x.split('_')[1][-90:] if (len(x.split('_')[0])>129)&(len(x.split('_')[1])>89) else x)
            else:
                rr_short=rr_trimmed.apply(lambda x: x.split('_')[0][:130] + '_' + x.split('_')[1][-90:] if (len(x.split('_')[0])>129)&(len(x.split('_')[1])>89) else x)
            rr=rr_short.value_counts()
        except:
            rr=rr_trimmed.value_counts()
            
        if len(rr[rr>=thresholdreads])>0:
            for test_full in rr[rr>=thresholdreads].index:
                if len(test_full.split('_'))==2:
                    test1=test_full.split('_')[0]
                    vspos=0
                    startpos=len(test1)-10
                    endpos=len(test1)
                    testsub=test1[startpos:endpos]
                    if varsq.find(testsub)==-1:
                        while True:
                            startpos-=1
                            endpos-=1
                            testsub=test1[startpos:endpos]
                            if (varsq.find(testsub)>-1)|(startpos<100):
                                break
                    while True:
                #        print(startpos)
                        vspos=varsq.find(testsub)
                        startpos-=1
                        testsub=test1[startpos:endpos]
                        if (varsq.find(testsub)==-1)|(startpos<0):
                            if (startpos<10):
                                vspos=0
                            break
                    intronends1.append(vspos)
                    if (vspos==0):
                        intronstarts1.append(0)
                        mutation1.append('')
                    elif (vspos==-1):
                        intronstarts1.append(-150)
                        mutation1.append('')
                    elif varsq.find(test1[startpos-10:startpos])==vspos-11:                    
                        intronstarts1.append(vspos)
                        mutation1.append(varsq[vspos-1]+test1[startpos])
                    else:
                        intronstarts1.append(varsq.find(test1[startpos-9:startpos+1])+10)
                        mutation1.append('')
    
                    test2=test_full.split('_')[1]
                    vspos=0
                    startpos=len(test2)-10
                    endpos=len(test2)
                    testsub=test2[startpos:endpos]
                    if varsq.find(testsub)==-1:
                        while True:
                            startpos-=1
                            endpos-=1
                            testsub=test2[startpos:endpos]
                            if (varsq.find(testsub)>-1)|(startpos<100):
                                break
                    while True:
                #        print(startpos)
                        vspos=varsq.find(testsub)
                        startpos-=1
                        testsub=test2[startpos:endpos]
                        if (varsq.find(testsub)==-1)|(startpos<0):
                            if (startpos<10):
                                vspos=0
                            break
                    intronends2.append(vspos)
                    if (vspos==0):
                        intronstarts2.append(0)
                        mutation2.append('')
                    elif (vspos==-1):
                        intronstarts2.append(-150)
                        mutation2.append('')
                    elif varsq.find(test2[startpos-10:startpos])==vspos-11:                    
                        intronstarts2.append(vspos)
                        mutation2.append(varsq[vspos-1]+test2[startpos])
                    else:
                        intronstarts2.append(varsq.find(test2[startpos-9:startpos+1])+10)
                        mutation2.append('')

            for i in range(len(intronstarts1)):
                if (intronstarts1[i]==intronstarts2[i])&(intronends1[i]==intronends2[i]):
                    intronstarts.append(str(intronstarts1[i]))
                    intronends.append(str(intronends1[i]))
                    mutation.append(mutation1[i])
                    reads.append(rr.keys()[i])
                    readcount.append(rr.values[i])
                elif ((intronstarts1[i]==0)|(intronstarts2[i]==0))&((intronends1[i]==0)|(intronends2[i]==0)):
                    intronstarts.append(str(np.max([intronstarts1[i],intronstarts2[i]])))
                    intronends.append(str(np.max([intronends1[i],intronends2[i]])))
                    mutation.append(mutation1[i]+mutation2[i])
                    reads.append(rr.keys()[i])
                    readcount.append(rr.values[i])
                else:
                    intronstarts.append(str(intronstarts1[i])+'_'+str(intronstarts2[i]))
                    intronends.append(str(intronends1[i])+'_'+str(intronends2[i]))
                    mutation.append(mutation1[i]+mutation2[i])
                    reads.append(rr.keys()[i])
                    readcount.append(rr.values[i])
                    
            
            if (len(intronstarts)>0)&(len(intronends)>0):
                introns=pd.DataFrame([intronstarts,intronends, mutation, reads, readcount]).transpose()
                introns.columns=['intronstart','intronend', 'mutation', 'read', 'readcount']
                
                introns['intlen']=introns.index.map(lambda x: int(introns.intronend[x])-int(introns.intronstart[x]) if len(introns.intronend[x])<4 else 300)
                
                introns['identifier']=introns.apply(lambda x: classifyread_ir(x), axis=1)
                
                for ident in introns.identifier.unique():
                    intron=introns[introns.identifier==ident]
                    rna.loc[str(var)+'_'+ident,'libindex']=var
                    rna.loc[str(var)+'_'+ident,'variant']=ident
                    rna.loc[str(var)+'_'+ident,'read']=intron.read.values[0]
                    rna.loc[str(var)+'_'+ident,'readcount']=intron.readcount.values[0]
                    rna.loc[str(var)+'_'+ident,'numberreads']=intron.readcount.sum()
                    if (kind=='bins')|(kind=='stop'):
                        rna.loc[str(var)+'_'+ident,'binbarcodes']=getbinbarcodes(intron.read.values, rr_raw.values)
                        rna.loc[str(var)+'_'+ident,'binbarcodecounts']=len(rna.loc[str(var)+'_'+ident,'binbarcodes'].split(' '))
                    
                        cov.loc[str(var)+'_'+ident]=0
                        for barcode in rna.loc[str(var)+'_'+ident,'binbarcodes'].split(' '): 
                            try:
                                cov.loc[str(var)+'_'+ident, binbc[barcode]]+=1
                            except:
                                pass        
    if (kind=='bins')|(kind=='stop'):
        return rna, cov
    else:
        return rna
    
    
def map_intronsdels_fs_difflib(rnareads, thresholdreads, startvariant, endvariant, kind='RNA', fstype='m1'):
    upstream='GAGCTGTACAAGCCGGACCG'
    if kind=='stop':
        upstream='GAGCTGTACAAGCCGGACTGATAGCTGACTAGTCGGACCGA'
    library=pd.read_pickle('/net/mraid08/export/genie/Runs/Martin/Library/alllibraries210.pkl')
    lib300=library[(library.library=='fs_designed')|(library.library=='fs_explorative')]
    lib300['barcode']=lib300.varseq.apply(lambda x: x[18:30])
    rna=pd.DataFrame()
    
    if kind=='bins':
        primers = SeqIO.to_dict(SeqIO.parse('/net/mraid08/export/genie/Runs/Martin/fsdesbins/primersmChRsr.txt',"tab"))
        binbc={}
        for x in primers.keys():
            binbc[str(primers[x].seq[6:14])]=x
        colnames=[]
        for i in range(1,17):
            for j in range(1,4):
                colnames.append('Bin#'+str(i)+'_'+str(j)+'_Fw')
        cov=pd.DataFrame(columns=colnames)

    if kind=='stop':       
        primers = {1:'HNHNHNCTGTGAAGgagctgtacaagccgGACtga',
        2	:'HNHNHNTGAGGTGGgagctgtacaagccgGACtga',
        3	:'HNHNHNTCCAAGGCgagctgtacaagccgGACtga',
        4	:'HNHNHNGTGAGTACgagctgtacaagccgGACtga',
        5	:'HNHNHNTCATCGGAgagctgtacaagccgGACtga',
        6	:'HNHNHNGCTTCTTCgagctgtacaagccgGACtga',
        7	:'HNHNHNTTAGTTCAgagctgtacaagccgGACtga',
        8	:'HNHNHNCCAGTATGgagctgtacaagccgGACtga',
        9	:'HNHNHNGCAGTCGAgagctgtacaagccgGACtga',
        10	:'HNHNHNAGGTATTGgagctgtacaagccgGACtga',
        11	:'HNHNHNACATGGAGgagctgtacaagccgGACtga',
        12	:'HNHNHNTTGAAGATgagctgtacaagccgGACtga'}
        
        binbc={}
        for x in primers.keys():
            binbc[primers[x][6:14]]='Bin'+str(x)
        
        colnames=[]
        for i in range(1,13):
            colnames.append('Bin'+str(i))
        cov=pd.DataFrame(columns=colnames)


    for var in range(startvariant, endvariant):

        varsq=lib300.varseq[var].upper()
        alal=difflib.SequenceMatcher(None)
        alal.set_seq1(varsq)
        intronstarts=[]
        intronends=[]
        mutation=[]
        intronstarts1=[]
        intronends1=[]
        intronstarts2=[]
        intronends2=[]
        mutation1=[]
        mutation2=[]
        classifier=[]
        reads=[]
        rr_raw=pd.Series(rnareads.loc[var].split(' '))        
        rr_raw.drop(0, inplace=True)

        if len(rr_raw)>1:
            for test_full in rr_raw.values:
                if len(test_full.split('_'))==2:
                    reads.append(test_full)
                    test1=test_full.split('_')[0]
                    alal.set_seq2(test1)    
                    matches=alal.get_matching_blocks()
                    if (len(matches)==2)&(matches[0][1]+matches[0][2]>144):
                        intronstarts1.append([])
                        intronends1.append([])
                        mutation1.append([])
                    elif len(matches)==1:
                        intronstarts1.append([-150])
                        intronends1.append([-1])
                        mutation1.append(['unknown'])
                    else: 
                        intronstarts1.append([matches[x][0]+matches[x][2] for x in range(len(matches)-1) if matches[x+1][1]<145])
                        intronends1.append([matches[x+1][0] - (matches[x+1][1]-(matches[x][2]+matches[x][1])) for x in range(len(matches)-1) if matches[x+1][1]<145])
                        if (intronstarts1[-1]==[])|(intronends1[-1]==[]):
                            mutation1.append([])
                        else:
                            mutation1.append(['del' if intronstarts1[-1][x]!=intronends1[-1][x] else 
                                          varsq[intronstarts1[-1][x]:intronstarts1[-1][x]+(matches[x+1][1]-(matches[x][2]+matches[x][1]))]
                                          +'mut'+ test1[matches[x][2]+matches[x][1]:matches[x+1][1]] 
                                          for x in range(len(intronstarts1[-1]))])

                    test2=test_full.split('_')[1]
                    alal.set_seq2(test2)    
                    matches=alal.get_matching_blocks()
                    if (len(matches)==2)&(matches[0][0]+matches[0][2]>204):
                        intronstarts2.append([0])
                        intronends2.append([0])
                        mutation2.append(['wt'])
                    elif len(matches)==1:
                        intronstarts2.append([-150])
                        intronends2.append([-1])
                        mutation2.append(['unknown'])
                    else: 
                        intronstarts2.append([matches[x][0]+matches[x][2] for x in range(len(matches)-2) if matches[x+1][1]>4])
                        intronends2.append([matches[x+1][0] - (matches[x+1][1]-(matches[x][2]+matches[x][1])) for x in range(len(matches)-2) if matches[x+1][1]>4])
                        if (intronstarts2[-1]==[])|(intronends2[-1]==[]):
                            mutation2.append(['wt'])
                        else:
                            mutation2.append(['del' if intronstarts2[-1][x]!=intronends2[-1][x] else 
                                          varsq[intronstarts2[-1][x]:intronstarts2[-1][x]+(matches[x+1][1]-(matches[x][2]+matches[x][1]))]
                                          +'mut'+ test1[matches[x][2]+matches[x][1]:matches[x+1][1]] 
                                          for x in range(len(intronstarts2[-1]))])
                                        

            for i in range(len(intronstarts1)):
                if ((intronstarts1[i]==[])&(intronstarts2[i]==[]))&((intronends1[i]==[])&(intronends2[i]==[])):
                    intronstarts.append(0)
                    intronends.append(0)
                    mutation.append([])
                    classifier.append('wt')
#                    mutation.append(mutation1[i]+mutation2[i])
#                    reads.append(rr.keys()[i])
#                    readcount.append(rr.values[i])
                else:
                    intronstarts.append(intronstarts1[i]+intronstarts2[i])
                    intronends.append(intronends1[i]+intronends2[i])
                    mutation.append(mutation1[i]+mutation2[i])
                    classifier.append('_'.join([str(intronstarts[i][x])+mutation[i][x]+str(intronends[i][x]) for x in range(len(intronstarts[i]))]))
                    if ('del' not in classifier[-1])&('mut' not in classifier[-1])&('wt' in classifier[-1]):
                        classifier[-1]='wt'
#                    reads.append(rr.keys()[i])
#                    readcount.append(rr.values[i])
                    
            
            if (len(intronstarts)>0)&(len(intronends)>0):
                introns=pd.DataFrame([intronstarts,intronends, mutation, reads,classifier]).transpose()
                introns.columns=['intronstart','intronend', 'mutation', 'read','classifier']
                introns['identifier']=introns.classifier.apply(lambda x: updateclassifier(x))
                introns['barcode']=introns.read.apply(lambda x: x[x.find(upstream)-8:x.find(upstream)])
                                
                for ident in introns.identifier.value_counts()[introns.identifier.value_counts()>thresholdreads].index:
                    intron=introns[introns.identifier==ident]
                    rna.loc[str(var)+'_'+ident,'libindex']=var
                    rna.loc[str(var)+'_'+ident,'variant']=ident
                    rna.loc[str(var)+'_'+ident,'read']=intron.read.values[0]
                    rna.loc[str(var)+'_'+ident,'numberreads']=len(intron)
                    if (kind=='bins')|(kind=='stop'):
                        rna.loc[str(var)+'_'+ident,'binbarcodes']=' '.join(intron.barcode)
                        rna.loc[str(var)+'_'+ident,'binbarcodecounts']=len(rna.loc[str(var)+'_'+ident,'binbarcodes'].split(' '))
                    
                        cov.loc[str(var)+'_'+ident]=0
                        for barcode in rna.loc[str(var)+'_'+ident,'binbarcodes'].split(' '): 
                            try:
                                cov.loc[str(var)+'_'+ident, binbc[barcode]]+=1
                            except:
                                pass        
    if (kind=='bins')|(kind=='stop'):
        return rna, cov
    else:
        return rna
    

def map_intronsdels_fs_difflibupdown(rnareads, thresholdreads, startvariant, endvariant, kind='RNA', fstype='m1'):
    upstream='GAGCTGTACAAGCCGGACCGA'
    if kind=='stop':
        upstream='GAGCTGTACAAGCCGGACTGATAGCTGACTAGTCGGACCGA'
    downstream='TTAGCTTAGGCGCGCCTGGTAGCTCTAGAGCTCGGACGGGTGCGCTC'
    if fstype=='p1':
        downstream='TTCATTGGGCGCGCCTGGTAGCTCTAGAGCTCGGACGGGTGCGCTC'
    if kind=='RNA':
        downstream='TT'
    library=pd.read_pickle('/net/mraid08/export/genie/Runs/Martin/Library/alllibraries210.pkl')
    lib300=library[(library.library=='fs_designed')|(library.library=='fs_explorative')]
    lib300['barcode']=lib300.varseq.apply(lambda x: x[18:30])
    rna=pd.DataFrame()
    
    if kind=='bins':
        primers = SeqIO.to_dict(SeqIO.parse('/net/mraid08/export/genie/Runs/Martin/fsdesbins/primersmChRsr.txt',"tab"))
        binbc={}
        for x in primers.keys():
            binbc[str(primers[x].seq[6:14])]=x
        colnames=[]
        for i in range(1,17):
            for j in range(1,4):
                colnames.append('Bin#'+str(i)+'_'+str(j)+'_Fw')
        cov=pd.DataFrame(columns=colnames)

    if kind=='stop':       
        primers = {1:'HNHNHNCTGTGAAGgagctgtacaagccgGACtga',
        2	:'HNHNHNTGAGGTGGgagctgtacaagccgGACtga',
        3	:'HNHNHNTCCAAGGCgagctgtacaagccgGACtga',
        4	:'HNHNHNGTGAGTACgagctgtacaagccgGACtga',
        5	:'HNHNHNTCATCGGAgagctgtacaagccgGACtga',
        6	:'HNHNHNGCTTCTTCgagctgtacaagccgGACtga',
        7	:'HNHNHNTTAGTTCAgagctgtacaagccgGACtga',
        8	:'HNHNHNCCAGTATGgagctgtacaagccgGACtga',
        9	:'HNHNHNGCAGTCGAgagctgtacaagccgGACtga',
        10	:'HNHNHNAGGTATTGgagctgtacaagccgGACtga',
        11	:'HNHNHNACATGGAGgagctgtacaagccgGACtga',
        12	:'HNHNHNTTGAAGATgagctgtacaagccgGACtga'}
        
        binbc={}
        for x in primers.keys():
            binbc[primers[x][6:14]]='Bin'+str(x)
        
        colnames=[]
        for i in range(1,13):
            colnames.append('Bin'+str(i))
        cov=pd.DataFrame(columns=colnames)


    for var in range(startvariant, endvariant):

        varsq=upstream + lib300.varseq[var].upper() + downstream
        alal=difflib.SequenceMatcher(None)
        alal.set_seq1(varsq)
        intronstarts=[]
        intronends=[]
        mutation=[]
        intronstarts1=[]
        intronends1=[]
        intronstarts2=[]
        intronends2=[]
        mutation1=[]
        mutation2=[]
        classifier=[]
        reads=[]
        rr_raw=pd.Series(rnareads.loc[var].split(' '))        
        rr_raw.drop(0, inplace=True)

        if len(rr_raw)>1:
            for test_full in rr_raw.values:
                if len(test_full.split('_'))==2:
                    reads.append(test_full)
                    test1=test_full.split('_')[0]
                    alal.set_seq2(test1)    
                    matches=alal.get_matching_blocks()
                    if (len(matches)==2)&(matches[0][1]+matches[0][2]>144):
                        intronstarts1.append([])
                        intronends1.append([])
                        mutation1.append([])
                    elif len(matches)==1:
                        intronstarts1.append([-150])
                        intronends1.append([-1])
                        mutation1.append(['unknown'])
                    else: 
#                        intronstarts1.append([matches[x][0]+matches[x][2] for x in range(len(matches)-1) if matches[x+1][1]<145])
#                        intronends1.append([matches[x+1][0] - (matches[x+1][1]-(matches[x][2]+matches[x][1])) for x in range(len(matches)-1) if matches[x+1][1]<145])
                        intronstarts1.append([matches[x][0]+matches[x][2] for x in range(len(matches)-2) if (matches[x+1][1]<len(test1)-9)&(matches[x][2]>9)])
                        intronends1.append([matches[x+1][0] - (matches[x+1][1]-(matches[x][2]+matches[x][1])) for x in range(len(matches)-2) if (matches[x+1][1]<len(test1)-9)&(matches[x][2]>9)])
                        if (intronstarts1[-1]==[])|(intronends1[-1]==[]):
                            mutation1.append([])
                        else:
                            muts=[]
                            for x in range(len(intronstarts1[-1])):
                                if intronstarts1[-1][x]!=intronends1[-1][x]:
                                    muts.append('del')
                                else:
                                    if (matches[x+1][1]>9)&(matches[x][2]>9):
                                        muts.append(varsq[intronstarts1[-1][x]:intronstarts1[-1][x]+(matches[x+1][1]-(matches[x][2]+matches[x][1]))]
                                              +'|mut|'+ test1[matches[x][2]+matches[x][1]:matches[x+1][1]])
                                    else: 
                                        muts.append(varsq[intronstarts1[-1][x]:intronstarts1[-1][x]+(matches[x+2][1]-(matches[x+1][2]+matches[x+1][1]))]
                                              +'|mut|'+ test1[matches[x+1][2]+matches[x+1][1]:matches[x+2][1]])
                            mutation1.append(muts)
                    test2=test_full.split('_')[1]
                    alal.set_seq2(test2)    
                    matches=alal.get_matching_blocks()
                    if (len(matches)==2)&(matches[0][0]+matches[0][2]>=200):
                        intronstarts2.append([0])
                        intronends2.append([0])
                        mutation2.append(['wt'])
                    elif len(matches)==1:
                        intronstarts2.append([-150])
                        intronends2.append([-1])
                        mutation2.append(['unknown'])
                    else: 
                        intronstarts2.append([matches[x][0]+matches[x][2] for x in range(len(matches)-2) if (matches[x+1][1]>9)&(matches[x][2]>9)])
                        intronends2.append([matches[x+1][0] - (matches[x+1][1]-(matches[x][2]+matches[x][1])) for x in range(len(matches)-2) if (matches[x+1][1]>9)&(matches[x][2]>9)])
                        if (intronstarts2[-1]==[])|(intronends2[-1]==[]):
                            mutation2.append(['wt'])
                        else:
                            muts2=[]
                            for x in range(len(intronstarts2[-1])):
                                if intronstarts2[-1][x]!=intronends2[-1][x]:
                                    muts2.append('del')
                                else:
                                    if (matches[x+1][1]>9)&(matches[x][2]>9):
                                        muts2.append(varsq[intronstarts2[-1][x]:intronstarts2[-1][x]+(matches[x+1][1]-(matches[x][2]+matches[x][1]))]
                                              +'|mut|'+ test2[matches[x][2]+matches[x][1]:matches[x+1][1]])
                                    else: 
                                        muts2.append(varsq[intronstarts2[-1][x]:intronstarts2[-1][x]+(matches[x+2][1]-(matches[x+1][2]+matches[x+1][1]))]
                                              +'|mut|'+ test2[matches[x+1][2]+matches[x+1][1]:matches[x+2][1]])
                            mutation2.append(muts2)

            for i in range(len(intronstarts1)):
                if ((intronstarts1[i]==[])&(intronstarts2[i]==[]))&((intronends1[i]==[])&(intronends2[i]==[])):
                    intronstarts.append(0)
                    intronends.append(0)
                    mutation.append([])
                    classifier.append('wt')
#                    mutation.append(mutation1[i]+mutation2[i])
#                    reads.append(rr.keys()[i])
#                    readcount.append(rr.values[i])
                elif ((intronstarts1[i]==[])&(intronstarts2[i]==[0]))&((intronends1[i]==[])&(intronends2[i]==[0])):
                    intronstarts.append(0)
                    intronends.append(0)
                    mutation.append([])
                    classifier.append('wt')
                else:
                    intronstarts.append([x-len(upstream) for x in intronstarts1[i]]+[x-len(upstream) for x in intronstarts2[i]])
                    intronends.append([x-len(upstream) for x in intronends1[i]]+[x-len(upstream) for x in intronends2[i]])
                    mutation.append(mutation1[i]+mutation2[i])
                    classifier.append('_'.join([str(intronstarts[i][x])+'|'+mutation[i][x]+'|'+str(intronends[i][x]) for x in range(len(intronstarts[i]))]))
                    if ('del' not in classifier[-1])&('mut' not in classifier[-1])&('wt' in classifier[-1]):
                        classifier[-1]='wt'
#                    reads.append(rr.keys()[i])
#                    readcount.append(rr.values[i])
                    
            
            if (len(intronstarts)>0)&(len(intronends)>0):
                introns=pd.DataFrame([intronstarts,intronends, mutation, reads,classifier]).transpose()
                introns.columns=['intronstart','intronend', 'mutation', 'read','classifier']
                introns['identifier']=introns.classifier.apply(lambda x: updateclassifier(x))
                introns['barcode']=introns.read.apply(lambda x: x[x.find(upstream)-8:x.find(upstream)])
                                
                for ident in introns.identifier.value_counts()[introns.identifier.value_counts()>thresholdreads].index:
                    intron=introns[introns.identifier==ident]
                    rna.loc[str(var)+'_'+ident,'libindex']=var
                    rna.loc[str(var)+'_'+ident,'variant']=ident
                    rna.loc[str(var)+'_'+ident,'read']=intron.read.values[0]
                    rna.loc[str(var)+'_'+ident,'numberreads']=len(intron)
                    if (kind=='bins')|(kind=='stop'):
                        rna.loc[str(var)+'_'+ident,'binbarcodes']=' '.join(intron.barcode)
                        rna.loc[str(var)+'_'+ident,'binbarcodecounts']=len(rna.loc[str(var)+'_'+ident,'binbarcodes'].split(' '))
                    
                        cov.loc[str(var)+'_'+ident]=0
                        for barcode in rna.loc[str(var)+'_'+ident,'binbarcodes'].split(' '): 
                            try:
                                cov.loc[str(var)+'_'+ident, binbc[barcode]]+=1
                            except:
                                pass        
    if (kind=='bins')|(kind=='stop'):
        return rna, cov
    else:
        return rna


def classify_variants(df, fslib, frame):
    def classify_variant(variant, varsq, frame):
        if variant=='wt':
            varclass='wt'
            gfpframe=frame
            gfpmadewithoutfs=False
            changeinnumberofstopsframe0=0
            changeinnumberofstopsframem1=0
            changeinnumberofstopsframep1=0
        elif '_' not in variant:
            if 'mut' in variant:
                startcoord=int(variant.split('|')[0])
                endcoord=int(variant.split('|')[-1])
                varsqmut=varsq[:startcoord]+variant.split('|')[-2]+varsq[endcoord+1:]
                varclass='pointmutation'
                gfpframe=frame
                gfpmadewithoutfs=False
                changeinnumberofstopsframe0=Seq(varsqmut).translate().count('*')-Seq(varsq).translate().count('*')
                changeinnumberofstopsframem1=Seq(varsqmut[2:]).translate().count('*')-Seq(varsq[2:]).translate().count('*')
                changeinnumberofstopsframep1=Seq(varsqmut[1:]).translate().count('*')-Seq(varsq[1:]).translate().count('*')
            elif 'del' in variant:
                startcoord=int(variant.split('|')[0])
                endcoord=int(variant.split('|')[-1])
                varsqmut=varsq[:startcoord]+varsq[endcoord:]
                
                varclass='deletion '+str(endcoord-startcoord)+'nt'
                gfpframe=(frame-endcoord+startcoord) % 3
                gfpmadewithoutfs=bool(Seq(varsqmut).translate().find('*')==-1)
                changeinnumberofstopsframe0=Seq(varsqmut).translate().count('*')-Seq(varsq).translate().count('*')
                changeinnumberofstopsframem1=Seq(varsqmut[2:]).translate().count('*')-Seq(varsq[2:]).translate().count('*')
                changeinnumberofstopsframep1=Seq(varsqmut[1:]).translate().count('*')-Seq(varsq[1:]).translate().count('*')
        else:
            varsqmut=''
            coordfromlast=0
            gfpframe_total=frame
            for varpart in variant.split('_'):
                varsqmut+=varsq[coordfromlast:int(varpart.split('|')[0])]
                coordfromlast=int(varpart.split('|')[-1])
                if 'mut' in varpart:
                    varsqmut+=varpart.split('|')[-2]
                    coordfromlast+=1
                if 'del' in varpart:
                    gfpframe_total=gfpframe_total-int(varpart.split('|')[-1])+int(varpart.split('|')[0])
                    
            varclass='multiple'
            gfpframe=gfpframe_total % 3
            gfpmadewithoutfs=bool(Seq(varsqmut).translate().find('*')==-1)
            changeinnumberofstopsframe0=Seq(varsqmut).translate().count('*')-Seq(varsq).translate().count('*')
            changeinnumberofstopsframem1=Seq(varsqmut[2:]).translate().count('*')-Seq(varsq[2:]).translate().count('*')
            changeinnumberofstopsframep1=Seq(varsqmut[1:]).translate().count('*')-Seq(varsq[1:]).translate().count('*')
        
        return pd.Series([varclass, 
                          gfpframe, 
                          gfpmadewithoutfs, 
                          changeinnumberofstopsframe0, 
                          changeinnumberofstopsframem1, 
                          changeinnumberofstopsframep1],
                index=['varclass','gfpframe','gfpmadewithoutfs','changeinnumberofstopsframe0','changeinnumberofstopsframem1','changeinnumberofstopsframep1'])
                
    return df.apply(lambda x: classify_variant(x.loc['variant'], fslib.varseq[int(x.loc['libindex'])], frame), axis=1)


def create_wt_combmuts_bincov(df, cov): 
    covwtmuts=pd.DataFrame(columns=cov.columns)  
    vardfwtmut=pd.DataFrame()  
    for var, group in df.groupby(by='libindex'):
        covwtmuts.loc[str(var)+'_wt']=cov.loc[group[(group.varclass=='wt')].index].sum()
        covwtmuts.loc[str(var)+'_del']=cov.loc[group[(group.varclass!='multiple')&(group.gfpframe==1)].index].sum()
        covwtmuts.loc[str(var)+'_ins']=cov.loc[group[(group.varclass!='multiple')&(group.gfpframe==0)&(group.gfpmadewithoutfs==True)].index].sum()
        covwtmuts.loc[str(var)+'_del_withmultiple']=cov.loc[group[(group.gfpframe==1)].index].sum()
        covwtmuts.loc[str(var)+'_ins_withmultiple']=cov.loc[group[(group.gfpframe==0)&(group.gfpmadewithoutfs==True)].index].sum()
        vardfwtmut.loc[str(var)+'_wt','libindex']=var
        vardfwtmut.loc[str(var)+'_wt','variant']='wt'
        vardfwtmut.loc[str(var)+'_wt','numberreads']=cov.loc[group[(group.varclass=='wt')].index].sum().sum()
        vardfwtmut.loc[str(var)+'_del','libindex']=var
        vardfwtmut.loc[str(var)+'_del','variant']='del'
        vardfwtmut.loc[str(var)+'_del','numberreads']=cov.loc[group[(group.varclass!='multiple')&(group.gfpframe==1)].index].sum().sum()
        vardfwtmut.loc[str(var)+'_ins','libindex']=var
        vardfwtmut.loc[str(var)+'_ins','variant']='ins'
        vardfwtmut.loc[str(var)+'_ins','numberreads']=cov.loc[group[(group.varclass!='multiple')&(group.gfpframe==0)&(group.gfpmadewithoutfs==True)].index].sum().sum()
        vardfwtmut.loc[str(var)+'_del_withmultiple','libindex']=var
        vardfwtmut.loc[str(var)+'_del_withmultiple','variant']='del_withmultiple'
        vardfwtmut.loc[str(var)+'_del_withmultiple','numberreads']=cov.loc[group[(group.gfpframe==1)].index].sum().sum()
        vardfwtmut.loc[str(var)+'_ins_withmultiple','libindex']=var
        vardfwtmut.loc[str(var)+'_ins_withmultiple','variant']='ins_withmultiple'
        vardfwtmut.loc[str(var)+'_ins_withmultiple','numberreads']=cov.loc[group[(group.gfpframe==0)&(group.gfpmadewithoutfs==True)].index].sum().sum()
    return covwtmuts, vardfwtmut



def create_wt_combmuts_bincov_p1(df, cov): 
    covwtmuts=pd.DataFrame(columns=cov.columns)  
    vardfwtmut=pd.DataFrame()  
    for var, group in df.groupby(by='libindex'):
        covwtmuts.loc[str(var)+'_wt']=cov.loc[group[(group.varclass=='wt')].index].sum()
        covwtmuts.loc[str(var)+'_ins']=cov.loc[group[(group.varclass!='multiple')&(group.gfpframe==2)].index].sum()
        covwtmuts.loc[str(var)+'_del']=cov.loc[group[(group.varclass!='multiple')&(group.gfpframe==0)&(group.gfpmadewithoutfs==True)].index].sum()
        covwtmuts.loc[str(var)+'_ins_withmultiple']=cov.loc[group[(group.gfpframe==2)].index].sum()
        covwtmuts.loc[str(var)+'_del_withmultiple']=cov.loc[group[(group.gfpframe==0)&(group.gfpmadewithoutfs==True)].index].sum()
        vardfwtmut.loc[str(var)+'_wt','libindex']=var
        vardfwtmut.loc[str(var)+'_wt','variant']='wt'
        vardfwtmut.loc[str(var)+'_wt','numberreads']=cov.loc[group[(group.varclass=='wt')].index].sum().sum()
        vardfwtmut.loc[str(var)+'_ins','libindex']=var
        vardfwtmut.loc[str(var)+'_ins','variant']='del'
        vardfwtmut.loc[str(var)+'_ins','numberreads']=cov.loc[group[(group.varclass!='multiple')&(group.gfpframe==2)].index].sum().sum()
        vardfwtmut.loc[str(var)+'_del','libindex']=var
        vardfwtmut.loc[str(var)+'_del','variant']='ins'
        vardfwtmut.loc[str(var)+'_del','numberreads']=cov.loc[group[(group.varclass!='multiple')&(group.gfpframe==0)&(group.gfpmadewithoutfs==True)].index].sum().sum()
        vardfwtmut.loc[str(var)+'_ins_withmultiple','libindex']=var
        vardfwtmut.loc[str(var)+'_ins_withmultiple','variant']='del_withmultiple'
        vardfwtmut.loc[str(var)+'_ins_withmultiple','numberreads']=cov.loc[group[(group.gfpframe==2)].index].sum().sum()
        vardfwtmut.loc[str(var)+'_del_withmultiple','libindex']=var
        vardfwtmut.loc[str(var)+'_del_withmultiple','variant']='ins_withmultiple'
        vardfwtmut.loc[str(var)+'_del_withmultiple','numberreads']=cov.loc[group[(group.gfpframe==0)&(group.gfpmadewithoutfs==True)].index].sum().sum()
    return covwtmuts, vardfwtmut
    
def check_binprimercorrelation(data):
    from scipy.stats import pearsonr
    colnames=[]
    for i in range(1,17):
        for j in range(1,4):
            colnames.append('Bin#'+str(i)+'_'+str(j)+'_Fw')
    print '1 vs 2'
    for i in range(1,17):
        print str(i) + '\t' + str(pearsonr(data['Bin#'+str(i)+'_1_Fw'],data['Bin#'+str(i)+'_2_Fw']))
    print 
    print '1 vs 3'        
    for i in range(1,17):
        print str(i) + '\t' + str(pearsonr(data['Bin#'+str(i)+'_1_Fw'],data['Bin#'+str(i)+'_3_Fw']))
    
    
def combine_primers(data, to_combine=[1,2,3]):
    df=pd.DataFrame(0, index=data.index, columns=['Bin'+str(x) for x in range(1,17)])
    for i in range(1,17):
        for j in to_combine:
            df['Bin'+str(i)]+=data['Bin#'+str(i)+'_'+str(j)+'_Fw']
    return df

def FilterMinBinCoverage(data, threshold):
    df=data.copy()
    df[df < threshold] = 0
    return df


def FilterPercentPerBin(data, threshold):
    df=data.copy()
    for idx in df.index:
        if df.loc[idx].sum()>0:
            for col in df.columns:
                if (np.float(df.loc[idx,col])/df.loc[idx].sum()*100 < threshold):
                    df.loc[idx,col]=0
    return df       


def FilterExtremeBins(data):
    df=data.copy()
    def checkrow(row):
        if row.sum()>0:
            if (row.iloc[0]>row.iloc[1]):
                row.iloc[0]=row.iloc[1]
            if (row.iloc[-1]>row.iloc[-2]):
                row.iloc[-1]=row.iloc[-2]
        return row
    df.apply(lambda x: checkrow(x), axis=1)
    return df
 
    
def FilterIsolated(data):
    df=data.copy()
    def checkrow(row):
        if row.sum()>0:
            if (row.iloc[1]==0):
                row.iloc[0]=0
            for col in range(1,15):
                if (row.iloc[col-1]==0)&(row.iloc[col+1]==0):
                    row.iloc[col]=0
            if (row.iloc[14]==0):
                row.iloc[15]=0
        return row
    df.apply(lambda x: checkrow(x), axis=1)
    return df

def FilterIsolated12bins(data):
    df=data.copy()
    def checkrow(row):
        if row.sum()>0:
            if (row.iloc[1]==0):
                row.iloc[0]=0
            for col in range(1,11):
                if (row.iloc[col-1]==0)&(row.iloc[col+1]==0):
                    row.iloc[col]=0
            if (row.iloc[10]==0):
                row.iloc[11]=0
        return row
    df.apply(lambda x: checkrow(x), axis=1)
    return df


def FilterPercentKept(original, filtered, threshold):
    df=filtered.copy()
    occur = 0
    for var, group in df.groupby(level=0):
        if (group.summedreads.sum()>0):
            if (group.summedreads.sum()/original.loc[var].summedreads.sum()*100 <threshold):
                df.loc[var]=0
                occur +=1
    print(occur)
    return df    
        

def NormalizeByBinPercentages(data, correctionfactor):
    df=data.copy()
    df=df.apply(lambda x: x*correctionfactor, axis=1)
    return df           

def NormalizeToReadNumbers(data):
    df=data.copy()
    def checkrow(row):
        rs=float(row.sum())
        if rs>0:
            for binn in range(0,16):
                row.iloc[binn] = float(row.iloc[binn])/rs
        return row
    df=df.apply(lambda x: checkrow(x), axis=1)
    return df

def NormalizeToReadNumbers12bins(data):
    df=data.copy()
    def checkrow(row):
        rs=float(row.sum())
        if rs>0:
            for binn in range(0,12):
                row.iloc[binn] = float(row.iloc[binn])/rs
        return row
    df=df.apply(lambda x: checkrow(x), axis=1)
    return df
        
def ExpCalc(df, estimate):
    from statsmodels.stats.weightstats import DescrStatsW
    expr=pd.DataFrame()
    def checkrow(row):
        st=DescrStatsW(estimate, weights=row)
        return pd.Series([st.mean, st.std, st.std_mean])
    expr=df.apply(lambda x: checkrow(x), axis=1)
    expr.columns=['wav','wstd','wstdm']
    return expr

def findpeaks(data, xval):
    from peakdetect import peakdet
    from savitzky_golay import savitzky_golay
    df=data.copy()
    def checkrow(row):
        rawpeaks=len(peakdet(row.values, 0.02, xval.values)[0])
        if (row.iloc[-1]>=row.iloc[-2])&(row.iloc[-1]>0):
            rawpeaks+=1
        smoothed = savitzky_golay(row.values, 3, 1)
        smoothedpeaks=len(peakdet(smoothed, 0.02, xval.values)[0])
        if (smoothed[-1]>=smoothed[-2])&(smoothed[-1]>0):
            smoothedpeaks+=1
        return pd.Series([rawpeaks, smoothedpeaks])
    peaks=df.apply(lambda x: checkrow(x), axis=1)
    peaks.columns=['rawpeaks','smoothedpeaks']
    return peaks

def smooth(data):
    from savitzky_golay import savitzky_golay
    df=data.copy()
    def checkrow(row):
        smoothed = savitzky_golay(row.values, 3, 1)
        return pd.Series(smoothed)
    dfsmoothed=df.apply(lambda x: checkrow(x), axis=1)
    dfsmoothed.columns=df.columns
    return dfsmoothed