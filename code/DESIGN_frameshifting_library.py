# -*- coding: utf-8 -*-

import pandas as pd

from Bio.Seq import Seq

from Bio import SeqIO

import random

import RNA


from Bio import SeqUtils

SynonymousCodons = { 
    'CYS': ['TGT', 'TGC'], 
    'ASP': ['GAT', 'GAC'], 
    'SER': ['TCT', 'TCG', 'TCA', 'TCC', 'AGC', 'AGT'], 
    'GLN': ['CAA', 'CAG'], 
    'MET': ['ATG'], 
    'ASN': ['AAC', 'AAT'], 
    'PRO': ['CCT', 'CCG', 'CCA', 'CCC'], 
    'LYS': ['AAG', 'AAA'], 
    'STOP': ['TAG', 'TGA', 'TAA'], 
    'THR': ['ACC', 'ACA', 'ACG', 'ACT'], 
    'PHE': ['TTT', 'TTC'], 
    'ALA': ['GCA', 'GCC', 'GCG', 'GCT'], 
    'GLY': ['GGT', 'GGG', 'GGA', 'GGC'], 
    'ILE': ['ATC', 'ATA', 'ATT'], 
    'LEU': ['TTA', 'TTG', 'CTC', 'CTT', 'CTG', 'CTA'], 
    'HIS': ['CAT', 'CAC'], 
    'ARG': ['CGA', 'CGC', 'CGG', 'CGT', 'AGG', 'AGA'], 
    'TRP': ['TGG'], 
    'VAL': ['GTA', 'GTC', 'GTG', 'GTT'], 
    'GLU': ['GAG', 'GAA'], 
    'TYR': ['TAT', 'TAC'] 
    } 

gencode={}
for i in SynonymousCodons.keys():
    gencode[SeqUtils.seq1(i)]=SynonymousCodons[i]

gencode['*']=['TAG', 'TGA', 'TAA']

codonmaxgc={
 'A': ['GCG'],
 'C': ['TGC'],
 'D': ['GAC'],
 'E': ['GAG'],
 'F': ['TTC'],
 'G': ['GGC'],
 'H': ['CAC'],
 'I': ['ATC'],
 'K': ['AAG'],
 'L': ['CTG'],
 'M': ['ATG'],
 'N': ['AAC'],
 'P': ['CCG'],
 'Q': ['CAG'],
 'R': ['CGG'],
 'S': ['TCC'],
 'T': ['ACG'],
 'V': ['GTC'],
 'W': ['TGG'],
 '*': ['TGA'],
 'Y': ['TAC']}

codonmingc={
 'A': ['GCA'],
 'C': ['TGT'],
 'D': ['GAT'],
 'E': ['GAA'],
 'F': ['TTT'],
 'G': ['GGA'],
 'H': ['CAT'],
 'I': ['ATA'],
 'K': ['AAA'],
 'L': ['TTA'],
 'M': ['ATG'],
 'N': ['AAT'],
 'P': ['CCT'],
 'Q': ['CAA'],
 'R': ['AGA'],
 'S': ['AGT'],
 'T': ['ACT'],
 'V': ['GTA'],
 'W': ['TGG'],
 '*': ['TAA'],
 'Y': ['TAT']}
#%%


#%% Library design - starting from previously reported PRF sites (fs, selected ones for systematic sequence alterations: fsmain)

fs=pd.read_csv('./design/fs.csv',index_col='FSDB ID')
fsmain=pd.read_csv('./design/fsmain.csv',index_col='name')


f=open('./design/varseqplus70nteachside.txt','w')
for i in fsmain.index:
    region_start = int(fsmain.loc[i, 'slippery site position'] + len(fsmain.loc[i, 'slippery']) - 112)
    region_end = int(fsmain.loc[i, 'slippery site position'] + len(fsmain.loc[i, 'slippery']) + 190)
    out=str('>' + i + '\n' + fsmain.loc[i, 'full length sequence'][region_start:region_end] + '\n')
    f.write(out)

f.close()
#%% exchange stop codons with TGG trp

fs_mutatedstops=pd.DataFrame()

for i in fs.index:
    fs_mutatedstops.loc[i,'organism']=fs.loc[i,'organism']
    fs_mutatedstops.loc[i,'gene']=fs.loc[i,'gene']
    fs_mutatedstops.loc[i,'sequence ID']=fs.loc[i,'sequence ID']    
    my_seq=Seq(fs.loc[i,'varseq']).tomutable()
    for pos in range(39, 160):
        if (my_seq[pos:pos+3]=='TAG') or (my_seq[pos:pos+3]=='TGA') or (my_seq[pos:pos+3]=='TAA'):
           my_seq[pos:pos+3]='TGG'
           
    fs_mutatedstops.loc[i,'varseq']=str(my_seq).upper();
    
fs_mutatedstops.to_csv('./design/fs_mutatedstops.csv')

for i in fsmain.index:
    fsmain.loc[i,'varseqstopmut']=fs_mutatedstops.varseq[fs_mutatedstops.loc[:,'sequence ID']==fsmain.loc[i,'sequence ID']].values[0]

'''
f=open(dropboxprefix + 'Dropbox/WIS/Library/fsmain_endogenous_and_stopmut_varseq.txt','w')
for i in fsmain.index:
    f.write(str(fsmain.varseq[i])+'\n' + str(fsmain.varseqstopmut[i])+'\n')

f.close()

f=open(dropboxprefix + 'Dropbox/WIS/Library/fsmain_endogenous_and_stopmut_varseq_afterslip.txt','w')
for i in fsmain.index:
    f.write('>'+str(i) + '\n' + str(fsmain.varseq[i])+'\n' + '>'+str(i) + 'stop\n' + str(fsmain.varseqstopmut[i])+'\n')

f.close()
'''

#%% permutations

fsmainnotm1=fsmain[(fsmain.type!=-1)]
fsmainm1perm=fsmain[(fsmain.type==-1)&(fsmain.stopsin0==0)]
fsmainm1stopperm=fsmain[(fsmain.type==-1)&(fsmain.stopsin0!=0)]
fsmainmixedperm=fsmain.loc[['hiv','oaz','influenza','nsp2f','herpes','srvl','sars']]
fsmainmixedpermnonstop=fsmainmixedperm.rename(columns={'varseq':'varseq_original','varseqstopmut':'varseq'})

def make_fs_perm(df,name):
    dfperm=pd.DataFrame()
    l=0
    for i in df.index:
        for j in df.index:
            for k in df.index:
                dfperm.loc[l,'varseq']=df.varseq[i][:30] + df.varseq[j][30:42] + df.varseq[k][42:]
                dfperm.loc[l,'upstream']=i
                dfperm.loc[l,'slippery']=j
                dfperm.loc[l,'downstream']=k
                dfperm.loc[l,'group']=str(name)
                l=l+1
                
    return dfperm
    
fsmainnotm1var=make_fs_perm(fsmainnotm1,'not m1')
fsmainm1permvar=make_fs_perm(fsmainm1perm, 'all m1 without stop in frame0')
fsmainm1stoppermvar=make_fs_perm(fsmainm1stopperm, 'all 3 m1 with stop in frame 0')
fsmainmixedpermvar=make_fs_perm(fsmainmixedperm, 'combination of types, endogenous')
fsmainmixedpermnonstopvar=make_fs_perm(fsmainmixedpermnonstop, 'combination of types, stop codons mutated')

fs_permutations=pd.concat([fsmainnotm1var,fsmainm1permvar,fsmainm1stoppermvar,fsmainmixedpermvar,fsmainmixedpermnonstopvar],ignore_index=True)
fs_permutations.to_csv('./design/fs_permutations.csv')

#%% SLIPPERY SITE

# create slippery sites
l=0
slipsites=pd.Series()
for i in ['A','C','G','T']:
    for j in ['A','C','G','T']:
        for k in ['A','C','G','T']:
            slipsites.loc[l]=Seq(i+i+i+j+j+j+k)
            l=l+1

fsmainm1=fsmain[(fsmain.type==-1)]

fs_slipperyvar=pd.DataFrame(columns=list(fsmain.columns))
k=0
for i in fsmainm1.index:
    for j in slipsites:
        fs_slipperyvar=fs_slipperyvar.append(fsmainm1.loc[i],ignore_index=True)
        fs_slipperyvar.varseq[k]=fsmainm1.varseq[i][:35] + j + fsmainm1.varseq[i][42:]
        fs_slipperyvar.slippery[k]=str(j)
        fs_slipperyvar.loc[k,'group']='subgroup 1'
        k=k+1
    if (fsmainm1.stopsin0[i]>0):
        for j in slipsites:
            fs_slipperyvar=fs_slipperyvar.append(fsmainm1.loc[i],ignore_index=True)
            fs_slipperyvar.varseq[k]=fsmainm1.varseq[i][:35] + j + fsmainm1.varseqstopmut[i][42:]
            fs_slipperyvar.slippery[k]=str(j)
            fs_slipperyvar.loc[k,'group']='subgroup 1 stops mutated'
            k=k+1



SynonymousCodons = { 
    'CYS': ['TGT', 'TGC'], 
    'ASP': ['GAT', 'GAC'], 
    'SER': ['TCT', 'TCG', 'TCA', 'TCC', 'AGC', 'AGT'], 
    'GLN': ['CAA', 'CAG'], 
    'MET': ['ATG'], 
    'ASN': ['AAC', 'AAT'], 
    'PRO': ['CCT', 'CCG', 'CCA', 'CCC'], 
    'LYS': ['AAG', 'AAA'], 
    'STOP': ['TAG', 'TGA', 'TAA'], 
    'THR': ['ACC', 'ACA', 'ACG', 'ACT'], 
    'PHE': ['TTT', 'TTC'], 
    'ALA': ['GCA', 'GCC', 'GCG', 'GCT'], 
    'GLY': ['GGT', 'GGG', 'GGA', 'GGC'], 
    'ILE': ['ATC', 'ATA', 'ATT'], 
    'LEU': ['TTA', 'TTG', 'CTC', 'CTT', 'CTG', 'CTA'], 
    'HIS': ['CAT', 'CAC'], 
    'ARG': ['CGA', 'CGC', 'CGG', 'CGT', 'AGG', 'AGA'], 
    'TRP': ['TGG'], 
    'VAL': ['GTA', 'GTC', 'GTG', 'GTT'], 
    'GLU': ['GAG', 'GAA'], 
    'TYR': ['TAT', 'TAC'] 
    } 

codons=[]
for i in SynonymousCodons.keys():
    codons=codons + SynonymousCodons[i]
    
for i in codons:
    fs_slipperyvar=fs_slipperyvar.append(fsmain.loc['oaz'],ignore_index=True)
    fs_slipperyvar.varseq[k]=fsmain.varseq['oaz'][:36] + Seq(i) + fsmain.varseq['oaz'][39:]
    fs_slipperyvar.slippery[k]=str(i) + str(fsmain.varseq['oaz'][39:42])
    fs_slipperyvar.loc[k,'group']='subgroup 2 - oaz'
    k=k+1
    
    fs_slipperyvar=fs_slipperyvar.append(fsmain.loc['oaz'],ignore_index=True)
    fs_slipperyvar.varseq[k]=fsmain.varseq['oaz'][:39] + Seq(i) + fsmain.varseq['oaz'][42:]
    fs_slipperyvar.slippery[k]=str(fsmain.varseq['oaz'][36:39]) + str(i)
    fs_slipperyvar.loc[k,'group']='subgroup 2 - oaz stops mutated'
    k=k+1
    
    fs_slipperyvar=fs_slipperyvar.append(fsmain.loc['oaz'],ignore_index=True)
    fs_slipperyvar.varseq[k]=fsmain.varseqstopmut['oaz'][:36] + Seq(i) + fsmain.varseqstopmut['oaz'][39:]
    fs_slipperyvar.slippery[k]=str(i) + str(fsmain.varseq['oaz'][39:42])
    fs_slipperyvar.loc[k,'group']='subgroup 2 - oaz'
    k=k+1
    
    fs_slipperyvar=fs_slipperyvar.append(fsmain.loc['oaz'],ignore_index=True)
    fs_slipperyvar.varseq[k]=fsmain.varseqstopmut['oaz'][:39] + Seq(i) + fsmain.varseqstopmut['oaz'][42:]
    fs_slipperyvar.slippery[k]=str(fsmain.varseq['oaz'][36:39]) + str(i)
    fs_slipperyvar.loc[k,'group']='subgroup 2 - oaz stops mutated'
    k=k+1
    
fsm=fsmain
#fsm.drop(['roadblock','paper','importance','FSsequence','FSseqlength','full length sequence','slippery site position','0stopbeforefs','0stopafterfs','m1stopafterfs','p1stopafterfs','structfullWT','MFEfullWT','structafterslipWT','MFEafterslipWT','structfullSTOPMUT','MFEfullSTOPMUT','structafterslipSTOPMUT','MFEafterslipSTOPMUT'],axis=1,inplace=True)

fsslip3=fsm[fsm.stopsinp1==0]

# create slippery sites
l=0
slipsites=pd.Series()
for i in ['A','C','G','T']:
    if (i=='T'):
        for j in ['C','G','T']:
            slipsites.loc[l]=Seq(i+i+j+j+j+j+j)
            l=l+1
    else:    
        for j in ['A','C','G','T']:
            slipsites.loc[l]=Seq(i+i+j+j+j+j+j)
            l=l+1

for i in fsslip3.index:
    for j in slipsites:
        fs_slipperyvar=fs_slipperyvar.append(fsslip3.loc[i],ignore_index=True)
        fs_slipperyvar.varseq[k]=fsslip3.varseq[i][:35] + j + fsslip3.varseq[i][42:]
        fs_slipperyvar.slippery[k]=str(j)
        fs_slipperyvar.loc[k,'group']='subgroup 3'
        k=k+1


l=0
slipsites=pd.Series()
for i in ['A','C','G','T']:
    for j in ['A','C','G','T']:
        for m in ['A','C','G','T']:
            slipsites.loc[l]=Seq(i+j+j+m+m+m+m)
            l=l+1

fsslip4=fsm.loc[['wnv','herpes','nsp2f','plrv','hiv','htcl']]

for i in fsslip4.index:
    for j in slipsites:
        fs_slipperyvar=fs_slipperyvar.append(fsslip4.loc[i],ignore_index=True)
        fs_slipperyvar.varseq[k]=fsslip4.varseq[i][:35] + j + fsslip4.varseq[i][42:]
        fs_slipperyvar.slippery[k]=str(j)
        fs_slipperyvar.loc[k,'group']='subgroup 4'
        k=k+1

    if (fsslip4.stopsinm1[i]>0):
        for j in slipsites:
            fs_slipperyvar=fs_slipperyvar.append(fsslip4.loc[i],ignore_index=True)
            fs_slipperyvar.varseq[k]=fsslip4.varseqstopmut[i][:35] + j + fsslip4.varseqstopmut[i][42:]
            fs_slipperyvar.slippery[k]=str(j)
            fs_slipperyvar.loc[k,'group']='subgroup 4 stopmut'
            k=k+1
    

fsslip5=fs[fs.type==0]

for i in fsslip5.index:
    fs_slipperyvar=fs_slipperyvar.append(fsmain.loc['influenza'],ignore_index=True)
    fs_slipperyvar.varseq[k]=fsslip5.varseq[i][:30] + fsmain.varseq['influenza'][30:]
    fs_slipperyvar.loc[k,'group']='subgroup 5 influenza upstream from ' + i
    k=k+1

    fs_slipperyvar=fs_slipperyvar.append(fsmain.loc['influenza'],ignore_index=True)
    fs_slipperyvar.varseq[k]=fsmain.varseq['influenza'][:30] + fsslip5.varseq[i][30:42] + fsmain.varseq['influenza'][42:]
    fs_slipperyvar.loc[k,'group']='subgroup 5 influenza slippery from ' + i
    k=k+1

    fs_slipperyvar=fs_slipperyvar.append(fsmain.loc['influenza'],ignore_index=True)
    fs_slipperyvar.varseq[k]=fsmain.varseq['influenza'][:42] + fsslip5.varseq[i][42:]
    fs_slipperyvar.loc[k,'group']='subgroup 5 influenza downstream from ' + i
    k=k+1

infslip=['TCCAAACGT','TCCCCCCGT','TCCGGGCGT','TCCTTTCGA','TCCTCTCGT',
         'TCCTTTTGT','TCCTTTGCT','TCCTTTCGG','TCCTTTCGC','TCCTATCGT']

for i in infslip:
    fs_slipperyvar=fs_slipperyvar.append(fsmain.loc['influenza'],ignore_index=True)
    fs_slipperyvar.varseq[k]=fsmain.varseq['influenza'][:33] + Seq(i) + fsmain.varseq['influenza'][42:]
    fs_slipperyvar.loc[k,'group']='subgroup 5 influenza slippery site mutations, codon 2 and 3'
    k=k+1
             
    fs_slipperyvar=fs_slipperyvar.append(fsmain.loc['influenza'],ignore_index=True)
    fs_slipperyvar.varseq[k]=fsmain.varseqstopmut['influenza'][:33] + Seq(i) + fsmain.varseqstopmut['influenza'][42:]
    fs_slipperyvar.loc[k,'group']='subgroup 5 influenza stop mutant slippery site mutations, codon 2 and 3'
    k=k+1
             
codonsnostop=['TGT', 'TGC', 'GAT', 'GAC', 'TCT', 'TCG', 'TCA', 'TCC', 'AGC', 'AGT', 'CAA',
 'CAG', 'AAG', 'AAA', 'GGT', 'GGG', 'GGA', 'GGC', 'CCT', 'CCG', 'CCA',
 'CCC', 'ACC', 'ACA', 'ACG', 'ACT', 'TTT', 'TTC', 'GCA',
 'GCC', 'GCG', 'GCT', 'CAT', 'CAC', 'ATG', 'ATC', 'ATA',
 'ATT', 'GAG', 'GAA', 'TTA', 'TTG', 'CTC', 'CTT', 'CTG',
 'CTA', 'CGA', 'CGC', 'CGG', 'CGT', 'AGG', 'AGA', 'TGG',
 'GTA', 'GTC', 'GTG', 'GTT', 'AAC', 'AAT', 'TAT', 'TAC']

for i in codonsnostop:
    fs_slipperyvar=fs_slipperyvar.append(fsmain.loc['influenza'],ignore_index=True)
    fs_slipperyvar.varseq[k]=fsmain.varseq['influenza'][:33] + Seq(i) + fsmain.varseq['influenza'][36:]
    fs_slipperyvar.loc[k,'group']='subgroup 5 influenza slippery site mutations, codon 1'
    k=k+1
             
    fs_slipperyvar=fs_slipperyvar.append(fsmain.loc['influenza'],ignore_index=True)
    fs_slipperyvar.varseq[k]=fsmain.varseqstopmut['influenza'][:33] + Seq(i) + fsmain.varseqstopmut['influenza'][36:]
    fs_slipperyvar.loc[k,'group']='subgroup 5 influenza stop mutant slippery site mutations, codon 1'
    k=k+1
    
fsslip6=fsmain[fsmain.type==-1]

for i in fsslip6.index:
    for j in range(0,6,1):
        for l in ['A','C','T','G']:
            seqtotest=fsslip6.varseq[i][:35+j] + Seq(l) + fsslip6.varseq[i][36+j:]
            if (seqtotest[:42].translate().find("*")==-1):
                fs_slipperyvar=fs_slipperyvar.append(fsslip6.loc[i],ignore_index=True)
                fs_slipperyvar.varseq[k]=seqtotest
                fs_slipperyvar.loc[k,'group']='subgroup 6, ' + l + ' at position ' + str(j+1) + ' in varseq'
                k=k+1
                


fs_slipperyvar.to_csv('./design/fs_slipperyvar.csv')


#%% SPACER

codonssafe=[]
for i in codonsnostop:
    if ((Seq('T' + i).translate().find('*')==-1)&(Seq('TA' + i).translate().find('*')==-1)):
        codonssafe.append(i)


for i in codonssafe:
    if (i[1:]=='TA')|(i[1:]=='TG')|(i[2:]=='T'):
        codonssafe.remove(i)

# create random spacer

randomnostop=pd.Series()

for i in range(0,10,1):
    randomnostop.loc[i]= random.choice(codonssafe) + random.choice(codonssafe) + random.choice(codonssafe) +random.choice(codonssafe)

'''
CGGCGGCAGCAA
TGCTACTCATGG
CGATACCGCTTC
TCACTCTGCCGG
TCCCAACGCCCG
'''

fs_spacervar = pd.DataFrame(columns=list(fsmain.columns))

fsspacer1=fsmain.loc[['nsp2f','plrv','hiv','herpes','srvl','htcl']]
m=0
for i in fsspacer1.index:
    for j in range(1,6):
        for k in range(1,13):
            fs_spacervar=fs_spacervar.append(fsmain.loc[i],ignore_index=True)
            if (fsmain.stopsin0[i]==0)&(fsmain.stopsinm1[i]==0)&(fsmain.stopsinp1[i]==0):
                fs_spacervar.varseq[m]=fsmain.varseq[i][:42] + Seq(randomnostop[j][:k]) + fsmain.varseq[i][42:-k]
                fs_spacervar.loc[m,'changes']='subgroup 1: ' + str(k) + 'nt from sequence ' + randomnostop[j] + ' inserted after slip'
            else:
                fs_spacervar.varseq[m]=fsmain.varseqstopmut[i][:42] + Seq(randomnostop[j][:k]) + fsmain.varseqstopmut[i][42:-k]
                fs_spacervar.loc[m,'changes']='subgroup 1: ' + str(k) + 'nt from sequence ' + randomnostop[j] + ' inserted after slip, stopmutant'
            m=m+1

for i in fsspacer1.index:
    for j in range(1,4):
        for k in range(1,7):
            fs_spacervar=fs_spacervar.append(fsmain.loc[i],ignore_index=True)
            if (fsmain.stopsin0[i]==0)&(fsmain.stopsinm1[i]==0)&(fsmain.stopsinp1[i]==0):
                fs_spacervar.varseq[m]=fsmain.varseq[i][:42] + fsmain.varseq[i][42+k:] + Seq(randomnostop[j][:k]) 
                fs_spacervar.loc[m,'changes']='subgroup 1: ' + str(k) + 'nts deleted after slip'
            else:
                fs_spacervar.varseq[m]=fsmain.varseqstopmut[i][:42] + fsmain.varseqstopmut[i][42+k:] + Seq(randomnostop[j][:k]) 
                fs_spacervar.loc[m,'changes']='subgroup 1: ' + str(k) + 'nts deleted after slip, stopmutant'
            m=m+1





fsspacer2=fsmain.loc[['oaz','hcv','influenza','wnv','siv','peg10','rous','sars','hcv','herv']]

for i in fsspacer2.index:
    for j in range(1,6):
        for k in range(1,5):
            fs_spacervar=fs_spacervar.append(fsmain.loc[i],ignore_index=True)
            fs_spacervar.varseq[m]=fsmain.varseq[i][:42] + Seq(randomnostop[j][:k*3]) + fsmain.varseq[i][42:-k*3]
            fs_spacervar.loc[m,'changes']='subgroup 2: ' + str(k*3) + 'nt from sequence ' + randomnostop[j] + ' inserted after slip'
            m=m+1

for i in fsspacer1.index:
    for j in range(1,4):
        for k in range(1,3):
            fs_spacervar=fs_spacervar.append(fsmain.loc[i],ignore_index=True)
            fs_spacervar.varseq[m]=fsmain.varseq[i][:42] + fsmain.varseq[i][42+k*3:] + Seq(randomnostop[j][:k*3]) 
            fs_spacervar.loc[m,'changes']='subgroup 2: ' + str(k*3) + 'nts deleted after slip'
            m=m+1


for i in fsmain.index:
    print (i)
    for k in ['A','C','G','T']:
        for l in range(0,6,1):
            seqtotest=fsmain.varseq[i][:42+j] + Seq(k) + fsmain.varseq[i][43+j:]
            if (seqtotest[42:48].translate().find("*")==-1)&(seqtotest[40:49].translate().find("*")==-1)&(seqtotest[41:50].translate().find("*")==-1):
                fs_spacervar=fs_spacervar.append(fsmain.loc[i],ignore_index=True)
                fs_spacervar.varseq[m]=seqtotest
                fs_spacervar.loc[m,'group']='subgroup 3: ' + k + ' at position ' + str(l+1) + ' after varseq'
                m=m+1

for i in fsmain.index:
    for k in codonssafe:
        fs_spacervar=fs_spacervar.append(fsmain.loc[i],ignore_index=True)
        fs_spacervar.varseq[m]=fsmain.varseq[i][:42] + Seq(k) + fsmain.varseq[i][45:]
        fs_spacervar.loc[m,'changes']='subgroup 3: Codon ' + k + ' to replace first codon after slip'
        m=m+1


fs_spacervar.to_csv('./design/fs_spacervar.csv')


#%% UPSTREAM REGION

fs_upstreamvar=pd.DataFrame(columns=list(fsmain.columns))
m=0
for i in fsmain.index:
    print (i)
    for j in range(1,3,1):
        
        for k in codonsnostop:
            fs_upstreamvar=fs_upstreamvar.append(fsmain.loc[i],ignore_index=True)
            fs_upstreamvar.varseq[m]=fsmain.varseq[i][:36-j*3] + Seq(k) + fsmain.varseq[i][36-(j-1)*3:]
            fs_upstreamvar.loc[m,'changes']='subgroup 1: Codon Number ' + str(j) + ' upstream of slipsite changed to ' + k 
            m=m+1



for i in fsmain.index:
    aaseq=Seq(fsmain.varseq[i][:30]).translate()
    for count in range (0,6,1):
        fs_upstreamvar=fs_upstreamvar.append(fsmain.loc[i],ignore_index=True)
        recod=str()
        for j in aaseq:
            recod=recod + random.choice(gencode[j])
        fs_upstreamvar.varseq[m]=Seq(recod) + fsmain.varseq[i][30:]
        fs_upstreamvar.loc[m,'changes']='subgroup 2: upstream sequence recoded randomly'
        m=m+1


for i in fsmain.index:
    [fsmain.loc[i,'structupstreamWT'],fsmain.loc[i,'MFEupstreamWT']]=RNA.fold(str(fsmain.loc[i,'varseq'][:42]))

import subprocess
import RNA


for i in fsmain.index:
    numberreturned=10
    inputstring=fsmain.loc[i,'structupstreamWT']
    slippery=str(fsmain.loc[i,'slippery'])
    gccontent=str('0.5')
    namevar=str(i + "_vienna_antarna")
    if len(slippery)==7:
        startseq=str('NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN' + slippery)
    if len(slippery)==9:
        startseq=str('NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN' + slippery)
    if len(slippery)==6:
        startseq=str('NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN' + slippery)

    print (namevar)
    
    commandlinestring=str("python ./antaRNA.py -Cseq '" + startseq + "' -Cstr '" + inputstring + "' -tGC 0.4 -tGCmax 0.8")
    
    currentcount=0
  
    while (currentcount < numberreturned):
        outputstring=subprocess.check_output(commandlinestring, shell=True).splitlines()
        seqtotest=Seq(outputstring[1]).back_transcribe()
        if (seqtotest[:36].translate().find("*")==-1):
            fs_upstreamvar=fs_upstreamvar.append(fsmain.loc[i],ignore_index=True)
            fs_upstreamvar.varseq[m]=seqtotest + fsmain.varseq[i][42:]
            fs_upstreamvar.loc[m,'changes']='subgroup 3: upstream antarna inverse fold'
            m=m+1
            currentcount=currentcount+1
            print(currentcount)
        else:
            print('-')

nostruct=[]
numberreturned=6
inputstring='..............................'

commandlinestring=str("python */antaRNA.py -Cstr '" + inputstring + "' -tGC 0.4 -tGCmax 0.8")

currentcount=0
  
while (currentcount < numberreturned):
    outputstring=subprocess.check_output(commandlinestring, shell=True).splitlines()
    seqtotest=Seq(outputstring[1]).back_transcribe()
    if (seqtotest.translate().find("*")==-1):
        nostruct.append(str(seqtotest))
        currentcount=currentcount+1
        print(currentcount)
    else:
        print('-')

hairpin3pr=[]
numberreturned=6
inputstring='.........(((((((((...)))))))))'

commandlinestring=str("python */antaRNA.py -Cstr '" + inputstring + "' -tGC 0.4 -tGCmax 0.8")

currentcount=0
  
while (currentcount < numberreturned):
    outputstring=subprocess.check_output(commandlinestring, shell=True).splitlines()
    seqtotest=Seq(outputstring[1]).back_transcribe()
    if (seqtotest.translate().find("*")==-1):
        hairpin3pr.append(str(seqtotest))
        currentcount=currentcount+1
        print(currentcount)
    else:
        print('-')

hairpin5pr=[]
numberreturned=6
inputstring='(((((((((...))))))))).........'

commandlinestring=str("python */antaRNA.py -Cstr '" + inputstring + "' -tGC 0.4 -tGCmax 0.8")

currentcount=0
  
while (currentcount < numberreturned):
    outputstring=subprocess.check_output(commandlinestring, shell=True).splitlines()
    seqtotest=Seq(outputstring[1]).back_transcribe()
    if (seqtotest.translate().find("*")==-1):
        hairpin5pr.append(str(seqtotest))
        currentcount=currentcount+1
        print(currentcount)
    else:
        print('-')


'''

nostruct=['CTACAACACAAGCGTTCGTTACCCTGCCAC',
 'CCCCTCGTTGTTCGCCTCCATTTCTCTCCT',
 'TGCCACCGTGCCATTGAGTCCATTTGTCCG',
 'AAATTCTCTCTATCGCCCTTGTTCATTAAC',
 'AGATGCAGCACCGACCCCTCCCACCACTCG',
 'AGCCCACATTCAATTTCTCTTGCACCTCCC']
 
 hairpin3pr=['ACGAAGGTACCACCCGTTCTCAACGGGTGG',
 'CACACGCGTCTTGCTAGAATGTCTAGCAAG',
 'ACAAATGGCTCGATTTCCCGGGGAAATCGA',
 'CCCGCTTGCAGCGCATGGGTCCCATGCGCT',
 'TCCGGGGCCAGAGCTCCGCTACGGAGCTCT',
 'AATTATCTATTGGTGAGGAGTCCTCACCAA']
 
 hairpin5pr=['CCTACACCAGCCTGGTGTAGGGGATCCCTA',
 'TGCCGTCCGCCACGGACGGCATTGTCACTA',
 'GCTACTCGCCGGGCGAGTAGCTGTCCCCCC',
 'GCGCATGTAGTTTACATGCGCGGGTCCTCC',
 'TATGCAAGCTACGCTTGCATACTCAGTCAC',
 'AGCAGGGGACCATCCCCTGCTACGTGCCGT']
'''

for i in fsmain.index:
    for j in nostruct:
        fs_upstreamvar=fs_upstreamvar.append(fsmain.loc[i],ignore_index=True)
        fs_upstreamvar.varseq[m]=Seq(j) + fsmain.varseq[i][30:]
        fs_upstreamvar.loc[m,'changes']='subgroup 4: upstream sequence no secondary structure'
        m=m+1
    for j in hairpin3pr:
        fs_upstreamvar=fs_upstreamvar.append(fsmain.loc[i],ignore_index=True)
        fs_upstreamvar.varseq[m]=Seq(j) + fsmain.varseq[i][30:]
        fs_upstreamvar.loc[m,'changes']='subgroup 4: upstream sequence hairpin at 3prime end, until pos 30'
        m=m+1
    for j in hairpin5pr:
        fs_upstreamvar=fs_upstreamvar.append(fsmain.loc[i],ignore_index=True)
        fs_upstreamvar.varseq[m]=Seq(j) + fsmain.varseq[i][30:]
        fs_upstreamvar.loc[m,'changes']='subgroup 4: upstream sequence hairpin at 5prime end, from pos 0'
        m=m+1
    

altstarts=['gaagctgctgcaagagaagctgcagctagggaggctgcagctagggaggctgctgcaaga',
 'CTACAACACAAGCGTTCGTTACCCTGCCAC',
 'TGCCACCGTGCCATTGAGTCCATTTGTCCG',
 'AGATGCAGCACCGACCCCTCCCACCACTCG']

for i in fsmain.index:
    for j in altstarts:
        for k in range(1,6,1):
            fs_upstreamvar=fs_upstreamvar.append(fsmain.loc[i],ignore_index=True)
            fs_upstreamvar.varseq[m]=Seq(j[:6*k]) + fsmain.varseq[i][6*k:]
            fs_upstreamvar.loc[m,'changes']='subgroup 5: replace beginning with constant region, 1xEAAAR and 3xfrom nostruct'
            m=m+1


fs_upstreamvar.to_csv('./design/fs_upstreamvar.csv')

#%%
twelve=pd.Series()
m=0

for i in codonsnostop:
    print(i)
    for j in codonsnostop:
        for k in codonsnostop:
            for l in codonsnostop:
                twelve.loc[m]=i+j+k+l
                m=m+1

for i in twelve:
    if ((Seq('T' + i).translate().find('*')>-1)|(Seq('TA' + i).translate().find('*')>-1)):
        twelve.drop(i,inplace=True)


nine=pd.Series()
m=0

for i in codonsnostop:
    for j in codonsnostop:
        for k in codonsnostop:
            nine.loc[m]=i+j+k
            m=m+1

for i in nine:
    if ((Seq('T' + str(i)).translate().find('*')>-1)|(Seq('TA' + str(i)).translate().find('*')>-1)):
        nine.drop(i,inplace=True)


#%% structure - from forFSstruct.py


f=open('./design/varseqplus70nteachside.txt','r')
fsmainstruct=list(SeqIO.parse(f,'fasta'))
f.close()

fsstruct=pd.DataFrame()

for i in range(len(fsmainstruct)):
    fsstruct.loc[fsmainstruct[i].id,'maxseq']=str(fsmainstruct[i].seq)

for i in fsstruct.index:
    [fsstruct.loc[i,'structfull'],fsstruct.loc[i,'MFEfull']]=RNA.fold(str(fsstruct.loc[i,'maxseq']))
    [fsstruct.loc[i,'struct20'],x]=RNA.fold(str(fsstruct.loc[i,'maxseq'][20:-20]))
    [fsstruct.loc[i,'struct40'],x]=RNA.fold(str(fsstruct.loc[i,'maxseq'][40:-40]))
    [fsstruct.loc[i,'struct60'],x]=RNA.fold(str(fsstruct.loc[i,'maxseq'][60:-60]))
    [fsstruct.loc[i,'structafterslip'],fsstruct.loc[i,'MFEfull']]=RNA.fold(str(fsstruct.loc[i,'maxseq'][112:-70]))    
    fsstruct.loc[i,'structfromfull']=fsstruct.loc[i,'structfull'][112:-70]
    fsstruct.loc[i,'structfrom20']=fsstruct.loc[i,'struct20'][92:-50]
    fsstruct.loc[i,'structfrom40']=fsstruct.loc[i,'struct40'][72:-30]
    fsstruct.loc[i,'structfrom60']=fsstruct.loc[i,'struct60'][52:-10]

fsstruct.to_csv('./design/fsstruct_forfsmain_fordifferent windows.csv')

f=open('./design/structures_from_different_windows.txt','w')

for i in fsstruct.index:
    f.write('>'+ i + '\n' + fsstruct.loc[i, 'structafterslip'] + '\n' + fsstruct.loc[i,'structfromfull'] + '\n' + fsstruct.loc[i,'structfrom20'] + '\n' + fsstruct.loc[i,'structfrom40'] + '\n' + fsstruct.loc[i,'structfrom60'] +'\n')
    
f.close()

#%% structures for stop codon mutants

for i in fsmain.index:
    [fsmain.loc[i,'structfullWT'],fsmain.loc[i,'MFEfullWT']]=RNA.fold(str(fsmain.loc[i,'varseq']))
    [fsmain.loc[i,'structafterslipWT'],fsmain.loc[i,'MFEafterslipWT']]=RNA.fold(str(fsmain.loc[i,'varseq'][42:]))    
    [fsmain.loc[i,'structfullSTOPMUT'],fsmain.loc[i,'MFEfullSTOPMUT']]=RNA.fold(str(fsmain.loc[i,'varseqstopmut']))
    [fsmain.loc[i,'structafterslipSTOPMUT'],fsmain.loc[i,'MFEafterslipSTOPMUT']]=RNA.fold(str(fsmain.loc[i,'varseqstopmut'][42:]))    

f=open('./design/fsmain_structures_withandwithout_nostop.txt','w')
for i in fsmain.index:
    f.write('>'+str(i) + '\n' + str(fsmain.structfullWT[i])+'\n'  + str(fsmain.structfullSTOPMUT[i])+'\n' + '>'+str(i) + 'afterslip\n' +  str(fsmain.structafterslipWT[i])+'\n' + str(fsmain.structafterslipSTOPMUT[i])+'\n')

f.close()

#%% DOWNSTREAM
        
nostructlong=[]
numberreturned=20
inputstring='........................................................................................................................'

commandlinestring=str("python ./antaRNA.py -Cstr '" + inputstring + "' -tGC 0.4 -tGCmax 0.8")

currentcount=0
  
while (currentcount < numberreturned):
    outputstring=subprocess.check_output(commandlinestring, shell=True).splitlines()
    seqtotest=Seq(outputstring[1]).back_transcribe()
    if (seqtotest.translate().find("*")==-1)&(seqtotest.translate()[1:].find("*")==-1)&(seqtotest.translate()[2:].find("*")==-1):
        nostructlong.append(str(seqtotest))
        currentcount=currentcount+1
        print(currentcount)
    else:
        print('-')

### select sequences to use.
#%%
fs_downstrvar=pd.DataFrame(columns=list(fsmain.columns))
m=0

forminus1=['AACCAGACCTTCTATATCGCCCTCTCCCCATCCCATCATCGACCTCCCGTTACAACCCATCACAAACCCACCAACCTCAAACATCCGATCCCCCAACCATCCCTGATTACCCACCTACAG',
           'ATCCCCGCCTTCACCCAACAACCACCCCCACTCCCCACCAATCCCACCCCTACCGTCTATCTCCCCACACACGTCACCCGCCACCACCCCGAACATCTTCCACCCACCACCCCGTCCCGA',
           'ccgaagctgctgcaagagaagctgcagctagggaggctgcagctagggaggctgctgcaa']
forplus1=['ATCCCCGCCCCAGCCCAACAACCACCCACCCTCCCCACCAATCCCACCCCTACCGTCTATCTCCCCACACACGTCACCCGCCACCACCCCGAACATCTTCCACCCACCACCCCGTCCCGA',
          'ATACCAGTCCACTCCACTCCCCCTCCTAACCCCGCCACCCCACCCTCCAGATCCATGCCTCCTCCCGCCCCCTTCCACACCCCTCCCACCACACCCACCCACACCCTTCCGTCAACACAC',
          'gaagctgctgcaagagaagctgcagctagggaggctgcagctagggaggctgctgcaaga']
          
for i in fsmain.index:
    if (fsmain.type[i]==-1):
        for j in forminus1:
            for k in range(1,len(j)//12+1,1):
                fs_downstrvar=fs_downstrvar.append(fsmain.loc[i],ignore_index=True)
                fs_downstrvar.varseq[m]=fsmain.varseq[i][:-k*12] + Seq(j[-k*12:])
                fs_downstrvar.loc[m,'changes']= 'subgroup 1: add constant sequence from the back, replace last ' + str(k*12) + 'nts'
                m=m+1
    if (fsmain.type[i]!=-1):
        for j in forplus1:
            for k in range(1,len(j)//12+1,1):
                fs_downstrvar=fs_downstrvar.append(fsmain.loc[i],ignore_index=True)
                fs_downstrvar.varseq[m]=fsmain.varseq[i][:-k*12] + Seq(j[-k*12:])
                fs_downstrvar.loc[m,'changes']= 'subgroup 1: add constant sequence from the back, replace last ' + str(k*12) + 'nts'
                m=m+1
                    
# scanning mutagenesis
                    
for i in fsmain.index:
    for j in range(0,30,1):
        fs_downstrvar=fs_downstrvar.append(fsmain.loc[i],ignore_index=True)
        fs_downstrvar.varseq[m]=fsmain.varseq[i][:42+6+j] + Seq('C') + fsmain.varseq[i][42+6+j+1:]
        fs_downstrvar.loc[m,'changes']= 'subgroup 2: C at position ' + str(6+j+1) + ' after slippery site'
        m=m+1
        fs_downstrvar=fs_downstrvar.append(fsmain.loc[i],ignore_index=True)
        fs_downstrvar.varseq[m]=fsmain.varseq[i][:42+6+j] + Seq('A') + fsmain.varseq[i][42+6+j+1:]
        fs_downstrvar.loc[m,'changes']= 'subgroup 2: A at position ' + str(6+j+1) + ' after slippery site'
        m=m+1
        
                                 

viennaanta=pd.read_pickle('./design/fs_afterfs_vienna_antarna.pkl')


for i in fsmain.index:    
    for j in viennaanta.index:
        if i in j:
            fs_downstrvar=fs_downstrvar.append(fsmain.loc[i],ignore_index=True)
            fs_downstrvar.varseq[m]=fsmain.varseq[i][:42] + viennaanta.varseq_afterfs[j]
            fs_downstrvar.loc[m,'changes']= 'subgroup 3: vienna-antarna ' + j
            m=m+1
     
     
     
pkissanta=pd.read_table('./design/fs_afterfs_pkiss_antarna_combined.txt', header=None)
pkissanta.columns=['ind','varseq_afterfs']


for i in fsmain.index:    
    for j in pkissanta.index:
        if i in pkissanta.ind[j]:
            fs_downstrvar=fs_downstrvar.append(fsmain.loc[i],ignore_index=True)
            fs_downstrvar.varseq[m]=fsmain.varseq[i][:42] + Seq(str(pkissanta.varseq_afterfs[j])).back_transcribe()
            fs_downstrvar.loc[m,'changes']= 'subgroup 3: pkiss-antarna ' + pkissanta.ind[j]
            m=m+1
 

hivstructvarnoGU=pd.read_table('./design/fs_hiv_afterfs_pkiss_antarna_structvar.txt', header=None)
hivstructvarnoGU.columns=['ind','varseq_afterfs']

for j in hivstructvarnoGU.index:
    fs_downstrvar=fs_downstrvar.append(fsmain.loc['hiv'],ignore_index=True)
    fs_downstrvar.varseq[m]=fsmain.varseq[i][:42] + Seq(str(hivstructvarnoGU.varseq_afterfs[j])).back_transcribe()
    fs_downstrvar.loc[m,'changes']= 'subgroup 4: hiv struct vars (no GU pairing) pkiss-antarna ' + hivstructvarnoGU.ind[j]
    m=m+1

hivstructvar=pd.read_pickle('./design/fs_changedstruct_hiv.pkl')

for i in fsmain.index: 
    if fsmain.type[i]==-1:
        for j in hivstructvar.index:
            fs_downstrvar=fs_downstrvar.append(fsmain.loc[i],ignore_index=True)
            fs_downstrvar.varseq[m]=fsmain.varseq[i][:42] + Seq(str(hivstructvar.varseq[j])).back_transcribe()
            fs_downstrvar.loc[m,'changes']= 'subgroup 4: hiv struct vars pkiss-antarna ' + hivstructvar.structure[j]
            m=m+1


srvlstructvar=pd.read_pickle('./design/fs_changedstruct_srvl.pkl')

for i in fsmain.index: 
    if fsmain.type[i]==-1:
        for j in srvlstructvar.index:
            fs_downstrvar=fs_downstrvar.append(fsmain.loc[i],ignore_index=True)
            fs_downstrvar.varseq[m]=fsmain.varseq[i][:42] + Seq(str(srvlstructvar.varseq[j])).back_transcribe()
            fs_downstrvar.loc[m,'changes']= 'subgroup 4: srvl struct vars pkiss-antarna ' + srvlstructvar.structure[j]
            m=m+1


nsp2fstructvar=pd.read_pickle('./design/fs_changedstruct_nsp2f.pkl')

for i in fsmain.index: 
    if fsmain.type[i]!=-1:
        for j in nsp2fstructvar.index:
            fs_downstrvar=fs_downstrvar.append(fsmain.loc[i],ignore_index=True)
            fs_downstrvar.varseq[m]=fsmain.varseq[i][:42] + Seq(str(nsp2fstructvar.varseq[j])).back_transcribe()
            fs_downstrvar.loc[m,'changes']= 'subgroup 4: nsp2f struct vars pkiss-antarna ' + nsp2fstructvar.structure[j]
            m=m+1



#  Check for Restriction enzyme binding sites

for i in fs_downstrvar.index:
    if (fs_downstrvar.varseq[i][41:].translate().find("*")>-1)&(fs_downstrvar.type[i]==-1):
        print(i)
        print(fs_downstrvar.loc[i,'sequence ID'])
        print(fs_downstrvar.changes[i])
        fs_downstrvar.drop(i,inplace=True)

for i in fs_downstrvar.index:
    if (fs_downstrvar.varseq[i][40:].translate().find("*")>-1)&(fs_downstrvar.type[i]!=-1):
        print(i)
        print(fs_downstrvar.loc[i,'sequence ID'])
        print(fs_downstrvar.changes[i])
        fs_downstrvar.drop(i,inplace=True)


for i in fs_downstrvar.index:
    if (len(fs_downstrvar.varseq[i])!=162):
        print(i)

from Bio.Restriction import RsrII
from Bio.Restriction import AscI
        
for i in fs_downstrvar.index:
    if (RsrII.search(fs_downstrvar.varseq[i])!=[])|(AscI.search(fs_downstrvar.varseq[i])!=[]):
        print(i)
        print(fs_downstrvar.changes[i])
        fs_downstrvar.drop(i,inplace=True)
       
fs_downstrvar.to_csv('./design/fs_downstrvar.csv')
        
        
#%% SYNTHETIC

fs_synthvar=pd.DataFrame(columns=list(fsmain.columns))
m=0

artificial1=pd.read_pickle('./design/fs_changedstruct_artificial1.pkl')
artificial2=pd.read_pickle('./design/fs_changedstruct_artificial2.pkl')

artificial=pd.concat([artificial1,artificial2],ignore_index=True)

for i in fsmain.index: 
    for j in artificial.index:
        fs_synthvar=fs_synthvar.append(fsmain.loc[i],ignore_index=True)
        fs_synthvar.varseq[m]=fsmain.varseq[i][:42] + Seq(str(artificial.varseq[j])).back_transcribe()
        fs_synthvar.loc[m,'changes']= 'endogenous upstream and slippery with synthetic downstream struct vars antarna ' + artificial.structure[j]
        m=m+1

altstarts=['gaagctgctgcaagagaagctgcagctagggaggctgcagctagggaggctgctgcaaga',
 'CTACAACACAAGCGTTCGTTACCCTGCCAC',
 'TGCCACCGTGCCATTGAGTCCATTTGTCCG',
 'AGATGCAGCACCGACCCCTCCCACCACTCG']

altslips=['GCTAATTTTTTT','ATTGATCCTTTT','acatgggttttt']

fsmain.loc['synth','organism']='synthetic'
for i in altstarts:
    for j in altslips:
        for k in artificial.index:
            fs_synthvar=fs_synthvar.append(fsmain.loc['synth'],ignore_index=True)
            fs_synthvar.varseq[m]=Seq(i[:30]) + Seq(j) + Seq(str(artificial.varseq[k])).back_transcribe()
            fs_synthvar.slippery[m]=j
            fs_synthvar.loc[m,'changes']= 'completely synthetic, downstream struct vars antarna ' + artificial.structure[k]
            m=m+1

fsmain.drop('synth',inplace=True)



from Bio.Restriction import RsrII
from Bio.Restriction import AscI
        
for i in fs_synthvar.index:
    if (RsrII.search(fs_synthvar.varseq[i])!=[])|(AscI.search(fs_synthvar.varseq[i])!=[]):
        print(i)
        print(RsrII.search(fs_synthvar.varseq[i]))
        print(AscI.search(fs_synthvar.varseq[i]))
        print(fs_synthvar.changes[i])
        fs_synthvar.drop(i,inplace=True)

      
fs_synthvar.to_csv('./design/fs_synthvar.csv')


#%% fs mimicking, test slippery in different frame

for i in fs.index:
    if fs.type[i]==0:
        fs.type[i]=1

fs_frameshifts=pd.DataFrame(columns=list(fs.columns))
m=0


for i in fs.index:
    if (fs.type[i]==-1) & (fs.m1stopafterfs[i]==True):
        fs_frameshifts=fs_frameshifts.append(fs.loc[i],ignore_index=True)
        fs_frameshifts.varseq[m]=fs.varseq[i][:39] + fs.varseq[i][38:-1]
        fs_frameshifts.loc[m,'changes']='frameshift mimicking'
        m=m+1
    elif (fs.type[i]==1) & (fs.p1stopafterfs[i]==True):
        fs_frameshifts=fs_frameshifts.append(fs.loc[i],ignore_index=True)
        fs_frameshifts.varseq[m]=fs.varseq[i][:39] + fs.varseq[i][40:] + Seq('C')
        fs_frameshifts.loc[m,'changes']='frameshift mimicking'
        m=m+1
    elif (fs.type[i]==-2) & (fs.p1stopafterfs[i]==True):
        fs_frameshifts=fs_frameshifts.append(fs.loc[i],ignore_index=True)
        fs_frameshifts.varseq[m]=fs.varseq[i][:39] + fs.varseq[i][37:-2]
        fs_frameshifts.loc[m,'changes']='frameshift mimicking'
        m=m+1


for i in fsmain.index:
    if (fsmain.type[i]==-1) & (fsmain.m1stopafterfs[i]==True):
        fs_frameshifts=fs_frameshifts.append(fsmain.loc[i],ignore_index=True)
        fs_frameshifts.varseq[m]=fsmain.varseq[i][:39] + fsmain.varseq[i][38:-1]
        fs_frameshifts.loc[m,'changes']='frameshift mimicking'
        m=m+1
    elif (fsmain.type[i]==1) & (fsmain.p1stopafterfs[i]==True):
        fs_frameshifts=fs_frameshifts.append(fsmain.loc[i],ignore_index=True)
        fs_frameshifts.varseq[m]=fsmain.varseq[i][:39] + fsmain.varseq[i][40:] + Seq('C')
        fs_frameshifts.loc[m,'changes']='frameshift mimicking'
        m=m+1
    elif (fsmain.type[i]==-2) & (fsmain.p1stopafterfs[i]==True):
        fs_frameshifts=fs_frameshifts.append(fsmain.loc[i],ignore_index=True)
        fs_frameshifts.varseq[m]=fsmain.varseq[i][:39] + fsmain.varseq[i][37:-2]
        fs_frameshifts.loc[m,'changes']='frameshift mimicking'
        m=m+1

for i in fsmain.index:
    if (fsmain.type[i]==-1) & (fsmain.m1stopafterfs[i]==True):
        fs_frameshifts=fs_frameshifts.append(fsmain.loc[i],ignore_index=True)
        fs_frameshifts.varseq[m]=fsmain.varseq[i][:39] + fsmain.varseq[i][38:-1]
        fs_frameshifts.loc[m,'changes']='frameshift mimicking'
        m=m+1
    elif (fsmain.type[i]==1) & (fsmain.p1stopafterfs[i]==True):
        fs_frameshifts=fs_frameshifts.append(fsmain.loc[i],ignore_index=True)
        fs_frameshifts.varseq[m]=fsmain.varseq[i][:39] + fsmain.varseq[i][40:] + Seq('C')
        fs_frameshifts.loc[m,'changes']='frameshift mimicking'
        m=m+1
    elif (fsmain.type[i]==-2) & (fsmain.p1stopafterfs[i]==True):
        fs_frameshifts=fs_frameshifts.append(fsmain.loc[i],ignore_index=True)
        fs_frameshifts.varseq[m]=fsmain.varseq[i][:39] + fsmain.varseq[i][37:-2]
        fs_frameshifts.loc[m,'changes']='frameshift mimicking'
        m=m+1

for i in fsmain.index:
    if (fsmain.type[i]==-1) & (fsmain.m1stopafterfs[i]==True):
        fs_frameshifts=fs_frameshifts.append(fsmain.loc[i],ignore_index=True)
        fs_frameshifts.varseq[m]=fsmain.varseq[i][:29] + fsmain.varseq[i][30:42] + Seq('C') +  fsmain.varseq[i][42:]
        fs_frameshifts.loc[m,'changes']='slippery shifted -1'
        m=m+1
    elif (fsmain.type[i]!=-1) & (fsmain.p1stopafterfs[i]==True):
        fs_frameshifts=fs_frameshifts.append(fsmain.loc[i],ignore_index=True)
        fs_frameshifts.varseq[m]=fsmain.varseq[i][:29] + fsmain.varseq[i][30:42] + Seq('C') +  fsmain.varseq[i][42:]
        fs_frameshifts.loc[m,'changes']='slippery shifted -1'
        m=m+1

    if (fsmain.type[i]==-1) & (fsmain.m1stopafterfs[i]==True):
        fs_frameshifts=fs_frameshifts.append(fsmain.loc[i],ignore_index=True)
        fs_frameshifts.varseq[m]=fsmain.varseq[i][:28] + fsmain.varseq[i][30:42] + Seq('CA') +  fsmain.varseq[i][42:]
        fs_frameshifts.loc[m,'changes']='slippery shifted -2'
        m=m+1
    elif (fsmain.type[i]!=-1) & (fsmain.p1stopafterfs[i]==True):
        fs_frameshifts=fs_frameshifts.append(fsmain.loc[i],ignore_index=True)
        fs_frameshifts.varseq[m]=fsmain.varseq[i][:28] + fsmain.varseq[i][30:42] + Seq('CA') +  fsmain.varseq[i][42:]
        fs_frameshifts.loc[m,'changes']='slippery shifted -2'
        m=m+1

    if (fsmain.type[i]==-1) & (fsmain.m1stopafterfs[i]==True):
        fs_frameshifts=fs_frameshifts.append(fsmain.loc[i],ignore_index=True)
        fs_frameshifts.varseq[m]=fsmain.varseq[i][:27] + fsmain.varseq[i][30:42] + Seq('CAG') +  fsmain.varseq[i][42:]
        fs_frameshifts.loc[m,'changes']='slippery shifted -3'
        m=m+1
    elif (fsmain.type[i]!=-1) & (fsmain.p1stopafterfs[i]==True):
        fs_frameshifts=fs_frameshifts.append(fsmain.loc[i],ignore_index=True)
        fs_frameshifts.varseq[m]=fsmain.varseq[i][:27] + fsmain.varseq[i][30:42] + Seq('CAG') +  fsmain.varseq[i][42:]
        fs_frameshifts.loc[m,'changes']='slippery shifted -3'
        m=m+1

    if (fsmain.type[i]==-1) & (fsmain.m1stopafterfs[i]==True):
        fs_frameshifts=fs_frameshifts.append(fsmain.loc[i],ignore_index=True)
        fs_frameshifts.varseq[m]=fsmain.varseq[i][:30]+ Seq('C') + fsmain.varseq[i][30:42]  +  fsmain.varseq[i][43:]
        fs_frameshifts.loc[m,'changes']='slippery shifted +1'
        m=m+1
    elif (fsmain.type[i]!=-1) & (fsmain.p1stopafterfs[i]==True):
        fs_frameshifts=fs_frameshifts.append(fsmain.loc[i],ignore_index=True)
        fs_frameshifts.varseq[m]=fsmain.varseq[i][:30] + Seq('C')  + fsmain.varseq[i][30:42]+  fsmain.varseq[i][43:]
        fs_frameshifts.loc[m,'changes']='slippery shifted +1'
        m=m+1

    if (fsmain.type[i]==-1) & (fsmain.m1stopafterfs[i]==True):
        fs_frameshifts=fs_frameshifts.append(fsmain.loc[i],ignore_index=True)
        fs_frameshifts.varseq[m]=fsmain.varseq[i][:30] + Seq('CA') + fsmain.varseq[i][30:42] +  fsmain.varseq[i][44:]
        fs_frameshifts.loc[m,'changes']='slippery shifted +2'
        m=m+1
    elif (fsmain.type[i]!=-1) & (fsmain.p1stopafterfs[i]==True):
        fs_frameshifts=fs_frameshifts.append(fsmain.loc[i],ignore_index=True)
        fs_frameshifts.varseq[m]=fsmain.varseq[i][:30] + Seq('CA')  + fsmain.varseq[i][30:42]+  fsmain.varseq[i][44:]
        fs_frameshifts.loc[m,'changes']='slippery shifted +2'
        m=m+1

    if (fsmain.type[i]==-1) & (fsmain.m1stopafterfs[i]==True):
        fs_frameshifts=fs_frameshifts.append(fsmain.loc[i],ignore_index=True)
        fs_frameshifts.varseq[m]=fsmain.varseq[i][:30] + Seq('CAG') + fsmain.varseq[i][30:42] +  fsmain.varseq[i][45:]
        fs_frameshifts.loc[m,'changes']='slippery shifted +3'
        m=m+1
    elif (fsmain.type[i]!=-1) & (fsmain.p1stopafterfs[i]==True):
        fs_frameshifts=fs_frameshifts.append(fsmain.loc[i],ignore_index=True)
        fs_frameshifts.varseq[m]=fsmain.varseq[i][:30] + Seq('CAG') + fsmain.varseq[i][30:42] +  fsmain.varseq[i][45:]
        fs_frameshifts.loc[m,'changes']='slippery shifted +3'
        m=m+1


for i in fsmain.index:
    fsmain.loc[i,'varseqstopmut']=fs_mutatedstops.varseq[fs_mutatedstops.loc[:,'sequence ID']==fsmain.loc[i,'sequence ID']][0]


for i in fsmain.index:
    if (fsmain.type[i]==-1) & (fsmain.m1stopafterfs[i]==True):
        fs_frameshifts=fs_frameshifts.append(fsmain.loc[i],ignore_index=True)
        fs_frameshifts.varseq[m]=fsmain.varseq[i][:29] + fsmain.varseq[i][30:42] + Seq('C') +  fsmain.varseqstopmut[i][42:]
        fs_frameshifts.loc[m,'changes']='slippery shifted -1'
        m=m+1
    elif (fsmain.type[i]!=-1) & (fsmain.p1stopafterfs[i]==True):
        fs_frameshifts=fs_frameshifts.append(fsmain.loc[i],ignore_index=True)
        fs_frameshifts.varseq[m]=fsmain.varseq[i][:29] + fsmain.varseq[i][30:42] + Seq('C') +  fsmain.varseqstopmut[i][42:]
        fs_frameshifts.loc[m,'changes']='slippery shifted -1'
        m=m+1

    if (fsmain.type[i]==-1) & (fsmain.m1stopafterfs[i]==True):
        fs_frameshifts=fs_frameshifts.append(fsmain.loc[i],ignore_index=True)
        fs_frameshifts.varseq[m]=fsmain.varseq[i][:28] + fsmain.varseq[i][30:42] + Seq('CA') +  fsmain.varseqstopmut[i][42:]
        fs_frameshifts.loc[m,'changes']='slippery shifted -2'
        m=m+1
    elif (fsmain.type[i]!=-1) & (fsmain.p1stopafterfs[i]==True):
        fs_frameshifts=fs_frameshifts.append(fsmain.loc[i],ignore_index=True)
        fs_frameshifts.varseq[m]=fsmain.varseq[i][:28] + fsmain.varseq[i][30:42] + Seq('CA') +  fsmain.varseqstopmut[i][42:]
        fs_frameshifts.loc[m,'changes']='slippery shifted -2'
        m=m+1

    if (fsmain.type[i]==-1) & (fsmain.m1stopafterfs[i]==True):
        fs_frameshifts=fs_frameshifts.append(fsmain.loc[i],ignore_index=True)
        fs_frameshifts.varseq[m]=fsmain.varseq[i][:27] + fsmain.varseq[i][30:42] + Seq('CAG') +  fsmain.varseqstopmut[i][42:]
        fs_frameshifts.loc[m,'changes']='slippery shifted -3'
        m=m+1
    elif (fsmain.type[i]!=-1) & (fsmain.p1stopafterfs[i]==True):
        fs_frameshifts=fs_frameshifts.append(fsmain.loc[i],ignore_index=True)
        fs_frameshifts.varseq[m]=fsmain.varseq[i][:27] + fsmain.varseq[i][30:42] + Seq('CAG') +  fsmain.varseqstopmut[i][42:]
        fs_frameshifts.loc[m,'changes']='slippery shifted -3'
        m=m+1

    if (fsmain.type[i]==-1) & (fsmain.m1stopafterfs[i]==True):
        fs_frameshifts=fs_frameshifts.append(fsmain.loc[i],ignore_index=True)
        fs_frameshifts.varseq[m]=fsmain.varseq[i][:30]+ Seq('C') + fsmain.varseq[i][30:42]  +  fsmain.varseqstopmut[i][43:]
        fs_frameshifts.loc[m,'changes']='slippery shifted +1'
        m=m+1
    elif (fsmain.type[i]!=-1) & (fsmain.p1stopafterfs[i]==True):
        fs_frameshifts=fs_frameshifts.append(fsmain.loc[i],ignore_index=True)
        fs_frameshifts.varseq[m]=fsmain.varseq[i][:30] + Seq('C')  + fsmain.varseq[i][30:42]+  fsmain.varseqstopmut[i][43:]
        fs_frameshifts.loc[m,'changes']='slippery shifted +1'
        m=m+1

    if (fsmain.type[i]==-1) & (fsmain.m1stopafterfs[i]==True):
        fs_frameshifts=fs_frameshifts.append(fsmain.loc[i],ignore_index=True)
        fs_frameshifts.varseq[m]=fsmain.varseq[i][:30] + Seq('CA') + fsmain.varseq[i][30:42] +  fsmain.varseqstopmut[i][44:]
        fs_frameshifts.loc[m,'changes']='slippery shifted +2'
        m=m+1
    elif (fsmain.type[i]!=-1) & (fsmain.p1stopafterfs[i]==True):
        fs_frameshifts=fs_frameshifts.append(fsmain.loc[i],ignore_index=True)
        fs_frameshifts.varseq[m]=fsmain.varseq[i][:30] + Seq('CA')  + fsmain.varseq[i][30:42]+  fsmain.varseqstopmut[i][44:]
        fs_frameshifts.loc[m,'changes']='slippery shifted +2'
        m=m+1

    if (fsmain.type[i]==-1) & (fsmain.m1stopafterfs[i]==True):
        fs_frameshifts=fs_frameshifts.append(fsmain.loc[i],ignore_index=True)
        fs_frameshifts.varseq[m]=fsmain.varseq[i][:30] + Seq('CAG') + fsmain.varseq[i][30:42] +  fsmain.varseqstopmut[i][45:]
        fs_frameshifts.loc[m,'changes']='slippery shifted +3'
        m=m+1
    elif (fsmain.type[i]!=-1) & (fsmain.p1stopafterfs[i]==True):
        fs_frameshifts=fs_frameshifts.append(fsmain.loc[i],ignore_index=True)
        fs_frameshifts.varseq[m]=fsmain.varseq[i][:30] + Seq('CAG') + fsmain.varseq[i][30:42] +  fsmain.varseqstopmut[i][45:]
        fs_frameshifts.loc[m,'changes']='slippery shifted +3'
        m=m+1


for i in fs.index:
    for j in range(0,3,1):
        fs_frameshifts=fs_frameshifts.append(fs.loc[i],ignore_index=True)
        fs_frameshifts.loc[m,'changes']='endogenous sequence'
        m=m+1
    
for i in fs.index:
    for j in range(0,3,1):
        fs_frameshifts=fs_frameshifts.append(fs.loc[i],ignore_index=True)
        fs_frameshifts.varseq[m]=fs_mutatedstops.varseq[i]
        fs_frameshifts.loc[m,'changes']='stop codons mutated'
        m=m+1
    


for i in fs_frameshifts.index:
    if (fs_frameshifts.varseq[i].translate().find("*")>-1):
        print(i)
        print(fs_frameshifts.loc[i,'sequence ID'])
 #       fs_frameshifts.drop(i,inplace=True)

for i in fs_frameshifts.index:
    if (fs_frameshifts.varseq[i][41:].translate().find("*")>-1)&(fs_frameshifts.type[i]==-1):
        print('index' + str(i))
        print(fs_frameshifts.varseq[i][41:].translate().find('*'))
        print(fs_frameshifts.loc[i,'sequence ID'])
 #       fs_frameshifts.drop(i,inplace=True)

for i in fs_frameshifts.index:
    if (fs_frameshifts.varseq[i][40:].translate().find("*")>-1):
        print(i)
 #       fs_frameshifts.drop(i,inplace=True)


for i in fs_frameshifts.index:
    if (len(fs_frameshifts.varseq[i])!=162):
        print(i)

      
fs_frameshifts.to_csv('./design/fs_frameshifts.pkl')



#%% Select HIV natural isolates (from HIV database)

toopen=open('./design/hiv-db.fasta')

hiv1db=list(SeqIO.parse(toopen,"fasta"))

hivdb=pd.DataFrame(columns=['varseq','id'])

for i in range(0,len(hiv1db)-1,1):
    start=hiv1db[i].seq.find('taatttttta')
    hivdb.loc[i,'varseq']=hiv1db[i].seq[start-32:start+130]
    hivdb.id[i]=hiv1db[i].id

hivdb.drop_duplicates('varseq',inplace=True)    

import re

for i in hivdb.index:
    if (re.findall(r'[qweryuiopsdfhjklzxvbnm]', hivdb.loc[i,'varseq'].tostring())!=[]):     
        hivdb.drop(i,inplace=True)

for i in hivdb.index:
    hivdb.loc[i,'fssite']=hivdb.varseq[i][24:84]

import random

hivdbvar=pd.DataFrame(columns=list(hivdb.columns))

for i in range(0,581,1):
    k=random.choice(hivdb.index)
    hivdbvar=hivdbvar.append(hivdb.loc[k])
    hivdb.drop(k,inplace=True)        

hivdbvar.drop_duplicates('varseq',inplace=True)

# Check for absence of restriction enzyme recognition sequences and stop codons in the relevant frames

from Bio.Restriction import RsrII
from Bio.Restriction import AscI

for i in hivdbvar.index:
    if (RsrII.search(hivdbvar.varseq[i])!=[])|(AscI.search(hivdbvar.varseq[i])!=[]):
        hivdbvar.drop(i,inplace=True)
        
for i in hivdbvar.index:
    if (len(hivdbvar.varseq[i])!=162):
        hivdbvar.drop(i,inplace=True)
        
for i in hivdbvar.index:
    if (hivdbvar.varseq[i].translate().find("*")>-1)|(hivdbvar.varseq[i][38:].translate().find("*")>-1):
        hivdbvar.drop(i,inplace=True)

for i in range(0,6,1):
    k=random.choice(hivdb.index)
    hivdbvar=hivdbvar.append(hivdb.loc[k])
    hivdb.drop(k,inplace=True)        
    
hivdbvar.to_csv('./design/hivdbvar.csv')

#%% COMBINE LIBRARIES


fs_permutations=pd.read_csv('./design/fs_permutations.pkl')
fs_slipperyvar=pd.read_csv('./design/fs_slipperyvar.pkl')
fsmain=pd.read_csv('./design/fsmain.pkl')
fs_spacervar=pd.read_csv('./design/fs_spacervar.pkl')
fs_upstreamvar=pd.read_csv('./design/fs_upstreamvar.pkl')
fs_downstrvar=pd.read_csv('./design/fs_downstrvar.pkl')
fs_synthvar=pd.read_csv('./design/fs_synthvar.pkl')
fs_frameshifts=pd.read_csv('./design/fs_frameshifts.pkl')
hivdbvar=pd.read_csv('./design/hivdbvar.csv')

fs_permutations.loc[:,'library']='fs_designed'
fs_permutations.loc[:,'subset']='fs_permutations'

fs_slipperyvar.loc[:,'library']='fs_designed'
fs_slipperyvar.loc[:,'subset']='fs_slipperyvar'

fsmain.loc[:,'library']='fs_designed'
fsmain.loc[:,'subset']='fsmain'

fs_spacervar.loc[:,'library']='fs_designed'
fs_spacervar.loc[:,'subset']='fs_spacervar'

fs_upstreamvar.loc[:,'library']='fs_designed'
fs_upstreamvar.loc[:,'subset']='fs_upstreamvar'

fs_downstrvar.loc[:,'library']='fs_designed'
fs_downstrvar.loc[:,'subset']='fs_downstrvar'

fs_synthvar.loc[:,'library']='fs_designed'
fs_synthvar.loc[:,'subset']='fs_synthvar'

fs_frameshifts.loc[:,'library']='fs_designed'
fs_frameshifts.loc[:,'subset']='fs_frameshifts'

hivdbvar.loc[:,'library']='fs_explorative'
hivdbvar.loc[:,'subset']='hivdbvar'


fslibrary=pd.concat([fs_permutations,fs_slipperyvar,fs_spacervar,fs_upstreamvar,
                     fs_downstrvar,fs_synthvar,fs_frameshifts,fsmain,fsmain,fsmain],ignore_index=True)

for i in fslibrary.index:
    if ('U' in fslibrary.varseq[i]):
        fslibrary.varseq[i]=fslibrary.varseq[i].back_transcribe()

from Bio.Restriction import RsrII
from Bio.Restriction import AscI

for i in fslibrary.index:
    if (RsrII.search(fslibrary.varseq[i])!=[])|(AscI.search(fslibrary.varseq[i])!=[]):
        fslibrary.drop(i,inplace=True)
        
for i in fslibrary.index:
    if (len(fslibrary.varseq[i])!=162):
        fslibrary.drop(i,inplace=True)
        
for i in fslibrary.index:
    if (fslibrary.varseq[i][:39].translate().find("*")>-1):
        fslibrary.drop(i,inplace=True)



# load barcodes
'''
f=open('./design/barcodes_pool_12bps_3SNPs_filtered.tab",'r')

bcfs=f.read().splitlines()    

f.close()

# All primers used for library amplification (FS and splicing)

primers=['GCCCCACGGAGGTGCCAC','CCTCCTCACGGCGACGCG','CTCCCGGGCATGCGAATT','TCAACCAGTCGCGGTCCA',
         'TGCGAGTTAGGGGACGGT','ACGGACGCGGGTATAGCA','CGAAATGGGCCGCATTGC','CACTGCGGCTGATGACGA',
            'GACAGATGCGCCGTGGAT','AGCCACCCGATCCAATGC','ATGGGGTTCGGTATGCGC','AAGGCTCCCCGAGACGAT']

def addbcfsdes(vs):
    while True:
        bc=random.choice(bcfs)
        vsnew=primers[0] + bc+str(vs) + primers[1]
        if (RsrII.search(Seq(vsnew))==[])&(AscI.search(Seq(vsnew))==[]):
            bcfs.remove(bc)
            break
    return vsnew;
    
fslibrary.varseq=fslibrary.varseq.apply(lambda x: addbcfsdes(x))

def addbcfsexp(vs):
    while True:
        bc=random.choice(bcfs)
        vsnew=primers[2] + bc+str(vs) + primers[3]
        if (RsrII.search(Seq(vsnew))==[])&(AscI.search(Seq(vsnew))==[]):
            bcfs.remove(bc)
            break
    return vsnew;
    
hivdbvar.varseq=hivdbvar.varseq.apply(lambda x: addbcfsexp(x))


fslibraryall=pd.concat([fslibrary,hivdbvar],ignore_index=True)

fslibraryall.to_csv('./design/fslibraryall.csv')
'''