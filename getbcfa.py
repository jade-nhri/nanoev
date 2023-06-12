#!/usr/bin/env python3
import sys, os, subprocess
import pandas as pd
import numpy as np
import io,gzip
import argparse
import multiprocessing as mp

threads=8
parser = argparse.ArgumentParser()
parser.add_argument('-m', help='the paf file produced by minimap2')
parser.add_argument('-i', help='the reads in fastq')
parser.add_argument('-o', help='an output folder')
parser.add_argument('-t', help='threads (default=8)')
args = parser.parse_args()

argv=sys.argv
if '-m' in argv:
    mappingfile=argv[argv.index('-m')+1]
if '-i' in argv:
    fqfile=argv[argv.index('-i')+1]
if '-o' in argv:
    outdir=argv[argv.index('-o')+1]
if '-t' in argv:
    threads=int(argv[argv.index('-t')+1])


mappingfile=os.path.abspath(mappingfile)
fqfile=os.path.abspath(fqfile)

if not os.path.exists(outdir):
    os.mkdir(outdir)
outdir=os.path.abspath(outdir)



df=pd.read_table(mappingfile,names=['Qname','Qlen','Qstart','Qend','strand','Tname','Tlen','Tstart','Tend','Nmatch','Alen','MQ1','MQ2','MQ3','MQ4','MQ5','MQ6','MQ7'],low_memory=False)
#print (df)
unibarcodes=np.unique(df['Tname'].values)

readc=df.groupby(['Qname']).size().reset_index(name='counts')

read1=readc[readc['counts']==1]['Qname'].values
#print (len(read1))

read2=readc[readc['counts']==2]['Qname'].values
#print (len(read2))

readm=readc[readc['counts']==3]['Qname'].values
#print (len(readm))

###One hit
bcr=dict()
dfset1=df[df['Qname'].isin(read1)]
dfset1=dfset1[dfset1['Nmatch']>=54]
#print (dfset1)
readid1=dfset1['Qname'].values
strand1=dfset1['strand'].values
tmpstart=dfset1['Qstart'].values
tmpend=dfset1['Qend'].values
bcs1=dfset1['Tname'].values

for i,j in zip(bcs1,readid1):
    if i not in bcr.keys():
        bcr[i]=[]
        bcr[i].append(j)
    else:
        bcr[i].append(j)
#print (bcr)

pos1=[]
for i,j,k in zip(strand1,tmpstart,tmpend):
    if '+' in i:
        pos1.append(k)
    if '-' in i:
        pos1.append(j)


###From multiple hits to get two hits
dfsetm=df[df['Qname'].isin(readm)]
#print (dfsetm)
dfsetm['QTname']=dfsetm['Qname']+dfsetm['Tname']
#print (dfsetm)
readc=dfsetm.groupby(['QTname']).size().reset_index(name='counts')
#print (readc)
qtname=readc[readc['counts']==2]['QTname'].values
adddfset2=dfsetm[dfsetm['QTname'].isin(qtname)]
#print (adddfset2)
add=adddfset2.drop(columns=['QTname'])
#print (add)


###Combine two hits
dfset2=df[df['Qname'].isin(read2)]
dfset2=dfset2.append(add,ignore_index=True)
#print (dfset2)

dfset2['QTname']=dfset2['Qname']+dfset2['Tname']
readc=dfset2.groupby(['QTname']).size().reset_index(name='counts')
qtname=readc[readc['counts']==2]['QTname'].values
dfset2f=dfset2[dfset2['QTname'].isin(qtname)]
dfset2=dfset2f.drop(columns=['QTname'])
#print (dfset2)
dfset2['QTname']=dfset2['Qname']+dfset2['strand']
readc=dfset2.groupby(['QTname']).size().reset_index(name='counts')
qtname=readc[readc['counts']==1]['QTname'].values
dfset2f=dfset2[dfset2['QTname'].isin(qtname)]
dfset2=dfset2f.drop(columns=['QTname'])
#print (dfset2)



dfset2tmp=dfset2.drop_duplicates(subset=['Qname'])
#print (dfset2tmp)
readid2=dfset2['Qname'].values
strand2=dfset2['strand'].values
tmpstart=dfset2['Qstart'].values
tmpend=dfset2['Qend'].values

readid22=dfset2tmp['Qname'].values
bcs2=dfset2tmp['Tname'].values

for i,j in zip(bcs2,readid22):
    j2=j+'_2'
    if i not in bcr.keys():
        bcr[i]=[]
        bcr[i].append(j2)
    else:
        bcr[i].append(j2)
#print (bcr)

pos2=[]
for i,j,k in zip(strand2,tmpstart,tmpend):
    if '+' in i:
        pos2.append(k)
    if '-' in i:
        pos2.append(j)
#print (pos2)

d=dict()
if '.gz' in fqfile:
    print ('Reading fq in gz...')
    f=io.TextIOWrapper(gzip.open(fqfile,'r'))

else:
    f=open(fqfile)

while True:
    h=f.readline()
    if not h: break
    h=h.replace('\n','')
    #print (h)
    rID=h.split()[0]
    rID=rID.replace('@','')
    #print (rID)
    seq=f.readline().replace('\n','')
    qh=f.readline().replace('\n','')
    qual=f.readline().replace('\n','')
    d[rID]=[]
    d[rID].append(h)
    d[rID].append(seq)
    d[rID].append(qual)
f.close()

def getread1(dfq,readID,strand,pos):
    d1=dict()
    for i,j,k in zip(readID,strand,pos):
        if j=='+':
            d1[i]=dfq[i][1][k:]
        if j=='-':
            d1[i]=dfq[i][1][:k]
    return(d1)

fa1=getread1(d,readid1,strand1,pos1)
print (len(fa1.keys()))

def getread2(dfq,readID,pos):
    d2=dict()
    for i in range(0,len(readID)):
        j=2*i+1
        start=min([pos[j],pos[j-1]])
        end=max([pos[j],pos[j-1]])
        #print ('{0}\t{1}\t{2}'.format(readID[i],start,end))
        j2=readID[i]+'_2'
        d2[j2]=dfq[readID[i]][1][start:end]
    return(d2)
            
fa2=getread2(d,readid22,pos2)
fa={**fa1,**fa2}
print (len(fa2.keys()))
print (len(fa.keys()))

def writefa(bc,fadict,readlist):
    fw=open(bc+'.fasta','w')
    print ('Writing reads for {0} ......'.format(bc))
    for i in readlist:
        midseq=fadict[i][100:-100]
        if 'AAAAAAAAAAAAAAAAAAAA' in midseq or 'TTTTTTTTTTTTTTTTTTTT' in midseq or len(midseq)<800:
            continue
        else:
            fw.write('>'+i+'\n')
            fw.write(fadict[i]+'\n')
    fw.close()

def main():
    po=mp.Pool(threads)
    for i in sorted(bcr.keys()):
        readlist=bcr[i]
        if len(readlist)>=80:
            po.apply_async(writefa,[i,fa,readlist])
    po.close()
    po.join()


if __name__=='__main__':
    os.chdir(outdir)
    main()

