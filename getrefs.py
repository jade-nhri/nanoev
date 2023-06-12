#!/usr/bin/python3
import sys, os, subprocess
import pandas as pd
import numpy as np
import io,gzip
import itertools


minreads=100

argv=sys.argv
if '-m' in argv:
    mappingfile=argv[argv.index('-m')+1]
if '-i' in argv:
    infile=argv[argv.index('-i')+1]
if '-n' in argv:
    minreads=int(argv[argv.index('-n')+1])


mappingfile=os.path.abspath(mappingfile)
infile=os.path.abspath(infile)


def read_fasta(fp):
    name, seq = None, []
    for line in fp:
        line = line.rstrip()
        if line.startswith(">"):
            print (line)
            tmp=line.split(' ')
            if name: yield (name, ''.join(seq))
            name, seq = tmp[0], []
        else:
            seq.append(line)
    if name: yield (name, ''.join(seq))


d=dict()
f=open(infile)
with f as fp:
    for name, seq in read_fasta(fp):
        d[name]=seq
f.close()


print (mappingfile)
df=pd.read_table(mappingfile,names=['Qname','Qlen','Qstart','Qend','strand','Tname','Tlen','Tstart','Tend','Nmatch','Alen','MQ1','MQ2','MQ3','MQ4','MQ5','MQ6'],low_memory=False)
print (df)
df=df.sort_values(by=['Qlen'],ascending=False)
print (df)
tempq=df[['Qname']].values.tolist()
tempq=list(itertools.chain.from_iterable(tempq))
print (tempq)
qid=[]
for i in tempq:
    if i not in qid:
       qid.append(i)


#print (qid)


fw=open(infile.replace('.fa','_ref.fa'),'w')
setb=set()
for id in qid:
    dfset=df[df['Qname']==id]
    Tname=dfset['Tname'].values.tolist()
    Nseq=len(Tname)
    print ('{0} with {1} sequences'.format(id,Nseq))
    #if Nseq>=minreads and '_2' in id:
    if Nseq>=minreads:
        Tname.append(id)
        seta=set(Tname)
        #print (seta)
        #print (seta&setb)
        if seta & setb ==set():
            print (id)
            setb=seta|setb
            #print (setb)
            fw.write('>'+id+'\n')
            fw.write(d['>'+id]+'\n')
fw.close()

