#!/usr/bin/env python3
import os,sys
import pandas as pd
import numpy as np
import pysam
from collections import Counter

indir=sys.argv[1]
indir=os.path.abspath(indir)

if not sys.argv[2]:
    minlen=6500
else:
    minlen=int(sys.argv[2])


outname=indir.split('/')[-1]
print (outname)

infa=os.path.join(indir,'consensus.fasta')
inbam=os.path.join(indir,'calls_to_draft.bam')
bamfile=pysam.AlignmentFile(inbam,'rb')

d=dict()
f=open(infa)
def read_fasta(fp):
    name, seq = None, []
    for line in fp:
        line = line.rstrip()
        if line.startswith(">"):
            if name: yield (name, ''.join(seq))
            name, seq = line, []
        else:
            seq.append(line)
    if name: yield (name, ''.join(seq))

with f as fp:
    for name, seq in read_fasta(fp):
        name=name.replace('>','')
        d[name]=seq
f.close()
print (d)

dout=dict()
for i in d.keys():
    #print (i)
    pos=[]
    for pileupcolumn in bamfile.pileup(i,0,):
        #print (pileupcolumn.n)
        if pileupcolumn.n>=10:
            pos.append(int(pileupcolumn.pos))
    print (len(pos))
    if len(pos)>=minlen:
        print (len(pos))
        minpos=min(pos)
        maxpos=max(pos)
        print ('get {0} from {1} to {2}'.format(i,minpos,maxpos))
        dout[i]=d[i][minpos:maxpos]
#print (dout)

fw=open(os.path.join(indir,outname+'_final.fasta'),'w')
for i in dout.keys():
    fw.write ('>'+i+'\n')
    fw.write (dout[i]+'\n')
fw.close()

