#!/usr/bin/env python3
import os,sys
import subprocess
import pandas as pd
import numpy as np

indir=os.path.abspath(sys.argv[1])
infile=sys.argv[2]
refgenome='/opt/nanoev/EV.fasta.concor_20210303.fasta'
refUTR='/opt/nanoev/EV71_UTR.fa'
outdir=infile.replace('_final.fasta','_homopolish')

comm='grep "TTTTTTTTTT" {0}/{1}/{2}'.format(indir,infile.replace('_final.fasta',''),infile)
#print (comm)
stdout=subprocess.getoutput(comm)
#print (stdout)
if len(stdout)>0:
    comm='seqkit seq {0}/{2}/{1} -r -p -w 0 > {0}/{1}'.format(indir,infile,infile.replace('_final.fasta',''))
    #print (comm)
    subprocess.getoutput(comm)


comm='python3 /opt/homopolish/homopolish.py polish -a {3}/{0} -l {1} -m /opt/homopolish/R9.4.pkl -o {3}/{2}'.format(infile,refgenome,outdir,indir)
print (comm)
subprocess.getoutput(comm)

comm='mv {0}/{1}/*_homopolished.fasta {0}/{1}/{2}'.format(indir,outdir,infile.replace('.fasta','_homopolished.fasta'))
print (comm)
subprocess.getoutput(comm)

comm='minimap2 {0} {1}/{2}/{3} > {1}/ori.paf'.format(refUTR,indir,outdir,infile.replace('.fasta','_homopolished.fasta'))
print (comm)
subprocess.getoutput(comm)

df=pd.read_csv('{0}/ori.paf'.format(indir),names=['Qname','Qlen','Qstart','Qend','strand','Tname','Tlen','Tstart','Tend','Nmatch','Alen','MQ1','MQ2','MQ3','MQ4','MQ5','MQ6','MQ7'],sep='\t') 
print (df) 



if not df.empty:
    strand=np.unique(df['strand'].values)
    print (strand)
    for i in strand:
        if '+' in i:
            comm='seqkit seq {0}/{1}/{2} -w 0 > {0}/{1}/{3}'.format(indir,outdir,infile.replace('.fasta','_homopolished.fasta'),infile)
            print (comm)
            subprocess.getoutput(comm)
        else:
            comm='seqkit seq {0}/{1}/{2} -r -p -w 0 > {0}/{1}/{3}'.format(indir,outdir,infile.replace('.fasta','_homopolished.fasta'),infile)
            print (comm)
            subprocess.getoutput(comm)
else:
    comm='seqkit seq {0}/{1}/{2} -w 0 > {0}/{1}/{3}'.format(indir,outdir,infile.replace('.fasta','_homopolished.fasta'),infile)
    print (comm)
    subprocess.getoutput(comm)

