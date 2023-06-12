#!/usr/bin/env python3
import os,sys
import subprocess
import pandas as pd

fqdir=os.path.abspath(sys.argv[1])
indir=os.path.abspath(sys.argv[2])
minlen=int(sys.argv[3])
t=int(sys.argv[4])

df=pd.read_csv(os.path.join(indir,'report.csv'))
#print (df)

dfset=df[df['Length']<minlen]
#print (dfset)
handeling=dfset['Sample'].values
print (handeling)

os.chdir(indir)
exlist=[]
for i in range(0,5):
    for j in handeling:
        if j in exlist:
            continue
        if os.path.exists(j+'_final.fasta'):
            comm='mv {0}_final.fasta {0}_{1}.final.fasta'.format(j,i)
            #print (comm)
            subprocess.getoutput(comm)
            comm='touch {0}_final.fasta'.format(j)
            #print (comm)
            subprocess.getoutput(comm)
    comm='rerun.py -i {0} -o {1} -l {2} -t {3}'.format(fqdir,indir,minlen,t)
    #print (comm)
    stdout=subprocess.getoutput(comm)
    print (stdout)

    if os.path.exists(os.path.join(indir,'report_rerun.csv')):
        dfr=pd.read_csv(os.path.join(indir,'report_rerun.csv'))
        dfr=dfr[dfr['Sample']!='Sample']
        dfr=dfr.astype({'Length':float})
        dfrset=dfr[dfr['Length']<minlen]
        dfrset1=dfr[dfr['Length']>=minlen]
        handelingnew=set(dfrset['Sample'].values)
        exlist=set(dfrset1['Sample'].values)
        hd=handelingnew-exlist
        print (hd)
        for k in hd:
            if k in handeling:
                continue
            if os.path.exists(k+'_final.fasta'):
                comm='mv {0}_final.fasta {0}_{1}.final.fasta'.format(k,i)
                print (comm)
                subprocess.getoutput(comm)
                comm='touch {0}_final.fasta'.format(k)
                print (comm)
                subprocess.getoutput(comm)



