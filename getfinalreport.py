#!/usr/bin/env python3
import os,sys
import pandas as pd
import numpy as np
import subprocess

indir=os.path.abspath(sys.argv[1])
df1=pd.read_csv(os.path.join(indir,'report.csv'))
print (df1)
df2=pd.read_csv(os.path.join(indir,'report_rerun.csv'))
print (df2)
df3=pd.read_csv(os.path.join(indir,'report2.csv'))
print (df3)

df=df1.append(df2)
df=df.append(df3)
print (df)

unisample=np.unique(df['Sample'].values)
sample=[]
for i in sorted(unisample):
    if 'Sample' not in i:
        sample.append(i)
print (sample)
print (len(sample))

count=0
uorf=[]
pp=[]
lenpp=[]
lenuorf=[]
for i in sample:
    count+=1
    dfset=df[df['Sample']==i]
    #print (dfset)
    if len(dfset)>1:
        dfset1=dfset.astype({'Length':float})
        temp=dfset1.sort_values(by=['Length'],ascending=False, na_position='last')
        tmp=temp.head(1)
    else:
        tmp=dfset
        #print (tmp)
    if count==1:
        dfout=tmp
    else:
        dfout=dfout.append(tmp)
    comm='/opt/getorf.py {0}/{1}_final.fasta'.format(indir,i)
    #print (comm)
    stdout=subprocess.getoutput(comm)
    print (stdout)
    tmp=stdout.split('\n')
    if '>' in tmp[0] and '>' in tmp[2]:
        pp.append(tmp[1])
        if len(tmp)==4:
            uorf.append(tmp[3])
        else:
            uorf.append('')
    else:
        pp.append('')
        uorf.append('')
#print (uorf)
#print (pp)
for i,j in zip(pp,uorf):
    lenpp.append(len(i))
    lenuorf.append(len(j))
#print (dfout)

dfout1=dfout.copy()
dfout1['Polyprotein']=pp
dfout1['Len_polyprotein']=lenpp
dfout1['uORF']=uorf
dfout1['Len_uORF']=lenuorf
print (dfout1)

dfout.to_csv(os.path.join(indir,'final_report.csv'),index=False)
dfout1.to_csv(os.path.join(indir,'final_report_new.csv'),index=False)

