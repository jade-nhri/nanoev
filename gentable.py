#!/usr/bin/env python3
import os,sys
import pandas as pd
import numpy as np

indir=sys.argv[1]
os.chdir(indir)
mycsv=[x for x in os.listdir() if '.csv' in x]
print (mycsv)

files=dict()
for i in sorted(mycsv):
    name=i.replace('.csv','')
    name=name.split('_')[0]
    if name not in files.keys():
        files[name]=i
    else:
        files[name]=i
print (files)


count=0
for i in files.keys():
    count+=1
    dfi=pd.read_table(files[i],sep=',')
    #print (dfi)
    if count==1:
        df=dfi
    else:
        df=pd.concat([df,dfi])
print (df)

df.to_csv('final_report.csv',index=False)
