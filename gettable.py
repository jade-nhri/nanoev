#!/usr/bin/env python3 
import pandas as pd 
import numpy as np 
import sys, os 
from Bio.Seq import Seq 
from Bio import SeqIO 
infile=sys.argv[1] 
outfile=sys.argv[2] 
dftype=pd.read_table(infile,sep='\t',names=['ReadID','Subject','Identity','Length']) 
print (dftype) 
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
def getorf(infile): 
    outseq=dict() 
    for record in SeqIO.parse(infile,'fasta'): 
        print (record.id) 
        seq='' 
        seq1='' 
        seq2='' 
        for j in list(record): 
            seq=seq+''.join(j) 
        myseq=Seq(seq) 
        #print (myseq) 
        seq='' 
        cds='' 
        for i in range (0,3): 
            seqf=myseq[i:] 
            cdsf=seqf.translate() 
            if (cdsf.count('*')<=40): 
                #print (i) 
                #print (cdsf) 
                seq=seqf 
                cds=cdsf 
        if seq=='' and cds=='': 
            myseqr=myseq.reverse_complement() 
            for i in range (0,3): 
                seqr=myseqr[i:] 
                cdsr=seqr.translate() 
                if (cdsr.count('*')<=40): 
                    seq=seqr 
                    cds=cdsr 
        #print (cds) 
        oricds=cds
        polyprotein=''
        for orf in cds.split('*'): 
            if len(orf)>1500: 
                polyprotein=orf 
                #print ('>'+record.id+'_polyprotein') 
                #print (polyprotein[polyprotein.find('M'):]) 
                seq1=str(polyprotein[polyprotein.find('M'):]) 
                #print (seq1) 
        if polyprotein!='':
            orfpos=cds.find(polyprotein)*3
        else:
            orfpos=201
        print (orfpos)
        if orfpos<200: 
            upstream=seq[orfpos-150:orfpos+150] 
        else: 
            upstream=seq[:orfpos+150] 
        if upstream!='':
            for i in range (0,3): 
                upseq=upstream[i:] 
                uorf=upseq.translate() 
                for upep in uorf.split('*'): 
                    if len(upep)>50 and 'MVT' in upep: 
                        #print (upep) 
                        upos=upep.find('MVT') 
                        UP=upep[upos:] 
                        seq2=str(UP)
        else:
            UP=''
        if record.id not in outseq.keys(): 
            outseq[record.id]=dict()
            outseq[record.id]['Seq']=seq
            outseq[record.id]['polyprotein']=seq1 
            outseq[record.id]['uorf']=seq2 
    #print (outseq) 
    return(outseq) 
dfuni=dftype.drop_duplicates(subset=['ReadID']) 
print (dfuni) 
subject=dfuni.Subject.values 
print (subject) 
evtype=[] 
for i in subject: 
    tmp=i.split('|') 
    temp='{0}|{1}'.format(tmp[1],tmp[2]) 
    print (temp) 
    evtype.append(temp) 
print (evtype) 
d=dict() 
fafile=infile.replace('_type.txt','_final.fasta') 
f=open(fafile) 
with f as fp: 
    for name, seq in read_fasta(fp): 
        name=name.replace('>','') 
        d[name]=seq 
f.close() 
print (d) 
name=[] 
seq=[] 
for i in d.keys(): 
    sample=os.path.split(infile)[-1] 
    name.append(sample.replace('_type.txt','')) 
    seq.append(d[i]) 
dfuni.insert(0,'Sample',name,True) 
dfuni.insert(2,'Type',evtype,True) 
dfuni.insert(6,'Sequence',seq,True) 
#dfout=dfuni[['Sample','ReadID','Type','Identity','Length','Sequence']] 
outseq=getorf(fafile) 
dfseq=pd.DataFrame.from_dict(outseq,orient='index') 
dfseq['ReadID']=dfseq.index 
dfseq['Len_polyprotein']=dfseq['polyprotein'].str.len() 
dfseq['Len_uORF']=dfseq['uorf'].str.len() 
dfout=pd.merge(dfuni,dfseq,left_on='ReadID',right_on='ReadID') 
dfout=dfout[['Sample','ReadID','Type','Identity','Length','Seq','polyprotein','Len_polyprotein','uorf','Len_uORF']] 
#dfout.to_csv(outfile,mode='a',header=False,index=False) 
dfout.to_csv(outfile,index=False) 
print (dfout)
