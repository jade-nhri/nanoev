#!/usr/bin/env python3
import os,sys
import argparse
import subprocess
import pandas as pd
import shutil
import time
import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument('-i', help='an input folder')
parser.add_argument('-o', help='an output folder')
parser.add_argument('-l', help='minimum length (default=6000)')
parser.add_argument('-t', help='threads (default=8)')
parser.add_argument('-c', help='read count (default=200)')

args = parser.parse_args()

refgenomes='/opt/nanoev/EV.fasta.concor_20210303.fasta'
comm='makeblastdb -in {0} -dbtype nucl'.format(refgenomes)
subprocess.getoutput(comm)


minlen=6000
t=8
readcount=200


argv=sys.argv
if '-i' in argv:
    indir=argv[argv.index('-i')+1]
if '-o' in argv:
    outdir=argv[argv.index('-o')+1]
if '-l' in argv: 
    minlen=int(argv[argv.index('-l')+1]) 
if '-t' in argv: 
    t=int(argv[argv.index('-t')+1]) 
if '-c' in argv:
    readcount=int(argv[argv.index('-c')+1])


indir=os.path.abspath(indir)
outdir=os.path.abspath(outdir)

if not os.path.exists(outdir):
    os.mkdir(outdir)

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

def preprocess(bc,minlen):
    i=bc+'.fasta'
    comm='cat {0} | seqkit sort -l -r -w0 | head -n {1} > {2}/{3}_clean.fa'.format(i,readcount*2,outdir,bc)
        #print (comm)
    subprocess.getoutput(comm)

    comm='minimap2 -x ava-ont -t {2} {0}/{1}_clean.fa {0}/{1}_clean.fa > {0}/{1}_map.paf'.format(outdir,bc,t)
    #print (comm)
    subprocess.getoutput(comm)

    comm="cat {0}/{1}_map.paf | awk '$1!=$6 && $11>={2}' > {0}/{1}_out.paf".format(outdir,bc,minlen)
    #print (comm)
    subprocess.getoutput(comm)

def reprocess(bc,minlen):
    comm='minimap2 -x ava-ont -t {2} {0}/{1}_clean.fa {0}/{1}_clean.fa > {0}/{1}_map.paf'.format(outdir,bc,t)
    #print (comm)
    subprocess.getoutput(comm)

    comm="cat {0}/{1}_map.paf | awk '$1!=$6 && $11>={2}' > {0}/{1}_out.paf".format(outdir,bc,minlen)
    #print (comm)
    subprocess.getoutput(comm)


def removerefprocess(bc,minlen):
    i=bc+'.fasta'

    remove_names=[]
    f=open(os.path.join(outdir,bc+'_clean_ref.fa'))
    with f as fp:
        for name,seq in read_fasta(fp):
            remove_names.append(name)
    f.close()
    #print (remove_names)

    comm='rm {0}/{1}_clean.ref.fa.*'.format(outdir,bc)
    subprocess.getoutput(comm)

    for j in remove_names:
        comm='sed -e "/{0}/,+1d" {2}/{1}_clean.fa > {2}/{1}_clean_tmp.fa'.format(j,bc,outdir)
        #print (comm)
        subprocess.getoutput(comm)

    shutil.move('{0}/{1}_clean_tmp.fa'.format(outdir,bc),'{0}/{1}_clean.fa'.format(outdir,bc))
   
    comm='minimap2 -x ava-ont -t {2} {0}/{1}_clean.fa {0}/{1}_clean.fa > {0}/{1}_map.paf'.format(outdir,bc,t)
    #print (comm)
    subprocess.getoutput(comm)

    comm="cat {0}/{1}_map.paf | awk '$1!=$6 && $11>={2}' > {0}/{1}_out.paf".format(outdir,bc,minlen)
    #print (comm)
    subprocess.getoutput(comm)


def getconsensus(bc):
    if os.path.exists(os.path.join(outdir,bc)):
        comm='rm {0}/{1} -rf'.format(outdir,bc)
        subprocess.getoutput(comm)
    if os.path.exists(os.path.join(outdir,bc+'_homopolish')):
        comm='rm {0}/{1}_homopolish'.format(outdir,bc)
        subprocess.getoutput(comm)

    ref=os.path.join(outdir,bc+'_clean_ref.fa')
    if os.path.exists(ref) and os.path.getsize(ref)>0:
        comm='medaka_consensus -i {0}/{2}_clean.fa -d {1} -o {0}/{2} -t {3}'.format(outdir,ref,bc,t)
        #print (comm)
        subprocess.getoutput(comm)

        for j in range(1,6):
            sseq1=''
            sseq2=''
            comm='cp {0}/{1}/consensus.fasta {0}/{1}_ref_{2}.fa'.format(outdir,bc,j)
            #print (comm)
            subprocess.getoutput(comm)
            comm='rm {0}/{1} -rf'.format(outdir,bc)
            #print (comm)
            subprocess.getoutput(comm)
            comm='rm {0}/{1}_ref_{2}.fa*'.format(outdir,bc,j-1)
            subprocess.getoutput(comm)

            comm='medaka_consensus -i {0}/{1}_clean.fa -d {0}/{1}_ref_{2}.fa -o {0}/{1} -t {3}'.format(outdir,bc,j,t)
            #print (comm)
            stdout=subprocess.getoutput(comm)
            #print (stdout)
            f1=open(os.path.join(outdir,bc+'_ref_'+str(j)+'.fa'))
            d1=dict()
            with f1 as fp:
                for name1, seq1 in read_fasta(fp):
                    #print (seq1)
                    sseq1=sseq1+''.join(seq1)
            #print (sseq1)
            f1.close()
            f2=open(os.path.join(outdir,bc+'/consensus.fasta'))
            d2=dict()
            with f2 as fp:
                for name2, seq2 in read_fasta(fp):
                    #print (seq2)
                    sseq2=sseq2+''.join(seq2)
            #print (sseq2)
            f2.close()
 
            if sseq1==sseq2:
                #print (j)
                break
        comm='getfinal.py {0}/{1} {2}'.format(outdir,bc,minlen)
        #print (comm)
        subprocess.getoutput(comm)


        comm='cp {0}/{1}/{1}_final.fasta {0}/{1}_final.fasta'.format(outdir,bc)
        #print (comm)
        subprocess.getoutput(comm)
        
        if os.path.exists('{0}/{1}_final.fasta'.format(outdir,bc)) and os.path.getsize('{0}/{1}_final.fasta'.format(outdir,bc))>0:
            comm='runpolish.py {0} {1}_final.fasta'.format(outdir,bc)
            #print (comm)
            subprocess.getoutput(comm)

            comm='cp {0}/{1}_homopolish/{1}_final.fasta {0}/{1}_final.fasta'.format(outdir,bc)
            #print (comm)
            subprocess.getoutput(comm)


def getevinfo(bc,iteri):
    comm='getrefs.py -i {0}/{1}_clean.fa -m {0}/{1}_out.paf'.format(outdir,bc) 
    #print (comm)
    stdout=subprocess.getoutput(comm)
    if os.path.exists('{0}/{1}_clean_ref.fa'.format(outdir,bc)) and os.path.getsize('{0}/{1}_clean_ref.fa'.format(outdir,bc))==0:
        preprocess(bc,minlen-1500)
        comm='getrefs.py -i {0}/{1}_clean.fa -m {0}/{1}_out.paf'.format(outdir,bc)
        #print (comm)
        stdout=subprocess.getoutput(comm)

    #print (stdout)
    getconsensus(bc)

    if os.path.exists('{0}/{1}_final.fasta'.format(outdir,bc)) and os.path.getsize('{0}/{1}_final.fasta'.format(outdir,bc))!=0:
        comm='blastn -query {0}/{1}_final.fasta -db {2} -out {0}/{1}_type.txt -num_threads 6 -max_target_seqs 1 -outfmt "6 qseqid stitle pident length"'.format(outdir,bc,refgenomes)
        #print (comm)
        subprocess.getoutput(comm)
        comm='cat {0}/{1}_type.txt'.format(outdir,bc)
        #print (comm)
        stdout=subprocess.getoutput(comm)
        print ('\t{0}'.format(stdout))

    if os.path.exists('{0}/{1}_type.txt'.format(outdir,bc)) and os.path.getsize('{0}/{1}_type.txt'.format(outdir,bc))!=0:
        if iteri==0:
            comm='gettable.py {0}/{1}_type.txt {0}/{1}.csv'.format(outdir,bc)
            #print (comm)
            stdout=subprocess.getoutput(comm)
        else:
            comm='gettable.py {0}/{1}_type.txt {0}/{1}_{2}.csv'.format(outdir,bc,iteri)
            #print (comm)
            stdout=subprocess.getoutput(comm)



os.chdir(indir)
mybc=[x for x in os.listdir() if '.fasta' in x]
for i in sorted(mybc):
    bc=i.replace('.fasta','')
    if os.path.exists('{0}/{1}_final.fasta'.format(outdir,bc)):
        continue
    else:
        print (bc)
        preprocess(bc,minlen)
        getevinfo(bc,0)
    j=0
    breakloop=False
    while not breakloop:
        j+=1
        if j==3:
            breakloop=True
    #for j in range (1,4) and not breakloop==True:

        if not os.path.exists('{0}/{1}.csv'.format(outdir,bc)) and os.path.getsize('{0}/{1}_clean_ref.fa'.format(outdir,bc))>0:
            print ('{0} try again'.format(bc))
            comm='cp {0}/{1}_clean.fa {0}/{1}_old_clean.fa'.format(outdir,bc)
            subprocess.getoutput(comm)

            removerefprocess(bc,minlen)
            getevinfo(bc,j)


        if os.path.exists('{0}/{1}.csv'.format(outdir,bc)) and os.path.getsize('{0}/{1}.csv'.format(outdir,bc))>0:
            dfcsv=pd.read_table('{0}/{1}.csv'.format(outdir,bc),sep=',')
            #print (dfcsv)
            temp1=dfcsv['Length'].values
            temp2=dfcsv['Len_polyprotein'].values
            temp3=dfcsv['Len_uORF'].values
            #print ('{0} {1} {2}'.format(temp1,temp2,temp3))
            for t1,t2,t3 in zip(temp1,temp2,temp3):
                #print ('{0} {1} {2}'.format(t1,t2,t3))
                if t1<6000 or t2<2000 or t3<40:
                    print ('{0} try again'.format(bc))
                    comm='cp {0}/{1}_clean.fa {0}/{1}_old_clean.fa'.format(outdir,bc)
                    subprocess.getoutput(comm)

                    removerefprocess(bc,minlen)
                    getevinfo(bc,j)
                    break


        if os.path.exists('{0}/{1}_{2}.csv'.format(outdir,bc,j)) and os.path.getsize('{0}/{1}_{2}.csv'.format(outdir,bc,j))>0:
            dfcsv=pd.read_table('{0}/{1}_{2}.csv'.format(outdir,bc,j),sep=',')
            #print (dfcsv)
            temp1=dfcsv['Length'].values
            temp2=dfcsv['Len_polyprotein'].values
            temp3=dfcsv['Len_uORF'].values
            for t1,t2,t3 in zip(temp1,temp2,temp3):
                if t1>6000 and t2>2000 and t3>40:
                    breakloop=True
                    break
                else:
                    temp4=dfcsv['polyprotein'].values[0]
                    #print (temp4)
                    if os.path.exists('{0}/{1}_{2}.csv'.format(outdir,bc,j-1)) and os.path.getsize('{0}/{1}_{2}.csv'.format(outdir,bc,j-1))>0:
                        predfcsv=pd.read_table('{0}/{1}_{2}.csv'.format(outdir,bc,j-1),sep=',')
                        pretemp4=predfcsv['polyprotein'].values[0]
                        if temp4==pretemp4 and temp4!='':
                            breakloop=True
                            break

os.chdir(outdir)
print (outdir)
rerunbcs=[x for x in os.listdir() if '_3.csv' in x]
print (rerunbcs)
for i in rerunbcs:
    mybc=i.replace('_3.csv','')
    comm='cat {0}/{1}.fasta | seqkit sort -l -r -w0 | head -n {2} > {3}/{1}_clean.fa'.format(indir,mybc,readcount*2,outdir)
    print (comm)
    subprocess.getoutput(comm)
    comm='binfa.py -i {0}/{1}_clean.fa -o {0}'.format(outdir,mybc)
    print (comm)
    subprocess.getoutput(comm)
    comm='rm {0}*.csv'.format(mybc)
    print (comm)
    subprocess.getoutput(comm)

myfiles=[x for x in os.listdir() if '-1_clean.fa' in x or '-2_clean.fa' in x]
print (myfiles)
for i in myfiles:
    bc=i.split('_')[0]
    print (bc)
    reprocess(bc,minlen)
    getevinfo(bc,0)
    j=0
    breakloop=False
    while not breakloop:
    #for j in range (1,4):
        j+=1
        if j==3:
            breakloop=True

        if not os.path.exists('{0}/{1}.csv'.format(outdir,bc)) and os.path.getsize('{0}/{1}_clean_ref.fa'.format(outdir,bc))>0:
            print ('{0} try again'.format(bc))
            comm='cp {0}/{1}_clean.fa {0}/{1}_old_clean.fa'.format(outdir,bc)
            subprocess.getoutput(comm)

            removerefprocess(bc,minlen)
            getevinfo(bc,j)


        if os.path.exists('{0}/{1}.csv'.format(outdir,bc)) and os.path.getsize('{0}/{1}.csv'.format(outdir,bc))>0:
            dfcsv=pd.read_table('{0}/{1}.csv'.format(outdir,bc),sep=',')
            #print (dfcsv)
            temp1=dfcsv['Length'].values
            temp2=dfcsv['Len_polyprotein'].values
            temp3=dfcsv['Len_uORF'].values
            if temp1<6000 or temp2<2000 or temp3<40:
                print ('{0} try again'.format(bc))
                comm='cp {0}/{1}_clean.fa {0}/{1}_old_clean.fa'.format(outdir,bc)
                subprocess.getoutput(comm)

                removerefprocess(bc,minlen)
                getevinfo(bc,j)
                break

        if os.path.exists('{0}/{1}_{2}.csv'.format(outdir,bc,j)) and os.path.getsize('{0}/{1}_{2}.csv'.format(outdir,bc,j))>0:
            dfcsv=pd.read_table('{0}/{1}_{2}.csv'.format(outdir,bc,j),sep=',')
            #print (dfcsv)
            temp1=dfcsv['Length'].values
            temp2=dfcsv['Len_polyprotein'].values
            temp3=dfcsv['Len_uORF'].values
            if temp1>6000 and temp2>2000 and temp3>40:
                breakloop=True
                break
            else:
                temp4=dfcsv['polyprotein'].values[0]
                #print (temp4)
                if os.path.exists('{0}/{1}_{2}.csv'.format(outdir,bc,j-1)) and os.path.getsize('{0}/{1}_{2}.csv'.format(outdir,bc,j-1))>0:
                    predfcsv=pd.read_table('{0}/{1}_{2}.csv'.format(outdir,bc,j-1),sep=',')
                    pretemp4=predfcsv['polyprotein'].values[0]
                    if temp4==pretemp4 and temp4!='':
                        breakloop=True
                        break


comm='gentable.py ./'
subprocess.getoutput(comm)
