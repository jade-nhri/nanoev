#!/usr/bin/python3
import sys,os,subprocess
import numpy as np
import scipy
import pandas as pd
import time
import math
from sklearn.cluster import AgglomerativeClustering
argv=sys.argv
if '-i' in argv:
    inputfile=argv[argv.index('-i')+1]
    inputfile=os.path.abspath(inputfile)
if '-o' in argv:
    outname=argv[argv.index('-o')+1]
    outname=os.path.abspath(outname)
    if not os.path.isdir(outname):
        os.mkdir(outname)
myTime=time.strftime('%Y%m%d_%H%M',time.localtime(time.time()))
mer=5
lthr=0.7
depthfile=''
ctgthr=0
folderName='%s_%dmer'%(myTime,mer)
if '-ct' in argv:
    ctgthr=argv[argv.index('-ct')+1]
    try:
        ctgthr=int(ctgthr)
    except:
        ctgthr=0
if os.path.exists(folderName):
    os.system('cd ./; rm -rf '+folderName+'')
    os.system('cd ./; mkdir '+folderName+'')
else:
    os.system('cd ./; mkdir '+folderName+'')
tmpNowdir=os.getcwd()
os.chdir(folderName)
if os.path.exists('1_Data'):
    os.system('cd ./; rm -rf 1_Data')
    os.system('cd ./; mkdir 1_Data')
else:
    os.system('cd ./; mkdir 1_Data')
os.chdir('1_Data')
#####################################################
def rename(infile,thr):
    f=open(infile)
    scf=dict()
    for i in f:
        i=i.strip()
        if '>' in i:
            h=i
            scf[h]=''
            continue
        scf[h]+=i.upper()
    f.close()
    fw=open('My.fa','w')
    fw1=open('Alias.txt','w')
    count=1
    for h in scf.keys():
        seq=scf[h]
        if len(seq)<thr:
            continue
        tmpH='>Seq_%d_len_%d\n'%(count,len(seq))
        fw.write(tmpH)
        fw1.write('%s\t%s\n'%(h,tmpH.strip()))
        fw.write('%s\n'%seq)
        count+=1
    fw.close()
    fw1.close()
    print('Seqs >= %d : %d'%(thr,count-1))
argv=sys.argv
inf=inputfile
thr=1000
rename(inf,thr)
os.chdir('../')
#####################################################
print('Get Feature.')
os.system('cd ./; mkdir 2_GetFeature')
os.chdir('2_GetFeature')
class r_comp:
    def __init__(self):
        self.c_nt=dict()
        self.c_nt['A']='T'
        self.c_nt['T']='A'
        self.c_nt['C']='G'
        self.c_nt['G']='C'
    def complement(self,seq):
        global c_nt
        s_list=list(seq)
        s_list.reverse()
        new_s_list=[]
        for i in range(len(s_list)):
            now_nt=s_list[i]
            if now_nt in['A','T','C','G']:new_s_list.append(self.c_nt[now_nt])
            else:new_s_list.append(now_nt)
        return ''.join(new_s_list)
def GetMer(myrc):
    my_mer=[]
    for i1 in ['A','T','C','G']:
        for i2 in ['A','T','C','G']:
            for i3 in ['A','T','C','G']:
                for i4 in ['A','T','C','G']:
                    for i5 in ['A','T','C','G']:
                        s=i1+i2+i3+i4+i5
                        ns=myrc.complement(s)
                        if (s in my_mer) or (ns in my_mer):
                            continue
                        my_mer.append(s)

    fw=open('2_Mer.txt','w')
    fw.write('\n'.join(my_mer))
    fw.close()
    return my_mer
def GetScfMer(myfile,myrc,myMerList):
    if not os.path.exists('Scf_layout'):
        os.mkdir('Scf_layout')
    if not os.path.exists('n_Scf_layout'):
        os.mkdir('n_Scf_layout')
    f=open(myfile)
    for i in f:
        i=i.strip()
        if '>' in i:
            filename=i.replace('>','')
            mylen=float(filename.split('_')[-1])
            continue
        seq=i
        fw=open('Scf_layout/'+filename+'.txt','w')
        fw1=open('n_Scf_layout/'+filename+'.txt','w')
        layout=[]
        for mer in myMerList:
            rmer=myrc.complement(mer)
            if mer==rmer:#palindrome coden
                tmpcount=seq.count(mer)
            else:
                tmpcountf=seq.count(mer)
                tmpcountr=seq.count(myrc.complement(mer))
                tmpcount=tmpcountr+tmpcountf
            if tmpcount==0:
                tmpcount=1
            layout.append(tmpcount)

        sumlayout=float(sum(layout))
        for tmpcount in layout:
            fw.write('%d\n'%tmpcount)
            fw1.write('%.20f\n'%(tmpcount/sumlayout))

        fw.close()
        fw1.close()
    f.close()
myrc=r_comp()
myMerList=GetMer(myrc)
print(len(myMerList))
GetScfMer('../1_Data/My.fa',myrc,myMerList)
os.chdir('../')
#####################################################
os.system('cd ./; mkdir 3_GetMatrix')
os.chdir('3_GetMatrix')
scfThr=1000
merList=[]
f=open('../2_GetFeature/2_Mer.txt')
for i in f:
    i=i.strip()
    merList.append(i)
f.close()

matrix=[]
t_matrix1=[]

for i in range(len(merList)):
    matrix.append([])
    t_matrix1.append([])


scfDict=dict()
scfSumDict=dict()
for i in os.listdir('../2_GetFeature/n_Scf_layout'):
    f=open('../2_GetFeature/n_Scf_layout/'+i)
    f1=open('../2_GetFeature/Scf_layout/'+i)
    data=[]
    data1=[]
    for j in f:
        j=j.strip()
        data.append(float(j))
    for x in f1:
        x=x.strip()
        data1.append(int(x))
    f.close()
    f1.close()
    scfDict[i]=data
    scfSumDict[i]=sum(data1)
scfList=sorted(scfSumDict, key=scfSumDict.get, reverse=True)

for i in scfList:
    data=scfDict[i]
    for j in range(len(data)):
        matrix[j].append(data[j])

fw=open('3_scf_order_by_column.txt','w')
fw.write('\n'.join(scfList))
fw.close()

fw=open('3_matrix_for_corr.txt','w')
for i in matrix:
    fw.write('\t'.join([str(x) for x in i])+'\n')
fw.close()


fw=open('3_matrix.txt','w')
for i in matrix:
    fw.write('\t'.join([str(x) for x in i])+'\n')
fw.close()

#clr Get Zi
for i in range(len(matrix)):
    data=matrix[i]
    denominator=0.0
    data1=[x for x in data if x != 0]
    denominator=sum([math.log(x) for x in data1])/float(len(data1))

    for j in data:
        if j==0:
            t_value=0.0
        else:
            t_value=math.log(j)-denominator
        t_matrix1[i].append(t_value)
fw=open('3_clr_matrix.txt','w')
for i in np.transpose(t_matrix1):
    #fw.write('\t'.join([str(x) for x in i])+'\n')
    fw.write('\t'.join(['%.20f'%x for x in i])+'\n')
fw.close()


count=0
fscf=open('3_scf_order_by_column.txt')
ffea=np.loadtxt("./3_matrix_for_corr.txt")
ffea=np.transpose(ffea)
fscfL=open('3_scf_order_by_column_long.txt','w')
fscfS=open('3_scf_order_by_column_short.txt','w')
ffeaL=open('3_matrix_long.txt','w')
ffeaS=open('3_matrix_short.txt','w')
for scf,fea in zip(fscf,ffea):
    scfLen=int(scf.split('_')[-1].replace('.txt',''))
    fea=('\t'.join([str(i) for i in fea]))+'\n'
    if (scfLen>=scfThr) or (scf.strip() in CogList):
        fscfL.write(scf)
        ffeaL.write(fea)
        count+=1
    else:
        fscfS.write(scf)
        ffeaS.write(fea)
fscf.close()
fscfL.close()
fscfS.close()
ffeaL.close()
ffeaS.close()
print(' %d contigs entering first stage clustering'%count)

fclr=open('3_clr_matrix.txt')
fscf=open('3_scf_order_by_column.txt')
fclrL=open('3_clr_matrix_long.txt','w')
fclrS=open('3_clr_matrix_short.txt','w')
for scf,clr in zip(fscf,fclr):
    scfLen=int(scf.split('_')[-1].replace('.txt',''))
    if (scfLen>=scfThr) or (scf.strip() in CogList):
        fclrL.write(clr)#fscfL.write(scf)
    else:
        fclrS.write(clr)#fscfS.write(scf)
fclr.close()
fclrL.close()
fclrS.close()
fscf.close()

#####################################################
os.chdir('../../')
d=dict()
f=open(folderName+'/1_Data/My.fa')
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
        d[name]=seq
f.close()
li = []
df = pd.read_csv(folderName+'/3_GetMatrix/3_scf_order_by_column.txt',sep='\t',names=['Data'])
df2 = pd.read_csv(folderName+'/1_Data/Alias.txt',sep='\t',names=['Before','After'])
df3 = pd.read_csv(folderName+'/3_GetMatrix/3_clr_matrix.txt',sep='\t',header=None)
n=0
for i in df['Data']:
    num = i.rfind('.')
    af = '>'+i[:num]
    nc = df2[(df2['After']==af)].index.tolist()
    bf = df2['Before'].iloc[nc[0]]
    li.append(bf)
    n+=1
lx=df3.to_numpy()
n=0
li2 = []
for i in li:
    li2.append(i)
    n+=1
model = AgglomerativeClustering(n_clusters=2, affinity='euclidean', linkage='ward').fit(lx)
#print(len(model.labels_))
m = list(set(model.labels_))
li3 = []
n = 0
for i in li2:
    k = model.labels_[n]
    li3.append([i,k])
    n+=1
dic = {}
strr = ''
for i in m:
    for j in li3:
        if j[1] == i:
            strr+=j[0]+','
    strr = strr[:-1]
    dic.setdefault(str(i),strr)
    strr = ''
for i in dic.keys():
    strr = ''
    fw = open(outname+'/group'+str(i)+'.fasta','w')
    print('label',i)
    sp = dic[i].split(',')
    print(sp)
    for j in sp:
        num = df2[df2['Before']==j].index.tolist()
        header = df2['After'].iloc[num[0]]
        strr += j + '\n' + d[header] + '\n'
    fw.write(strr)
    fw.close()
print('#####')
os.system('cd ./; rm -r '+folderName+'')

