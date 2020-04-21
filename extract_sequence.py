# -*- coding: utf-8 -*-
"""
Created on Thu Nov 28 17:10:21 2019

@author: liushang
"""
import numpy as np
import pandas as pd
import sys
gff_file=sys.argv[1]
genome_file=sys.argv[2]
express_file=sys.argv[3]
index_name=sys.argv[4]
with open(gff_file,'r') as file:
    list_gff=[]
    for line in file:
        line=line.strip()
        if line[0]=='#':
            pass
        else:
            line=line.split('\t')
            line[8]=line[8].split(';')[0][3:]
            list_gff.append(line)

list_utr_gene=[]
for i in list_gff:
    if (i[2]=='gene')&(i[0][0]!='S')&(i[8][5]!='S'):
        list_utr_gene.append(i)
utr_gene_nd=[]
for i in list_utr_gene:
    if i[8][:15] not in utr_gene_nd:
        utr_gene_nd.append(i[8][:15])
dic_len={}
for i in list_utr_gene:
    if int(i[3])==int(i[4]):
                continue
    if i[8][:15] in utr_gene_nd:
        if i[8][:15] in dic_len.keys():
            if (int(i[3]),int(i[4])) not in dic_len[i[8][:15]]:    
                dic_len[i[8][:15]].append((int(i[3]),int(i[4]),i[6]))
            else:
                continue
        else:
            dic_len[i[8][:15]]=[]
            dic_len[i[8][:15]].append((int(i[3]),int(i[4]),i[6]))
dic_max_len={}
for i in dic_len.keys():
    dic_max_len[i]=max([j[1]-j[0] for j in dic_len[i]])
length=[i for i in dic_max_len.values()]

upstream=1000
downstream=500
dic_result={}
mid_point=0
for i in dic_len.keys():
    temp=[]
    if dic_len[i][0][2]=='+':
        mid_point=min([j[0] for j in dic_len[i]])
        temp=[(max(0,mid_point-upstream),mid_point),(mid_point,mid_point+downstream)]
    else:
        mid_point=max([j[1] for j in dic_len[i]])
        temp=[(mid_point,mid_point+upstream),(max(0,mid_point-downstream),mid_point)]
    dic_result[i]=temp
#map
with open(genome_file,'r') as file:
    name,seq,temp=[],[],[]
    for i in file:
        i=i.strip()
        if i[0]=='>':
            name.append(i[1:])
            seq.append(temp)
            temp=[]
        else:
            temp.append(i)
    seq.append(temp)
    seq=seq[1:]
dic_chromosome={}
for i,j in zip(name,seq):
    dic_chromosome[i]=''.join(j)
del name,seq
dic_upstream={}
dic_downstream={}
def reverse_complete(s):
    s=s[::-1]
    table=str.maketrans('AGCT','TCGA')
    return s.translate(table)
for i in dic_result.keys():
    if dic_result[i][0][1]==dic_result[i][1][0]: 
        dic_upstream[i]=dic_chromosome[i[5:8]][dic_result[i][0][0]-1:dic_result[i][0][1]-1]
        dic_downstream[i]=dic_chromosome[i[5:8]][dic_result[i][1][0]-1:dic_result[i][1][1]-1]
    else:
        dic_upstream[i]=reverse_complete(dic_chromosome[i[5:8]][dic_result[i][0][0]:dic_result[i][0][1]])
        dic_downstream[i]=reverse_complete(dic_chromosome[i[5:8]][dic_result[i][1][0]:dic_result[i][1][1]])
del dic_chromosome
#no utr gene
'''
dic_gene_noutr={}
for i in list_gff:
    if (i[2]=='gene')&(i[0][0]!='S')&(i[8][5]!='S')&(i[8] not in utr_gene_nd):
       dic_gene_noutr[i[8]]=(int(i[3]),int(i[4]),i[6])
import pandas as pd
gene_expression=pd.read_csv('F:/UTR/all_fpkm.txt',sep='\t')
for i in dic_gene_noutr.keys():
    if i in dic_downstream.keys():
        print(i)
name_gene=[i for i in gene_expression['name'] if (i in dic_gene_noutr.keys())|(i in dic_downstream.keys())]

with open('F:/gene_name.txt','w') as file:
    for i in name_gene:
        file.write(i+'\n')
        '''
'''
all_ex=pd.read_csv('F:/UTR/all_express.txt',sep='\t')
all_ex_name=[i for i in all_ex['name']]
alter=pd.read_csv('F:/UTR/alter.txt',sep='\t')
alter_name=[i for i in alter['name']]
all_low=pd.read_csv('F:/UTR/all_low.txt',sep='\t')  
all_low_name=[i for i in all_low['name']] 
'''
'''
data=pd.read_csv('F:/UTR/all_express.txt',sep='\t')
'''
#with open('F:/deepfiber/all_ex.txt','r') as file:
#    all_ex_name=[i.strip() for i in file]
#with open('F:/deepfiber/all_low.txt','r') as file:
#    all_low_name=[i.strip() for i in file]
'''
with open('F:/deepfiber/alter.txt','r') as file:
    alter_name=[i.strip() for i in file]
'''
data=pd.read_csv(express_file,sep='\t')
def fetch_class(dataframe):
    list_feature=[i for i in dataframe.columns]
    list_feature=list_feature[1:]
    dic={}
    for i in list_feature:
        all_ex=[j for j in dataframe.loc[dataframe[i]==2][index_name]]
        all_low=[j for j in dataframe.loc[dataframe[i]==0][index_name]]
        dic[i]=(all_ex,all_low)
    return dic
dic_class=fetch_class(data)
for i in dic_class.keys():
    with open('%s_utr5.fa'%i,'w') as file:
        for x,j,z in zip(dic_upstream.keys(),dic_upstream.values(),dic_downstream.values()):
            if x in dic_class[i][0]:
                file.write('>'+x+'ae'+'\n')
                file.write(j+z+'\n')
            elif x in dic_class[i][1]:
                file.write('>'+x+'aw'+'\n')
                file.write(j+z+'\n')
            else:
                continue
'''            
with open('F:/deepfiber/init_utr5.fa','w') as file:
    for i,j,z in zip(dic_upstream.keys(),dic_upstream.values(),dic_downstream.values()):
        if i in all_ex_name:            
            file.write('>'+i+'ae'+'\n')
            file.write(j+z+'\n')
            
        elif i in alter_name:
            file.write('>'+i+'al'+'\n')
            
        elif i in all_low_name:
            file.write('>'+i+'aw'+'\n')
            file.write(j+z+'\n')
            
        else:
            continue        
'''    
