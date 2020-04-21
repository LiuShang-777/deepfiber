import sys
import random
#输入文件为提取了sep.txt的fa序列
#输出文件为进行shuffle过的fa序列
input_file_list=sys.argv[1]
output_file_list=sys.argv[2]
with open(input_file_list,'r') as file:
    input_list=[]
    for line in file:
        line=line.strip()
        input_list.append(line)
with open(output_file_list,'r') as file:
    output_list=[]
    for line in file:
        line=line.strip()
        output_list.append(line)
def get_list(input_file):
    with open(input_file,'r') as file:
        list_id,list_seq=[],[]
        for line in file:
            line=line.strip()
            if line[0]=='>':
                list_id.append(line)
            else:
                list_seq.append(line)
    return (list_id,list_seq)
def sishuffle(list_id,list_seq):
    list_shuffle=[]
    for i,j in zip(list_id,list_seq):
        chrlist=list(j)
        random.shuffle(chrlist)
        list_shuffle.append(''.join(chrlist))
    return (list_id,list_shuffle)
def shuffle(input_file,output_file):    
    list_id,list_seq=get_list(input_file)
    list_id,list_new_seq=sishuffle(list_id,list_seq)
    with open(output_file,'w') as file:
        for i,j in zip(list_id,list_new_seq):
            file.write(i+'\n')
            file.write(j+'\n')
for input_file,output_file in zip(input_list,output_list):
    shuffle(input_file,output_file)
