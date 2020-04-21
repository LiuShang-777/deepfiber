import sys
import random
input_file=sys.argv[1]
output_file=sys.argv[2]
with open(input_file,'r') as file:
    list_all=[]
    for line in file:
        line=line.strip()
        if line[0]=='>':
            list_all.append(line)
        else:
            continue
def allocate_list(list_all):
    list1,list2=[],[]
    for i in list_all:
        if 'ae' in i:
            list1.append(i)
        else:
            list2.append(i)
    if len(list1)<=len(list2):
        return (list1,list2,'ex')
    else:
        return (list2,list1,'un')
list_total=allocate_list(list_all)
list_min,list_max=list_total[0],list_total[1]
print(len(list_min),len(list_max),list_total[2])
list_max_new=random.sample(list_max,len(list_min))
for i in list_max_new:
    list_min.append(i)
del list_max
def get_sequence(list_min,input_file):
    with open(input_file,'r') as file:
        list_seq,list_id=[],[]
        for line in file:
            line=line.strip()
            if line[0]=='>':
                list_id.append(line)
            else:
                list_seq.append(line)
    list_select_id,list_select_seq=[],[]
    for i,j in zip(list_id,list_seq):
        if i in list_min:
            list_select_id.append(i)
            list_select_seq.append(j)
    return (list_select_id,list_select_seq)
list_select_id,list_select_seq=get_sequence(list_min,input_file)
with open(output_file,'w') as file:
    for i ,j in zip(list_select_id,list_select_seq):
        file.write(i+'\n')
        file.write(j+'\n')

