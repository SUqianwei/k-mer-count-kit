# -*- coding: utf-8 -*-
"""
Created on Thu Sep 17 21:43:55 2020

@author: suqianwei
"""
# =============================================================================
# #读入原始序列
# #对原始序列进行归一化，去除大小写和第一行的影响,去除换行符，去掉所有的 N和 M（2），R（1）
# =============================================================================

def getSeq(filename):  
    import re
    import sys
    raw_seq = open(filename,"r").read()
    
    #seq = re.match(r'^>.*$',raw_seq,re.M)
    #return seq.group()
    seq_num = len(re.findall(r'^>',raw_seq,flags=re.M))
    if seq_num > 1 or seq_num == 0 :
        print('您输入的文件不符合要求，请输入含有一条序列的fasta文件。')
        sys.exit(1)
    seq = re.sub(r'^>.*$','',raw_seq,flags=re.M)
    seq = seq.strip()
    seq = seq.upper()
    seq = seq.replace('\n','').replace('\r','')
    seq = seq.replace('N','')
    seq = seq.replace('M','')
    seq = seq.replace('R','')
    return str(seq),seq_num

# =============================================================================
# 获得反向互补序列 reverse sequence
# 注释部分是我自己写的循环，运行的太慢了；
# 网上的这个方法，运行的真的超快，不知道底层原理是什么，学到了学到了。'''
# =============================================================================

#def reverse_seq(seq):
#    pairing_list = {'A':'T','T':'A','C':'G','G':'C'}
#    rev_seq_list = []    
#    for i in seq:        
#        rev_seq_list.insert(0,pairing_list[i])
#    rev_seq = ''.join(rev_seq_list)
#    return rev_seq    
def reverseSeq(seq):
    transtab = str.maketrans('ATCG','TAGC')
    rev_seq = seq.translate(transtab)
    return rev_seq[::-1]

# =============================================================================
# #把 DNA或 RNA序列翻译为蛋白质序列
# =============================================================================

def transToProtein(seq):
    #标准密码子表
    codon_table = {'GCU':'A','GCC':'A','GCA':'A','GCG':'A',
                'CGU':'R','CGC':'R','CGA':'R','CGG':'R','AGA':'R','AGG':'R',
                'UCU':'S','UCC':'S','UCA':'S','UCG':'S','AGU':'S','AGC':'S',
                'AUU':'I','AUC':'I','AUA':'I','AUU':'I','AUC':'I','AUA':'I',
                'UUA':'L','UUG':'L','CUU':'L','CUC':'L','CUA':'L','CUG':'L',
                'GGU':'G','GGC':'G','GGA':'G','GGG':'G',
                'GUU':'V','GUC':'V','GUA':'V','GUG':'V',
                'ACU':'T','ACC':'T','ACA':'T','ACG':'T',
                'CCU':'P','CCC':'P','CCA':'P','CCG':'P',
                'AAU':'N','AAC':'N',
                'GAU':'D','GAC':'D',
                'UGU':'C','UGC':'C',
                'CAA':'Q','CAG':'Q',
                'GAA':'E','GAG':'E',
                'CAU':'H','CAC':'H',
                'AAA':'K','AAG':'K',
                'UUU':'F','UUC':'F',
                'UAU':'Y','UAC':'Y',
                'AUG':'M',
                'UGG':'W',
                'UAG':'STOP','UGA':'STOP','UAA':'STOP'}
    #把 DNA序列转为 RNA序列,并去除掉起始密码子 AUG 上游的序列
    start_codon = 'AUG'
    #stop_codon = ['UAG','UGA','UAA']
    seq = str(seq)
    rna_seq = seq.translate(str.maketrans('T','U'))
    #print(rna_seq)
    try:
        rna_seq = rna_seq[rna_seq.index(start_codon):]
    except ValueError:
        print('您输入的序列中不存在起始密码子AUG。请注意，您输入的DNA序列或RNA序列必须含有起始密码子AUG或者ATG。')
        return
    
    #print('您输入序列的RNA形式为：\n%s' % rna_seq)
    #把 RNA序列翻译为蛋白质序列
    protein = ''
    
    try:
        for i in range(0,len(rna_seq),3):
            if codon_table[rna_seq[i:i+3]] == 'STOP':
                protein += '*'
                break
            else:
                protein += codon_table[rna_seq[i:i+3]]
    except KeyError:
        print('您输入的编码序列格式有误，没有终止密码子并且不是3的倍数\
              下方输出的序列为前 3n 个碱基翻译成的蛋白质序列:')
          
    return protein


# =============================================================================
# k-mer counter 模块
# 把子序列映射为列表索引，这个idea 来源于 https://developer.51cto.com/art/201003/188595.html
# 输出一个列表，列表的索引表示不同的子序列，列表的值表示该子序列出现的次数。
# =============================================================================

# 把子序列映射为连续的非负整数，该非负整数可以被作为列表索引
def motif2int(motif):      
    m2iD = dict(A=0,C=1,G=2,T=3)      
    total = 0 
    for i, letter in enumerate(motif):  
        total += m2iD[letter]*4**(len(motif)-i-1)  
    return total


# 把 motif2int 得到的非负整数作索引，输出一个列表，列表的值表示该子序列出现的次数
def toDinucleotideList(seq,k):
    di_list = [0 for index in range(4**k)]
    for i in range(len(seq)-k+1):
        index = motif2int(seq[i:i+k])
        di_list[index] += 1
    return di_list


# 把索引值转化为对应的子序列
def indexToDinucleotide(num,k):
    l = []
    
    i2mD = {0:'A', 1:'C', 2:'G', 3:'T'}
    for i in range(k):
        #print(i)
        a = num // (4**(k-i-1))
        num = num % (4**(k-i-1))
        l.append(i2mD[a])
    nucle = ''.join(l)
    return nucle


# 打开一个待写入的文本
def openFile(fn): 
    obj = open(fn,'a+')
    obj.seek(0)
    obj.truncate()
    return obj

# =============================================================================
# #通过 k-mer counter 对四条序列的相关性作分析
# =============================================================================

#获得序列的 k-mer list，并获得列名
def getSeries(name,k):
    import pandas as pd
    import re
    seq = getSeq(name)[0]
    l = toDinucleotideList(seq,k)
    colna = re.sub(r'\..*$','',name)    
    return l,colna




# =============================================================================
# #主程序，用来实现不同的功能
# =============================================================================
import re
import argparse

parser = argparse.ArgumentParser()

parser.add_argument("--length",action='store_true', help="统计输入序列的碱基数")
parser.add_argument("--translate", action='store_true', help="把DNA或RNA序列翻译为蛋白质序列。")
parser.add_argument("-f", "--file", nargs='+', help="请输入要处理的文件。")
parser.add_argument("--complement", action='store_true', help="把输入的DNA序列按碱基互补配对原则转换成其反向互补序列。")
parser.add_argument("--kmer_counter", action='store_true', help="请输入进行kmer计数的k值。")
parser.add_argument("--compare", action='store_true', help="通过k-mer counter比较不同菌株基因组的相似性。")
parser.add_argument("-k", "--k_value", type=int, help="用于进行k-mer计算的mer的值。")
args = parser.parse_args()

if args.length:
    for item in args.file:
        seq, seq_num = getSeq(item)    
        print("The length of %s:\t%i" % (item, len(seq)))
elif args.translate:
    for item in args.file:
        seq = getSeq(item)[0]
        temp = re.sub("(.*)\.fa.*","\g<1>",item) + ".protien.fasta"
        f2 = openFile(temp)
        print(">",temp, "\n", transToProtein(seq),sep="")
        print('>%s' % item, file=f2)
        print(transToProtein(seq),file=f2)
        print('您也可以在输出的文本%s文件中查看。' % temp)
        f2.close()
elif args.complement:
    for item in args.file:
        temp = re.sub("(.*)\.fa.*","\g<1>",item) + ".complement.fasta"
        f1 = openFile(temp)
        seq = getSeq(item)[0]
        print(">", temp,"\n", reverseSeq(seq), sep="")
        print('>%s_complement' % item, file=f1)        
        print(reverseSeq(seq),file=f1)
        print('请到%s文件中查看结果。' % temp)
        f1.close()
elif args.kmer_counter:
    for item in args.file:
        temp = re.sub("(.*)\.fa.*","\g<1>",item) + "_kmc.out"
        #n = int(input('请输入您的k值：\n'))
        seq = getSeq(item)[0]
        b = toDinucleotideList(seq, args.k_value)
        #f3 = openFile('temp1')
        f4 = openFile(temp)
        for i,element in enumerate(b):
            if element>1:
                #print(i,element,file=f3)
                print(indexToDinucleotide(i,args.k_value), element, file=f4)
        #f3.close()
        f4.close()
        print('请到%s中查看结果' % temp)

elif args.compare:
    import pandas as pd
    k = args.k_value 
    df = pd.DataFrame()
    for item in args.file:
        temp1,temp2=getSeries(item,k)
        df[temp2]=temp1
    result = df.corr()
    print(result)

else :
    print("您输入的功能尚未被开发。")




    





     
    






            
            

            

