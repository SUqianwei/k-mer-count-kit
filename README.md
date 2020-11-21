# k-mer-count-kit
#2020年11月更新：  
1、 主程序变为improved1.py。  
2、	采用argparse来进行互动，极大的方便了程序的运行。  
3、	对可处理文件的数目做了优化，更新后可以处理多个文件。   
4、 新增加两个测试数据albA.fasta和Streptococcus_pneumoniae_R6.fasta，用来测试多个输入文件的处理过程。albA.fasta是一个蛋白的DNA序列，Streptococcus_pneumoniae_R6.fasta是肺炎链球菌的基因组序列。  
5、 缺点：对多序列的fasta文件不具备处理能力。  
请到 python程序运行记录.docx 中查看详细的更新说明和运行结果。  

###参数说明：  
usage: improved1.py [-h] [--length] [--translate] [-f FILE [FILE ...]]
                    [--complement] [--kmer_counter] [--compare] [-k K_VALUE]  

optional arguments:  
  -h, --help            show this help message and exit  
  --length              统计输入序列的碱基数  
  --translate           把DNA或RNA序列翻译为蛋白质序列。  
  -f FILE [FILE ...], --file FILE [FILE ...]
                        请输入要处理的文件。   
  --complement          把输入的DNA序列按碱基互补配对原则转换成其反向互补序列。  
  --kmer_counter        请输入进行kmer计数的k值。  
  --compare             通过k-mer counter比较不同菌株基因组的相似性。  
  -k K_VALUE, --k_value K_VALUE
                        用于进行k-mer计算的mer的值。  

###示例：  
获得输入序列文件的长度：  
	python improved1.py --length -f albA.fasta albB_dna.fasta  
翻译DNA或RNA序列为蛋白质序列：  
	python improved1.py --translate -f albA.fasta albB_dna.fasta  
转换反向互补序列：  
	python improved1.py --complement -f albA.fasta albB_dna.fasta  
kmer-counter计数：  
	python improved1.py --kmer_counter -k 3 -f EcoliK12.fasta EcoliO157.fasta coelicolor.fasta  
通过k-mer counter比较不同菌株基因组的相似性:  
	python improved1.py --compare -k 3 -f EcoliK12.fasta EcoliO157.fasta coelicolor.fastafasta Streptococcus_pneumoniae_R6.fasta caulobacterNA1000.fasta  

	
#程序说明(2020年10月14日):  
用python写的一个关于k-mer counter的程序。主程序为 improved.py 。

该程序有五个功能：  
	length :统计输入序列的碱基数。  
	translate :把DNA或RNA序列翻译为蛋白质序列。  
	complement :把输入的DNA序列按碱基互补配对原则转换成其反向互补序列。  
	k_mer_counter :对输入序列的k-mer进行基本的统计。  
	compare : 通过k-mer counter比较不同菌株基因组的相似性。  

有6个测试数据：  
albB_dna.fasta ：该序列用来测试translate和complement两个功能，它是我本科的时候用过的一个蛋白序列，这里提供的是它的DNA序列，需要先把该序列转变为反向互补序列，再对得到的反向互补序列进行翻译，才能得到蛋白序列。  
caulobacterNA1000.fasta，coelicolor.fasta，EcoliK12.fasta，EcoliO157.fasta：这四个数据均可用来进行k_mer_counter功能的测试，分别为四个细菌的基因组序列。  
第四个功能compare的实现需要把以上四个文件根据提示依次输入。  
另外还有一个人的基因组数据GRCh38.p13.fasta，用该基因组来做k-mer counter 功能测试，k=10时可以得到结果。详见 python程序运行记录.docx。  

鉴于程序的运行繁琐且冗长，建议您先在 python程序运行记录.docx 中查看具体的运行结果和参数的测试结果。
如果看完 python程序运行记录.docx 后您还有兴趣运行程序，那么您需要在一个具有python3，并且安装了pandas的环境下运行该程序。
在linux中，请执行 python improved.py ,随后根据程序提示输入相关数据。  

程序自评：  
	优点：把多个功能集成到了一个程序中，方便相关功能的调用。  
	      把子序列使用特定的数字表示，并存储为列表索引，一定程度上减少了对内存的使用。这个idea来自于https://developer.51cto.com/art/201003/188595.html。
	缺点：缺点那就有很多了，这里就说两个吧。  
	      需要交互，不能放到后台运行,不能用slurm系统提交。运行时建议重开一个会话。  
	      在 compare 功能中，有且只写了四个序列的相关性分析。  
	      算法和代码都不够精妙，运行的不快，耗时。  




