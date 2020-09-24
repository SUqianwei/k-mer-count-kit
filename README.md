# k-mer-count-kit
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