

### 配置所需的软件包

* hisat2 
* samtools 
* sratoolkit 
* fastqc
* trimmomatic


### 下载SRA
1. prefetch 
>prefetch SRA*** -O 输出目录

2. docker安装kingfisher(主要是用conda依赖一堆，github下载慢)

>docker run -v $PWD wwood/kingfisher:0.1.2 get -r SRR23330487 -m ena-ascp ena-ftp prefetch aws-http

3. ascp 下不了NCBI的SRA了

### SRA转化成fq文件

* 来源于sra-tools的fastq-dump命令，-split-3 把pair-end双端测序分成两个文件输出；
--gzip 输出压缩文件，节约储存
>fastq-dump --gzip --split-3 SRA*** -O 输出文件夹名

### fq比对到基因组

构建索引：
>extract_exons.py refseq_annotation.gtf > c_elegans.exons.gtf 

>extract_splice_sites.py refseq_annotation.gtf > c_elegans.splice_sites.gtf 


>hisat2-build --ss c_elegans.splice_sites.gtf --exon c_elegans.exon.gtf GCF_000002985.6_WBcel235_genomic.fna c.elegans

比对：

>hisat2 --time --threads 10 -x /home/luotao/c.elegans/RNA-seq/ref/c.elegans -1 /home/luotao/c.elegans/RNA-seq/sra/SRR23330487/clean_data/SRR23330487_1.fq -2 /home/luotao/c.elegans/RNA-seq/sra/SRR23330487/clean_data/SRR23330487_2.fq -S /home/luotao/c.elegans/RNA-seq/sra/SRR23330487/sam/SRR23330487.sam

sam转化为bam：





















### PS
被CSDN气哭，一群傻逼下个SRA的方法过时了还抄来抄去

