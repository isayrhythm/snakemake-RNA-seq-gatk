

### 配置所需的软件包

* hisat2 （比对到基因组）
* samtools （sam格式处理）
* sratoolkit （下SRA）
* fastqc （质控）
* fastp（自动化质控报告）


### 下载SRA
1. prefetch 
```
prefetch SRA*** -O 输出目录
# prefetch也很快，突然就很快
```


2. docker安装kingfisher(主要是用conda依赖一堆，github下载慢)

```
docker run -v $PWD wwood/kingfisher:0.1.2 get -r SRR23330487 -m ena-ascp ena-ftp prefetch aws-http
```

3. ascp 已经下不了NCBI的SRA了(2023-6-7)

### SRA转化成fq文件

* 来源于sra-tools的fastq-dump命令，-split-3 把pair-end双端测序分成两个文件输出；
--gzip 输出压缩文件，节约储存
```
fastq-dump --gzip --split-3 {SRA***} -O {outfile}
```

### fq比对到基因组

构建索引：
```
extract_exons.py refseq_annotation.gtf > c_elegans.exons.gtf 

extract_splice_sites.py refseq_annotation.gtf > c_elegans.splice_sites.gtf 
```

```
hisat2-build --ss c_elegans.splice_sites.gtf --exon c_elegans.exon.gtf GCF_000002985.6_WBcel235_genomic.fna c.elegans
```

**以下开始使用snakemake**

比对：

```
hisat2 --time --threads 10 -x /home/luotao/c.elegans/RNA-seq/ref/c.elegans -1 /home/luotao/c.elegans/RNA-seq/sra/SRR23330487/clean_data/SRR23330487_1.fq -2 /home/luotao/c.elegans/RNA-seq/sra/SRR23330487/clean_data/SRR23330487_2.fq -S /home/luotao/c.elegans/RNA-seq/sra/SRR23330487/sam/SRR23330487.sam
```
sam转化为bam：


```
samtools view -S {input} -b -@ {threads} > {output} 
# -S sam， -b bam ，-@ 线程
```

bam sorted，bam排序：
```
samtools sort {input} -o {output} -@ {threads}
```

### reads计数
```
htseq-count -r name -f bam -n {threads} {input} {elegans_gtf}
# -r 排序方式 -f bam 文件格式  
```

### RPKM、FPKM、TPM

* RPKM用于单端测序  
* FPKM双端测序：  
readcount消除掉样品的差异，所以除以样品的总的readcount（双端测序为fragments）  
再消除基因长度的差异，除以map到的基因的长度  
* TPM：
对每个基因的read数用基因长度进行校正，之后再用这个校正数与校正后的这个样本的所以reads数求商
(read/基因长度)/[(read/基因1长度)+(read/基因2长度)+...]  

总之是为了消除样品不一致，与基因长度不一致的影响

### 后续R语言分析
代码另附



### snakemake运行代码
```
snakemake -s {*.smk.py} --cores {thread} -pn
```
-pn 表示试运行

流程图
```
snakemake --dag -s rawfq2count.smk.py | dot -Tpdf > dag.pdf  
```
svg改pdf可以生成pdf

### PS
~~被CSDN气哭，一群傻逼下个SRA的方法过时了还抄来抄去~~

