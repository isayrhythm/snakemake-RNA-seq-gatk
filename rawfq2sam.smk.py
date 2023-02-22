# snakemake with rna-seq
# 2022-5-22


# smaple

SAMPLE_INDEX = ['SRR23330487','SRR23330488','SRR23330489','SRR23330490','SRR23330491','SRR23330492']

elegans_gtf = '/home/luotao/c.elegans/RNA-seq/ref/GCF_000002985.6_WBcel235_genomic.gtf'

rule all:
    input:
        #expand("/home/luotao/c.elegans/RNA-seq/sra/{SAMPLE}/{SAMPLE}_1.fastq",SAMPLE=SAMPLE_INDEX),
        #expand("/home/luotao/c.elegans/RNA-seq/sra/{SAMPLE}/{SAMPLE}_2.fastq",SAMPLE=SAMPLE_INDEX),
        #expand("/home/luotao/c.elegans/RNA-seq/sra/{SAMPLE}/clean_data/{SAMPLE}_1.fq",SAMPLE=SAMPLE_INDEX),
        #expand("/home/luotao/c.elegans/RNA-seq/sra/{SAMPLE}/clean_data/{SAMPLE}_2.fq",SAMPLE=SAMPLE_INDEX),
        #expand("/home/luotao/c.elegans/RNA-seq/sra/{SAMPLE}/clean_data/fastp.html", SAMPLE=SAMPLE_INDEX),
        #expand("/home/luotao/c.elegans/RNA-seq/sra/{SAMPLE}/sam/{SAMPLE}.sam",SAMPLE=SAMPLE_INDEX),
        #expand("/home/luotao/c.elegans/RNA-seq/sra/{SAMPLE}/{SAMPLE}.bam",SAMPLE=SAMPLE_INDEX)
        #expand("/home/luotao/c.elegans/RNA-seq/sra/{SAMPLE}/{SAMPLE}_sorted.bam",SAMPLE=SAMPLE_INDEX),
        expand("/home/luotao/c.elegans/RNA-seq/sra/{SAMPLE}.count",SAMPLE=SAMPLE_INDEX)
        

# 清洗数据
rule clean_data:
    input:
        i1 = "/home/luotao/c.elegans/RNA-seq/sra/{SAMPLE}/{SAMPLE}_1.fastq", 
        i2 = "/home/luotao/c.elegans/RNA-seq/sra/{SAMPLE}/{SAMPLE}_2.fastq"
    output:
        o1 = "/home/luotao/c.elegans/RNA-seq/sra/{SAMPLE}/clean_data/{SAMPLE}_1.fq",
        o2 = "/home/luotao/c.elegans/RNA-seq/sra/{SAMPLE}/clean_data/{SAMPLE}_2.fq"
    threads: 16
    log:
        "/home/luotao/c.elegans/RNA-seq/sra/{SAMPLE}/clean_data/fastp.html"
    
    shell:
        "fastp -i {input.i1} -I {input.i2} -o {output.o1} -O {output.o2} -w {threads} -h {log}"
    
# 回贴
rule align:
    input:
        o1 = "/home/luotao/c.elegans/RNA-seq/sra/{SAMPLE}/clean_data/{SAMPLE}_1.fq",
        o2 = "/home/luotao/c.elegans/RNA-seq/sra/{SAMPLE}/clean_data/{SAMPLE}_2.fq"
    output:
        "/home/luotao/c.elegans/RNA-seq/sra/{SAMPLE}/sam/{SAMPLE}.sam"
    params:
        genome = '/home/luotao/c.elegans/RNA-seq/ref/c.elegans'  #
    threads: 16
    shell:
        #单端测序去掉-1 -2 使用-U
        "hisat2 --time --threads {threads} -x {params.genome} -1 {input.o1} -2 {input.o2} -S {output}"

# BAM 文件是压缩的二进制文件，使得文件压缩比提高了
rule sam2bam:
    input:
        "/home/luotao/c.elegans/RNA-seq/sra/{SAMPLE}/sam/{SAMPLE}.sam"
    output:
        "/home/luotao/c.elegans/RNA-seq/sra/{SAMPLE}/{SAMPLE}.bam"
    threads:20
    shell:
        #参数-S 输入sam文件；参数-b 指定输出的文件为bam；最后重定向写入bam文件
        "samtools view -S {input} -b > {output} -@ {threads}"
# 排序
rule bam_sort:
    input:
        "/home/luotao/c.elegans/RNA-seq/sra/{SAMPLE}/{SAMPLE}.bam",
    output:
        "/home/luotao/c.elegans/RNA-seq/sra/{SAMPLE}/{SAMPLE}_sorted.bam"
    threads: 20
    shell:
        # -@：设置排序和压缩的线程数，默认单线程
        "samtools sort {input} -o {output} -@ {threads}"

#reads计数
rule htseq:
    input:
        "/home/luotao/c.elegans/RNA-seq/sra/{SAMPLE}/{SAMPLE}_sorted.bam"
    output:
        "/home/luotao/c.elegans/RNA-seq/sra/{SAMPLE}.count"
    threads:20  #似乎不能够多线程
    shell:  
        "htseq-count -r name -f bam -n {threads} {input} {elegans_gtf} "






    
    
    
    