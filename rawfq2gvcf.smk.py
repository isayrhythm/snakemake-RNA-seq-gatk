# snakemake with WGS/RNA-seq
# 2022-5-28

 
# smaple

SAMPLE_INDEX = ['WT1','WT2','WT3','sma_2_rax5__1','sma_2_rax5__2','sma_2_rax5__3','sma_4_rax3__1','sma_4_rax3__2','sma_4_rax3__3']

elegans_gtf = '/home/luotao/c_elegans/RNA-seq/ref/GCF_000002985.6_WBcel235_genomic.gtf'
genome = '/home/luotao/c_elegans/RNA-seq/ref/GCF_000002985.6_WBcel235_genomic.fna'    #参考库不能改


rule gvcf:
    input:
        #expand("/home/luotao/c_elegans/RNA-seq/SRA_yy/{SAMPLE}/clean_data/{SAMPLE}_1.fq",SAMPLE=SAMPLE_INDEX),
        #expand("/home/luotao/c_elegans/RNA-seq/SRA_yy/{SAMPLE}/clean_data/{SAMPLE}_2.fq",SAMPLE=SAMPLE_INDEX),
        #expand("/home/luotao/c_elegans/RNA-seq/SRA_yy/{SAMPLE}/clean_data/fastp.html", SAMPLE=SAMPLE_INDEX),
        #expand("/home/luotao/c_elegans/RNA-seq/SRA_yy/{SAMPLE}/sam/{SAMPLE}.sam",SAMPLE=SAMPLE_INDEX),
        #expand("/home/luotao/c_elegans/RNA-seq/SRA_yy/{SAMPLE}/{SAMPLE}.bam",SAMPLE=SAMPLE_INDEX),
        #expand("/home/luotao/c_elegans/RNA-seq/SRA_yy/{SAMPLE}/{SAMPLE}_sorted.bam",SAMPLE=SAMPLE_INDEX),
        expand("/home/luotao/c_elegans/RNA-seq/SRA_yy/out/{SAMPLE}.gvcf",SAMPLE=SAMPLE_INDEX)
        

# 质控
rule clean_data:
    input:
        i1 = "/home/luotao/c_elegans/RNA-seq/SRA_yy/{SAMPLE}/{SAMPLE}_1.fastq", 
        i2 = "/home/luotao/c_elegans/RNA-seq/SRA_yy/{SAMPLE}/{SAMPLE}_2.fastq"
    output:
        o1 = "/home/luotao/c_elegans/RNA-seq/SRA_yy/{SAMPLE}/clean_data_SNP/{SAMPLE}_1.fq",
        o2 = "/home/luotao/c_elegans/RNA-seq/SRA_yy/{SAMPLE}/clean_data_SNP/{SAMPLE}_2.fq"
    threads: 16
    log:
        "/home/luotao/c_elegans/RNA-seq/SRA_yy/{SAMPLE}/clean_data_SNP/fastp.html"
    
    shell:
        "fastp -i {input.i1} -I {input.i2} -o {output.o1} -O {output.o2} -w {threads} -h {log} --trim_front1 13 --trim_tail1 1 "
        # 切掉前方13个碱基，切掉末尾1个碱基

# fastq2bam
rule BWA_sam2bam:
    input:
        o1 = "/home/luotao/c_elegans/RNA-seq/SRA_yy/{SAMPLE}/clean_data_SNP/{SAMPLE}_1.fq",
        o2 = "/home/luotao/c_elegans/RNA-seq/SRA_yy/{SAMPLE}/clean_data_SNP/{SAMPLE}_2.fq",
        
    params:
        p1 = "{SAMPLE}"
    output:
        "/home/luotao/c_elegans/RNA-seq/SRA_yy/{SAMPLE}/{SAMPLE}_SNP.bam"
    threads: 16
    shell:
        # 单端测序去掉-1 -2 使用-U
        # 参数-S 输入sam文件；参数-b 指定输出的文件为bam；最后重定向写入bam文件
        "bwa mem {genome} {input.o1} {input.o2} -t {threads} -R '@RG\\tID:{params.p1}\\tSM:{params.p1}\\tPL:Illumina' | samtools view -@ {threads} -Sb - > {output} "
        #    | samtools view -Sb - > 19P0126636WES.bam
        # -R 为添加RG信息 -R '@RG\tID:sample1\tSM:sample1\tPL:Illumina'

        
rule bam_sort:
    input:
        "/home/luotao/c_elegans/RNA-seq/SRA_yy/{SAMPLE}/{SAMPLE}_SNP.bam",
    output:
        "/home/luotao/c_elegans/RNA-seq/SRA_yy/{SAMPLE}/{SAMPLE}_SNP_sorted.bam"
    threads: 20
    shell:
        # -@：设置排序和压缩的线程数，默认单线程
        "samtools sort {input} -o {output} -@ {threads}"
        
        
# 标记PCR重复  
rule MarkDuplicates:
    input:
        i1 = "/home/luotao/c_elegans/RNA-seq/SRA_yy/{SAMPLE}/{SAMPLE}_SNP_sorted.bam",
        i2 = "/home/luotao/c_elegans/RNA-seq/SRA_yy/{SAMPLE}"
    output:
        o1 = "/home/luotao/c_elegans/RNA-seq/SRA_yy/{SAMPLE}/{SAMPLE}_SNP.markdup.bam",
        o2 = "/home/luotao/c_elegans/RNA-seq/SRA_yy/{SAMPLE}/{SAMPLE}_SNP.markdup_metrics.txt"
    shell:
        # docker中调用
        # 创建比对索引文件
        # 将会生成.bai文件
        "docker run -it -v {input.i2}:{input.i2} --user $(id -u):$(id -g) broadinstitute/gatk:latest java -jar /gatk/gatk.jar MarkDuplicates \
        -I {input.i1} \
        -O {output.o1} \
        -M {output.o2} \
        ; samtools index {output.o1} "

#docker run -it -v /home/luotao/c_elegans/RNA-seq/SRA_yy/SRR5851337:/home/luotao/c_elegans/RNA-seq/SRA_yy/SRR5851337 broadinstitute/gatk:latest java -jar /gatk/gatk.jar MarkDuplicates         -I /home/luotao/c_elegans/RNA-seq/SRA_yy/SRR5851337/SRR5851337_SNP_sorted.bam         -O /home/luotao/c_elegans/RNA-seq/SRA_yy/SRR5851337/SRR5851337_SNP.markdup.bam         -M /home/luotao/c_elegans/RNA-seq/SRA_yy/SRR5851337/SRR5851337_SNP.markdup_metrics.txt         ; samtools index /home/luotao/c_elegans/RNA-seq/SRA_yy/SRR5851337/SRR5851337_SNP.markdup.bam

# 创建参考文件索引
rule ref_dict:
    input:
        i1 = genome,
        i2 = '/'.join(genome.split('/')[:-1])
        
    output:
        '/home/luotao/c_elegans/RNA-seq/ref/GCF_000002985.6_WBcel235_genomic.dict',
    shell:
        "docker run -it -v {input.i2}:{input.i2} --user $(id -u):$(id -g) broadinstitute/gatk:latest java -jar /gatk/gatk.jar CreateSequenceDictionary \
        -R {input.i1} \
        -O {output}"
        
        
        
rule HaplotypeCaller:
    input:
        i1 = "/home/luotao/c_elegans/RNA-seq/SRA_yy/{SAMPLE}/{SAMPLE}_SNP.markdup.bam",
        i2 = "/home/luotao/c_elegans/RNA-seq/SRA_yy/{SAMPLE}",
        i3 = '/home/luotao/c_elegans/RNA-seq/ref/GCF_000002985.6_WBcel235_genomic.fna',
        i4 = '/'.join(genome.split('/')[:-1])
    output:
        o1 = "/home/luotao/c_elegans/RNA-seq/SRA_yy/out/{SAMPLE}.gvcf",
    shell:
        "docker run -it \
        -v {input.i2}:{input.i2} \
        -v {input.i4}:{input.i4} \
        -v /home/luotao/c_elegans/RNA-seq/SRA_yy/out:/home/luotao/c_elegans/RNA-seq/SRA_yy/out \
        --user $(id -u):$(id -g) broadinstitute/gatk:latest java -jar /gatk/gatk.jar HaplotypeCaller \
        -R {input.i3} \
        -I {input.i1} \
        -O {output.o1} \
        --emit-ref-confidence GVCF"
        
        
        
        