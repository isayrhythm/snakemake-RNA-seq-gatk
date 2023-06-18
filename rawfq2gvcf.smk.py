# snakemake with WGS/RNA-seq
# 2022-5-28


# smaple
import os
smp_list = os.listdir('/home/luotao/c_elegans/WGS/glycine_max')
smp_list = [i for i in smp_list if 'Miss' in i]

SAMPLE_INDEX = ['Missingdata_17k243','Missingdata_17K0263']# smp_list

mark_list = ['2','3','4']

# 执行文件的上层目录

run_path = '/home/luotao/c_elegans/WGS/glycine_max/'


# refrence
elegans_gtf = '/home/luotao/c_elegans/WGS/ref/Glycine_max.Glycine_max_v2.1.56.gff3'
genome = '/home/luotao/c_elegans/WGS/ref/Glycine_max.Glycine_max_v2.1.dna.toplevel.fa'    #参考库不能改


# 创建索引
# 使用dict结尾可以自动识别, 默认识别尾缀替换成.dict

ref_dict = '/home/luotao/c_elegans/WGS/ref/Glycine_max.Glycine_max_v2.1.dna.toplevel.dict'

if os.path.exists(ref_dict):
    pass
else:
    genome_path = '/'.join(genome.split('/')[:-1])
    os.system(f"docker run -it -v {genome_path}:{genome_path} --user $(id -u):$(id -g) broadinstitute/gatk:latest java -jar /gatk/gatk.jar CreateSequenceDictionary \
        -R {genome} \
        -O {ref_dict}")





# rule start




rule gvcf:
    input:
        expand(run_path+"gvcf/{SAMPLE}_{mark}.gvcf",SAMPLE=SAMPLE_INDEX,mark=mark_list)


# 质控
rule clean_data:
    input:
        i1 = run_path+"{SAMPLE}/{SAMPLE}_1.fastq", 
        i2 = run_path+"{SAMPLE}/{SAMPLE}_2.fastq"
    output:
        o1 = run_path+"{SAMPLE}/clean_data/{SAMPLE}_1.fq",
        o2 = run_path+"{SAMPLE}/clean_data/{SAMPLE}_2.fq"
    threads: 16
    log:
        run_path+"{SAMPLE}/clean_data/fastp.html"
    
    shell:
        "fastp -i {input.i1} -I {input.i2} -o {output.o1} -O {output.o2} -w {threads} -h {log} --trim_front1 13"
        # --trim_front1 13 切掉前方13个碱基，--trim_tail1 1 切掉末尾1个碱基

# fastq2sorted_bam
rule BWA_sam2sorted_bam:
    input:
        o1 = run_path+"{SAMPLE}/clean_data/{SAMPLE}_1.fq",
        o2 = run_path+"{SAMPLE}/clean_data/{SAMPLE}_2.fq",
        
    params:
        p1 = "{SAMPLE}"
    output:
        run_path+"{SAMPLE}/{SAMPLE}_sorted.bam"
        
    threads: 40
    shell:
        # 单端测序去掉-1 -2 使用-U
        # 参数-S 输入sam文件；参数-b 指定输出的文件为bam；最后重定向写入bam文件
        # -@：设置排序和压缩的线程数，默认单线程
        #    | samtools view -Sb - > 19P0126636WES.bam
        # -R 为添加RG信息 -R '@RG\tID:sample1\tSM:sample1\tPL:Illumina'
        "bwa mem {genome} {input.o1} {input.o2} -t {threads} -M -R '@RG\\tID:{params.p1}\\tSM:{params.p1}\\tPL:Illumina' | \
        samtools view -@ {threads} -Sb - | \
        samtools sort -@ {threads} -o {output} - "
        
        
# 标记PCR重复  
rule MarkDuplicates:
    input:
        i1 = run_path+"{SAMPLE}/{SAMPLE}_sorted.bam",
        i2 = run_path+"{SAMPLE}"
    output:
        o1 = run_path+"{SAMPLE}/{SAMPLE}.mkdup.bam",
        o2 = run_path+"{SAMPLE}/{SAMPLE}.mkdup_metrics.txt"
    shell:
        # docker中调用
        # 创建比对索引文件
        # 将会生成.bai文件
        # 可以指定内存 –java-options "-Xmx4g"
        "docker run -it -v {input.i2}:{input.i2} --user $(id -u):$(id -g) broadinstitute/gatk:latest java -jar /gatk/gatk.jar MarkDuplicates \
        -I {input.i1} \
        -O {output.o1} \
        -M {output.o2} \
        ; samtools index {output.o1} "

#docker run -it -v /home/luotao/c_elegans/RNA-seq/SRA_yy/SRR5851337:/home/luotao/c_elegans/RNA-seq/SRA_yy/SRR5851337 broadinstitute/gatk:latest java -jar /gatk/gatk.jar MarkDuplicates         -I /home/luotao/c_elegans/RNA-seq/SRA_yy/SRR5851337/SRR5851337_SNP_sorted.bam         -O /home/luotao/c_elegans/RNA-seq/SRA_yy/SRR5851337/SRR5851337_SNP.markdup.bam         -M /home/luotao/c_elegans/RNA-seq/SRA_yy/SRR5851337/SRR5851337_SNP.markdup_metrics.txt         ; samtools index /home/luotao/c_elegans/RNA-seq/SRA_yy/SRR5851337/SRR5851337_SNP.markdup.bam


        
        
        
rule HaplotypeCaller:
    input:
        i1 = run_path+"{SAMPLE}/{SAMPLE}.mkdup.bam",
        i2 = run_path+"{SAMPLE}",
        i3 = genome,
        i4 = '/'.join(genome.split('/')[:-1]),
        
    params:
        p1 = run_path+"gvcf",  # 路径
        p2 = '{mark}'  # 
    output:
        o1 = run_path+"gvcf/{SAMPLE}_{mark}.gvcf",
    shell:
        "docker run -it \
        -v {input.i2}:{input.i2} \
        -v {input.i4}:{input.i4} \
        -v {params.p1}:{params.p1} \
        --user $(id -u):$(id -g) broadinstitute/gatk:latest java -jar /gatk/gatk.jar HaplotypeCaller \
        -R {input.i3} \
        -I {input.i1} \
        -O {output.o1} \
        -L {params.p2} \
        --native-pair-hmm-threads 4 \
        --emit-ref-confidence GVCF"
        
        
        
        