


sample = ['WT','sma-2','sma-4']


rule vcf:
    input:
        expand("/home/luotao/c_elegans/RNA-seq/SRA_yy/vcf/{SAMPLE}.filter.vcf",SAMPLE=sample)



rule SelectVariants_SNP:
    input:
        "/home/luotao/c_elegans/RNA-seq/SRA_yy/vcf/{SAMPLE}.vcf"
    output:
        "/home/luotao/c_elegans/RNA-seq/SRA_yy/vcf/{SAMPLE}.snp.vcf"
    shell:
        "docker run -it -v /home/luotao/c_elegans/RNA-seq/SRA_yy/vcf:/home/luotao/c_elegans/RNA-seq/SRA_yy/vcf \
        --user $(id -u):$(id -g) \
        broadinstitute/gatk:latest java -jar /gatk/gatk.jar \
        SelectVariants \
        -select-type SNP \
        -V {input} \
        -O {output}"


rule SelectVariants_INDEL:
    input:
        "/home/luotao/c_elegans/RNA-seq/SRA_yy/vcf/{SAMPLE}.vcf"
    output:
        "/home/luotao/c_elegans/RNA-seq/SRA_yy/vcf/{SAMPLE}.indel.vcf"
    shell:
        "docker run -it -v /home/luotao/c_elegans/RNA-seq/SRA_yy/vcf:/home/luotao/c_elegans/RNA-seq/SRA_yy/vcf \
        --user $(id -u):$(id -g) \
        broadinstitute/gatk:latest java -jar /gatk/gatk.jar \
        SelectVariants \
        -select-type INDEL \
        -V {input} \
        -O {output}"
        
        
# 为SNP作硬过滤
rule VariantFiltration_SNP:
    input:
        "/home/luotao/c_elegans/RNA-seq/SRA_yy/vcf/{SAMPLE}.snp.vcf"
    output:
        "/home/luotao/c_elegans/RNA-seq/SRA_yy/vcf/{SAMPLE}.filter.snp.vcf"
    shell:
        "docker run -it -v /home/luotao/c_elegans/RNA-seq/SRA_yy/vcf:/home/luotao/c_elegans/RNA-seq/SRA_yy/vcf \
        --user $(id -u):$(id -g) \
        broadinstitute/gatk:latest java -jar /gatk/gatk.jar \
        VariantFiltration \
        -V {input} \
        --filter-expression 'QD < 2.0 || MQ < 40.0 || FS > 60.0 || SOR > 3.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0' \
        --filter-name 'Filter' \
        -O {output} "

# 为INDEL作硬过滤
rule VariantFiltration_INDEL:
    input:
        "/home/luotao/c_elegans/RNA-seq/SRA_yy/vcf/{SAMPLE}.indel.vcf"
    output:
        "/home/luotao/c_elegans/RNA-seq/SRA_yy/vcf/{SAMPLE}.filter.indel.vcf"
    shell:
        "docker run -it -v /home/luotao/c_elegans/RNA-seq/SRA_yy/vcf:/home/luotao/c_elegans/RNA-seq/SRA_yy/vcf \
        --user $(id -u):$(id -g) \
        broadinstitute/gatk:latest java -jar /gatk/gatk.jar \
        VariantFiltration \
        -V {input} \
        --filter-expression 'QD < 2.0 || FS > 200.0 || SOR > 10.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0' \
        --filter-name 'Filter' \
        -O {output} "

rule MergeVcfs:
    input:
        i1 = "/home/luotao/c_elegans/RNA-seq/SRA_yy/vcf/{SAMPLE}.filter.snp.vcf",
        i2 = "/home/luotao/c_elegans/RNA-seq/SRA_yy/vcf/{SAMPLE}.filter.indel.vcf"
    output:
        "/home/luotao/c_elegans/RNA-seq/SRA_yy/vcf/{SAMPLE}.filter.vcf"
    shell:
        "docker run -it -v /home/luotao/c_elegans/RNA-seq/SRA_yy/vcf:/home/luotao/c_elegans/RNA-seq/SRA_yy/vcf \
        --user $(id -u):$(id -g) \
        broadinstitute/gatk:latest java -jar /gatk/gatk.jar \
        MergeVcfs \
        -I {input.i1} \
        -I {input.i2} \
        -O {output} "
