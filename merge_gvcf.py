

import os,csv,shutil
import threading

# smaple

# sma-2(rax5) replicate
SAMPLE_INDEX1 = ['sma_2_rax5__1','sma_2_rax5__2','sma_2_rax5__3','sma-2']

# sma-4(rax3) replicate
SAMPLE_INDEX2 = ['sma_4_rax3__1','sma_4_rax3__2','sma_4_rax3__3','sma-4']

# Wild type replicate
SAMPLE_INDEX3 = ['WT1','WT2','WT3','WT']


elegans_gtf = '/home/luotao/c_elegans/RNA-seq/ref/GCF_000002985.6_WBcel235_genomic.gtf'
genome = '/home/luotao/c_elegans/RNA-seq/ref/GCF_000002985.6_WBcel235_genomic.fna'    #参考库不能改



inpath = '/home/luotao/c_elegans/RNA-seq/SRA_yy/out/'



def gvcf2vcf(sample_list,name):
    cmd = f"docker run -it \
        -v {inpath}:{inpath} \
        -v /home/luotao/c_elegans/RNA-seq/ref:/home/luotao/c_elegans/RNA-seq/ref \
        -v /home/luotao/c_elegans/RNA-seq/SRA_yy:/home/luotao/c_elegans/RNA-seq/SRA_yy \
        --user $(id -u):$(id -g) broadinstitute/gatk:latest java -jar /gatk/gatk.jar CombineGVCFs \
        -R {genome} \
        -V {inpath+sample_list[0]}.gvcf \
        -V {inpath+sample_list[1]}.gvcf \
        -V {inpath+sample_list[2]}.gvcf \
        -O /home/luotao/c_elegans/RNA-seq/SRA_yy/{name[0]}.gvcf "
    print(cmd)
    os.system(cmd)
    
    cmd2 = f"docker run -it \
        -v {inpath}:{inpath} \
        -v /home/luotao/c_elegans/RNA-seq/ref:/home/luotao/c_elegans/RNA-seq/ref \
        -v /home/luotao/c_elegans/RNA-seq/SRA_yy:/home/luotao/c_elegans/RNA-seq/SRA_yy \
        --user $(id -u):$(id -g) broadinstitute/gatk:latest java -jar /gatk/gatk.jar GenotypeGVCFs \
        -R {genome} \
        -V /home/luotao/c_elegans/RNA-seq/SRA_yy/{name[0]}.gvcf \
        -O /home/luotao/c_elegans/RNA-seq/SRA_yy/{name[0]}.vcf"
    os.system(cmd2)

if __name__ == '__main__':
    for i in [SAMPLE_INDEX1,SAMPLE_INDEX2,SAMPLE_INDEX3]:
        task_thread = threading.Thread(target=gvcf2vcf, args=(i[:3],i[-1:]))
        task_thread.start()
        print(i)
        

