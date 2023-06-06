import os

import threading


inpath = '/home/luotao/c_elegans/RNA-seq/SRA_yy/'

nalist = os.listdir(inpath)

def sra2fastq(sra_path):
    if os.path.isdir(inpath+sra_path):
        for i in os.listdir(inpath+sra_path):
            sra = inpath+sra_path+os.sep+i
            print(f"fastq-dump --split-3 {sra} --outdir {inpath+sra_path}")
            os.system(f"fastq-dump --split-3 {sra} --outdir {inpath+sra_path}")



if __name__ == '__main__':
    for sra_path in nalist:
        task_thread = threading.Thread(target=sra2fastq, args=(sra_path,))
        task_thread.start()
        print(sra_path)  