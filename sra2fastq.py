import os

import threading


inpath = '/home/luotao/forwork/RNA-seq/sra/'

nalist = os.listdir(inpath)

def sra2fastq(sra_path):
    for i in os.listdir(inpath+sra_path):
        if '.sra' in i:
            sra = inpath+sra_path+os.sep+i
            os.system(f"fastq-dump --split-3 {sra} --outdir {inpath+sra_path}")


if __name__ == '__main__':
    for sra_path in nalist:
        task_thread = threading.Thread(target=sra2fastq, args=(sra_path,))
        task_thread.start()
        print(sra_path)