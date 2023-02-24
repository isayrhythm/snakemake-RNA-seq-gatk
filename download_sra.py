
import os
import threading


outpath = '/home/luotao/c.elegans/RNA-seq/sra2/'


args_task = [] 

with open('/home/luotao/c.elegans/RNA-seq/sra2/SraAccList.csv') as f:
    for i in f:
        args_task.append(i[:-1])



def download(i):
    os.system(f'prefetch {i} --output-directory {outpath}')


        
if __name__ == '__main__':
    for i in args_task:
        task_thread = threading.Thread(target=download, args=(i,))
        task_thread.start()
        print(i)  

