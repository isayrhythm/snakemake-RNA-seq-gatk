
import os
import threading


outpath = '/home/luotao/forwork/RNA-seq/sra/'


args_task = [] 

with open('/home/luotao/forwork/RNA-seq/sra/SraAccList.csv') as f:
    for i in f:
        args_task.append(i[:-1])



def download(i):
    os.system(f'prefetch {i} --output-directory {outpath}')


        
if __name__ == '__main__':
    for i in args_task:
        task_thread = threading.Thread(target=download, args=(i,))
        task_thread.start()
        print(i)  

