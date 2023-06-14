
import os,csv,shutil
import threading


outpath = '/home/luotao/c_elegans/WGS/glycine_max/'

# 下载的是summary
sra_result = '/home/luotao/c_elegans/WGS/glycine_max/sra_result.csv'

args_task = [] 

with open(sra_result) as f:
    reader = csv.reader(f)
    for i in reader:
        # 7是SRA  13是样品名
        if i[7] == 'Sample Accession':
            pass
        else:
            args_task.append([i[7],i[13].replace('(','_').replace(')','_').replace('-','_')])




def download(i):
    if os.path.exists(outpath+i[1]):
        pass
    else:
        os.makedirs(outpath+i[1])
    os.system(f'prefetch {i[0]} --output-directory {outpath+i[1]}')
    print(os.listdir(outpath+i[1]),'###############')
    folder_1 = os.listdir(outpath+i[1])[0]  # 文件夹
    file_2 = os.listdir(outpath+i[1]+os.sep+folder_1)[0]  # 文件夹中的文件
    shutil.move(outpath+i[1]+os.sep+folder_1+os.sep+file_2,outpath+i[1]) # 文件挪到外层
    os.rmdir(outpath+i[1]+os.sep+folder_1)  # 删除文件夹
    os.rename(outpath+i[1]+os.sep+file_2,outpath+i[1]+os.sep+i[1])        # 
    
    # 鲁棒性很低

        
if __name__ == '__main__':
    for i in args_task:
        task_thread = threading.Thread(target=download, args=(i,))
        task_thread.start()
        print(i)  

