import os
import sys
import argparse
import subprocess



def args_check():
    parser = argparse.ArgumentParser(description= """
     
     
     
    parametres:
    
    --path_fastq : path to the folder fastq
    --path_out : path for saving the result of the star(defaut : /home/storage_1/yuping/matrix_data/ )
    --path_inde : path to the folder references of the mapping
    -- e : threads
    
     """
    ,formatter_class = argparse.RawTextHelpFormatter)
    parser.add_argument("--path_fastq", type=str, help="path to the folder fastq")
    parser.add_argument("--path_out", type=str,default='/home/storage_1/yuping/matrix_data/', help="path for saving the result of the star")
    parser.add_argument("--path_index", type=str,default='/home/storage_1/yuping/index_star_mouseGRC39/',help="path to the folder references of the mapping")
    parser.add_argument("--e", type=int,default=20, help="threads")
    parameters = parser.parse_args()
    return parameters 

def read_dict(your_path):
    dict_files={}
    n=1
    files = os.listdir(your_path)
    for file in files:
        test2=file.split('.')
        test3=test2[0].split('_')
        if len(test3)==1:
            dict_files[n]=test3
            n=n+1
        if len(test3)==2 and test3[1]=='1':
            test4=test3[0]+'_2'
            dict_files[n]=[test2[0],test4]
            n=n+1    
    return dict_files


def star(path_out,your_path,dict_files,path_index,e):
    """
    this fonction will use STAR automatically.
    
    """   
    
    for i in dict_files:
        if len(dict_files[i])==2:
            print('pair_end')
            name=dict_files[i][0].split('_')[0]
            star = f' STAR --runMode alignReads --outSAMtype BAM SortedByCoordinate --genomeDir {path_index} --readFilesIn {your_path}/{dict_files[i][0]}.fastq {your_path}/{dict_files[i][1]}.fastq --runThreadN {e} --outFileNamePrefix {path_out}{name} --quantMode GeneCounts'
            subprocess.call(star, shell=True)
        else:
            print('single_end')
            star = f' STAR --runMode alignReads --outSAMtype BAM SortedByCoordinate --genomeDir {path_index} --readFilesIn {your_path}/{dict_files[i][0]}.fastq --runThreadN {e} --outFileNamePrefix {path_out}{dict_files[i][0]} --quantMode GeneCounts'
            subprocess.call(star, shell=True)


def main():
    parameters = args_check()
    dict_files = read_dict(parameters.path_fastq)
    star(parameters.path_out,parameters.path_fastq,dict_files,parameters.path_index,parameters.e)
    print('------success----------')


if __name__ == "__main__":
    main() 