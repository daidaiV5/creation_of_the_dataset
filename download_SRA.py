import os
import sys
import argparse
import subprocess

def check_file_path(path):
    """ This function checks that a file exist and return an absolute
        path towards the file if it exists.
        Parameter:
            - path : a string representing a path.
    """
    true_path = os.path.expanduser(path)
    if os.path.isfile(true_path):
        return os.path.abspath(true_path)
    msg = "The path: {} is not a file or does not exist.\n".format(path)
    raise argparse.ArgumentTypeError(msg)

def args_check():
    parser = argparse.ArgumentParser(description= """
     it's a program for Automating downloads using Python
     Since there are lots of SRA files associated with our samples, it would take a long time to manually run prefetch and fastq-dump for all the files. 
     Install:
     SRA tools version : > 2.10(because of the fasterq-dump)
     
    parametres:
    
    --list_SRA, type=check_file_path, help="list of SRA:SRRXXXXX(.txt)"
    --fastq, type=str,default='fastq', help="folder name"
    --new_path, type=str,default='/home/storage_1/yuping/raw_data/',help="path to save the SRA fill
    
    (echo '/repository/user/main/public/root= new_path ' > $HOME/.ncbi/user-settings.mkfg) before run the script"
    
    --e, type=int,default='20', help="threads"
 
     """
    ,formatter_class = argparse.RawTextHelpFormatter)
    parser.add_argument("--list_SRA", type=check_file_path, help="list of SRA:SRRXXXXX(.txt)")
    parser.add_argument("--fastq", type=str,default='fastq', help="folder name")
    parser.add_argument("--new_path", type=str,default='/home/storage_1/yuping/raw_data/',help="path to save the SRA fill(echo '/repository/user/main/public/root= new_path ' > $HOME/.ncbi/user-settings.mkfg) before run the script")
    parser.add_argument("--e", type=int,default=20, help="threads")
    parameters = parser.parse_args()
    return parameters    

def open_list(fichier):
    """
    This function for traite the fille fasta, including open the fasta,
    Parameter:
    fichier: str,path of fasta
    vitesse: int
    
    output:
    
    """
    dicte=[]
    with open(fichier) as file_one:
        for line in file_one:
            line = line.strip()
            dicte.append(line)
    return dicte
    
def download_tools(sra_numbers,fastq,new_path,e):
    """ Automating downloads using Python
    To automate this process,the fonction will download each SRA file using prefetch and then run fastq-dump. 
    Parameter:
            - sra_number : a list of all sra_number.
            - path : path to save the rawdata of the SRA (dafault : ~/ncbi/ncbi_public/sra)
            - name_of_document: name of the fold(default: fastq)
            - e : int how many thread (default=6)
    """
    for sra_id in sra_numbers:
        print ("Currently downloading: " + sra_id)
        prefetch = "prefetch " + sra_id
        print ("The command used was: " + prefetch)
        subprocess.call(prefetch, shell=True)

    # this will extract the .sra files from above into a folder named 'fastq'
    for sra_id in sra_numbers:
        print ("Generating fastq for: " + sra_id)
        fastq_dump = f'fasterq-dump {new_path}sra/{sra_id}.sra --outdir {fastq} -e {e} '
        print ("The command used was: " + fastq_dump)
        subprocess.call(fastq_dump, shell=True)    


def main():
    parameters = args_check()
    sra_numbers = open_list(parameters.list_SRA)
    download_tools(sra_numbers,parameters.fastq,parameters.new_path,parameters.e)
    print('------success----------')


if __name__ == "__main__":
    main() 
