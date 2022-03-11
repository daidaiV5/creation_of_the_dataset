# creation_of_the_dataset
1. GSM to SRR : 

2.download the 
install package:
1.SRA tools :
	version : > 2.10(because of the fasterq-dump)
	https://github.com/ncbi/sra-tools/wiki/HowTo:-fasterq-dump
		
Script for download the fastq 
Automating downloads using Python ,Since there are lots of SRA files associated with our samples,it would take a long time to manually run prefetch and fastq-dump for all the files. 
To automate this process,the fonction will download each SRA file using prefetch and then run fasterq-dump. 
    Parameter:
            - sra_number : a list of all sra_number.
            - path : path to save the rawdata of the SRA (dafaut : ~/ncbi/ncbi_public/sra)
            - name_of_document: name of the fold(defaut: fastq)
  
