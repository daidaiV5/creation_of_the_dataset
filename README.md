# creation_of_the_dataset
## 1. GSM to SRR : 

### package : 
#### pysradb :
version : 0.9.7
https://github.com/saketkc/pysradb


## 2.download the fastq from list of SRA

### Package:
#### SRA tools : 
version : > 2.10(because of the fasterq-dump)
https://github.com/ncbi/sra-tools/wiki/HowTo:-fasterq-dump

impotant :change the default path to which SRA files are downloaded :  https://www.biostars.org/p/159950/

	path : echo \'/repository/user/main/public/root= "{new_path}"\' > $HOME/.ncbi/user-settings.mkfg' 

(new_path wil use for the next step: parameter path)
				
#### download_SRA.py : Script for download the fastq 

Automating downloads using Python ,Since there are lots of SRA files associated with our samples,it would take a long time to manually run prefetch and fastq-dump for all the files. 

To automate this process,the fonction will download each SRA file using prefetch and then run fasterq-dump. 
    
    Parameter:
    
		--list_SRA, type=check_file_path, help="list of SRA:SRRXXXXX(.txt)"
    	--fastq, type=str,default='fastq', help="folder name"
		--new_path, type=str,default='/home/storage_1/yuping/raw_data/',help="path to save the SRA fill (echo '/repository/user/main/public/root= new_path ' > $HOME/.ncbi/user-settings.mkfg) before run the script"
    	--e, type=int,default='20', help="threads"
		
## 3.STAR 

1.creation of the index(reference)

	STAR --runThreadN 20 --runMode genomeGenerate --genomeDir /home/storage_1/yuping/index_star_mouseGRC39/ --genomeFastaFiles /home/storage_1/yuping/index_star_mouseGRC39/GRCm39.primary_assembly.genome.fa --sjdbGTFfile /home/storage_1/yuping/index_star_mouseGRC39/gencode.vM28.annotation.gtf
	
--genomeFastaFiles :GENCODE: files marked with PRI (primary). Strongly recommended for mouse and human:
http://www.gencodegenes.org/. Genome sequence, primary assembly (GRCm39)

--sjdbGTFfile : gencode.vM28.annotation.gtf

2.mapping star
#### star.py : Script for mapping automatically for a lot of fastq

the script will open a folder which contain a lot of fastq, and it will automatically determine whether it is single-end or paired-end, and use the STAR for mapping.

	--path_fastq, type=str, help="path to the folder fastq
	--path_out, type=str,default='/home/storage_1/yuping/matrix_data/', help="path for saving the result of the star
	--path_index, type=str,default='/home/storage_1/yuping/index_star_mouseGRC39/',help="path to the folder references of the mapping
	--e, type=int,default=20, help="threads"


## 4.matrix count


for counting reads or fragments within R/Bioconductor is the Rsubread package which contains the featureCounts function (Liao, Smyth, and Shi 2014). This is very simple to use and very fast, and returns the count matrix as part of the result. 

package use : Rsubread


#### count.R : the script will use for counting read from the fastq, and return the count matrix 
 

## 5.filtering:

#### matrix_remove_bruit.R : the script will use for filtering RNA protein non-coding and remove batch effect





