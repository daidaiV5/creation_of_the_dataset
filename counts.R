library("Rsubread")

#count the fastq to matrix
# paired end vs single end

#input:annotation of the sample 

annotation_sample=read.csv("/home/storage_1/yuping/git_selection/hypocampus_final.csv", sep = ";",row.names = 'SRR',header = TRUE)
annotation_sample=annotation_sample[c('end')]
my_data_single_end=c(row.names(annotation_sample)[which(annotation_sample$end=='single_end')])
my_data <- my_data_single_end


dir='/home/storage_1/yuping/index_star_mouseGRC39'
gtffile <- file.path(dir,"gencode.vM28.annotation.gtf")

dir_bam='/home/storage_1/yuping/raw_data/liver_BAM'
list_name=list.files(dir_bam,pattern = '*.bam')

filenames=file.path(dir_bam,my_data)
file.exists(filenames)


fc_paired_end <- featureCounts(files=filenames, 
                    annot.ext=gtffile, 
                    isGTFAnnotationFile=TRUE,
                    isPairedEnd=TRUE,
                    nthreads = 30)


x=fc_paired_end$counts
write.csv(x,"/home/storage_1/yuping/raw_data/fa_paired_end1.csv")



a=read.csv("/home/storage_1/yuping/raw_data/fc_single_end.csv")
b=read.csv("/home/storage_1/yuping/raw_data/fa_paired_end1.csv")
c=read.csv("/home/storage_1/yuping/raw_data/fa_paired_end2.csv")
d=read.csv("/home/storage_1/yuping/raw_data/fa_paired_end3.csv")
df <-cbind(a,b,c,d,x) 


write.csv(df,"/home/storage_1/yuping/raw_data/matrice_count.csv")


