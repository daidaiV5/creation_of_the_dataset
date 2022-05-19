### package used
library("Rsubread")
library('ggplot2')
library("DESeq2")
library("limma")

### remove differents type RNA keep only mRNA


#### load raw data
matrice_count=read.csv("/home/storage_1/yuping/git_selection/matrice_count_hypo.csv",row.names = 1)

file_mat='/home/storage_1/yuping/git_selection/matrice_count_hypo.csv'
matrice_count=read.csv(file=file_mat,header=TRUE, sep=",", row.names = 1)

#### read annotation of gene

annotation_grf=read.csv("/home/storage_1/yuping/git_selection/annotation_grf.csv",row.names = 1)


#### remove differents type RNA keep only mRNA
annotation_total=c(row.names(annotation_grf))
annotation_mRNA=c(row.names(annotation_grf)[which(annotation_grf$gene_type=='protein_coding')])
matrice_mRNA=matrice_count[annotation_mRNA,]
annotation_autre=setdiff(annotation_total,annotation_mRNA)
matrice_non_mRNA=matrice_count[annotation_autre,]

### read sample annotation
annotation_sample=read.csv("/home/storage_1/yuping/git_selection/hypocampus_final.csv", sep = ";",row.names = 'SRR',header = TRUE)
annotation_sample=annotation_sample[c('Sex','age_class','series_id')]

cts <- as.matrix(matrice_mRNA)
cts= cts[,order(names(matrice_mRNA))]

coldata=annotation_sample[order(row.names(annotation_sample)),]


all(rownames(coldata) == colnames(cts))

  
dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = coldata,
                              design = ~ Sex+age_class)


dds <- DESeq(dds)
resMF <- results(dds)

keep <- rowSums(counts(dds)) >= 20							# Parameter: keep only counts >= 20, if keep only counts >=10,removeBatchEffect will return value negetif
dds <- dds[keep,]
vsd <- vst(dds, blind=FALSE)
### batch effect


assay(vsd)  <- removeBatchEffect(assay(vsd), 
                                 design = model.matrix(~ coldata$age_class+coldata$Sex), 
                                 batch = coldata$series_id)

which(assay(vsd)==min(assay(vsd)),arr.ind=TRUE) ## make sure all value is positif

matrix_final=assay(vsd)

write.csv(matrix_final,"/home/storage_1/yuping/git_selection/matrice_final_Hypo_remove.csv")
write.csv(coldata,"/home/storage_1/yuping/git_selection/annotation_final_Hypo_remove.csv")

  ### plot with series_id(GEO)
pcaData <- plotPCA(vsd, intgroup=c("series_id"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=series_id,label = rownames(coldata))) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  coord_fixed()+ggtitle('Principal component plot of the sample-to-sample distances')



### plot with age_class and sex
library("ggplot2")
pcaData <- plotPCA(vsd, intgroup=c("age_class","Sex"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=age_class, shape=Sex,label = rownames(coldata))) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed()+ggtitle('Principal component plot of the sample-to-sample distances')



###seperation of the matrix for each class of age
library("do")

class_age=unique(coldata$age_class)
sex=unique(coldata$Sex)

for(j in sex){
  x=strsplit(j,": ")[[1]]
  x=x[2]
  print(x)
}


dir.create('/home/storage_1/yuping/GCN/dataset_test/hypo_mouse/matrice_hypo')
dir.create('/home/storage_1/yuping/GCN/dataset_test/hypo_mouse/annotation_hypo')
for (i in class_age) {
for(j in sex){
    print(i)
    print(j)
    row_name=c(row.names(annotation_sample)[which(annotation_sample$Sex==j)])
    cts=matrix_final[,row_name]
    coldata=annotation_sample[row_name,]
    cts_sex_class=cts[,c(row.names(coldata)[which(coldata$age_class==i)])]
    coldata_sex_class=coldata[c(row.names(coldata)[which(coldata$age_class==i)]),]
    x=strsplit(j,": ")[[1]]
    j=x[2]
    name_of_cts=paste('/home/storage_1/yuping/GCN/dataset_test/hypo_mouse/matrice_hypo/hypo_',i,'_',j,'.csv',sep="")
    name_of_coldata=paste('/home/storage_1/yuping/GCN/dataset_test/hypo_mouse/annotation_hypo/annotation_hypo_',i,'_',j,'.csv',sep="")
    print(name_of_cts)
    print(name_of_coldata)
    write.csv2(cts_sex_class,name_of_cts,sep=',')
    write.csv2(coldata_sex_class,name_of_coldata,sep=',')
    }
}



### annotation of ontology
###BiocManager::install('org.Mm.eg.db',force = TRUE)
###BiocManager::install("clusterProfiler")  

library('org.Mm.eg.db')
library("clusterProfiler")

keytypes(org.Mm.eg.db)

a=c('ENSMUSG00000040660')



AnnotationDbi::select(org.Mm.eg.db, keys = a, columns = c("SYMBOL","GENENAME","ONTOLOGY","GO"), keytype = "ENSEMBL")



