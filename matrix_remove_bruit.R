### package used
library("Rsubread")
library('ggplot2')
library("DESeq2")

### remove differents type RNA keep only mRNA


#### load raw data
matrice_count=read.csv("/Users/daiyuping/Desktop/gitstage/dataselection/creation_of_the_dataset/matrice_count.csv",row.names = 1)

file_mat='/Users/daiyuping/Desktop/gitstage/dataselection/creation_of_the_dataset/matrice_count.csv'
matrice_mRNA=read.csv(file=file_mat,header=TRUE, sep=",", row.names = 1)

#### read annotation of gene

annotation_grf=read.csv("/Users/daiyuping/Desktop/gitstage/dataselection/creation_of_the_dataset/annotation_grf.csv",row.names = 1)


#### remove differents type RNA keep only mRNA
annotation_total=c(row.names(annotation_grf))
annotation_mRNA=c(row.names(annotation_grf)[which(annotation_grf$gene_type=='protein_coding')])
matrice_mRNA=matrice_count[annotation_mRNA,]
annotation_autre=setdiff(annotation_total,annotation_mRNA)
matrice_non_mRNA=matrice_count[annotation_autre,]

### read sample annotation
annotation_sample=read.csv("/Users/daiyuping/Desktop/gitstage/dataselection/creation_of_the_dataset/liver_final.csv", sep = ";",row.names = 'SRR',header = TRUE)
annotation_sample=annotation_sample[c('Sex','age_class','series_id')]

cts <- as.matrix(matrice_mRNA)
cts= cts[,order(names(matrice_mRNA))]

coldata=annotation_sample[order(row.names(annotation_sample)),]

all(rownames(coldata) == colnames(cts))


dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = coldata,
                              design = ~ series_id)


dds <- DESeq(dds)
resMF <- results(dds)

keep <- rowSums(counts(dds)) >= 20							# Parameter: keep only counts >= 20, if keep only counts >=10,removeBatchEffect will return value negetif ***
dds <- dds[keep,]
vsd <- vst(dds, blind=FALSE)

### batch effect


assay(vsd)  <- removeBatchEffect(assay(vsd), 
                                 design = model.matrix(~ coldata$age_class+coldata$Sex), 
                                 batch = coldata$series_id)

which(assay(vsd)==min(assay(vsd)),arr.ind=TRUE) ## make sure all value is positif

matrix_final=assay(vsd)

write.csv(matrix_final,"/Users/daiyuping/Desktop/gitstage/dataselection/creation_of_the_dataset/matrice_final_liver_remove.csv")

### plot with series_id(GEO)
# pcaData <- plotPCA(vsd, intgroup=c("series_id"), returnData=TRUE)
# percentVar <- round(100 * attr(pcaData, "percentVar"))
# ggplot(pcaData, aes(PC1, PC2, color=series_id,label = rownames(coldata))) +
#   geom_point(size=3) +
#   xlab(paste0("PC1: ",percentVar[1],"% variance")) +
#   ylab(paste0("PC2: ",percentVar[2],"% variance")) +
#   coord_fixed()+ggtitle('Principal component plot of the sample-to-sample distances')
# 
# 

### plot with age_class and sex
library("ggplot2")
pcaData <- plotPCA(vsd, intgroup=c("age_class","Sex"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC2, PC1, color=age_class, shape=Sex,label = rownames(coldata))) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed()+ggtitle('Principal component plot of the sample-to-sample distances')



###seperation of the matrix for each class of age
class_age=unique(coldata$age_class)
sex=unique(coldata$Sex)

for (i in class_age) {
  for( j in sex){
    print(i)
    print(j)
    row_name=c(row.names(annotation_sample)[which(annotation_sample$Sex==j)])
    cts=matrix_final[,row_name]
    coldata=annotation_sample[row_name,]
    cts_sex_class=cts[,c(row.names(coldata)[which(coldata$age_class==i)])]
    coldata_sex_class=coldata[c(row.names(coldata)[which(coldata$age_class==i)]),]
    name_of_cts=paste('/Users/daiyuping/Desktop/gitstage/dataselection/creation_of_the_dataset/matrice_liver2/liver_',i,'_',j,'.csv',sep="")
    name_of_coldata=paste('/Users/daiyuping/Desktop/gitstage/dataselection/creation_of_the_dataset/annotation_liver2/annotation_liver_',i,'_',j,'.csv',sep="")
    print(name_of_cts)
    print(name_of_coldata)
    write.csv(cts_sex_class,name_of_cts)
    write.csv(coldata_sex_class,name_of_coldata)
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



