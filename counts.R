library("Rsubread")
library('ggplot2')


#

my_data <- readLines("/home/storage_1/yuping/git_selection/paired_end_list.txt")


dir='/home/storage_1/yuping/index_star_mouseGRC39'
gtffile <- file.path(dir,"gencode.vM28.annotation.gtf")

dir_bam='/home/storage_1/yuping/raw_data/liver_BAM'
list_name=list.files(dir_bam,pattern = '*.bam')


my_data=my_data[120:132]

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


matrice_count=read.csv("/Users/daiyuping/Desktop/gitstage/dataselection/creation_of_the_dataset/matrice_count.csv",row.names = 1)

file_mat='/Users/daiyuping/Desktop/gitstage/dataselection/creation_of_the_dataset/matrice_count.csv'
matrice_mRNA=read.csv(file=file_mat,header=TRUE, sep=",", row.names = 1)
### read annotation of gene

annotation_grf=read.csv("/Users/daiyuping/Desktop/gitstage/dataselection/creation_of_the_dataset/annotation_grf.csv",row.names = 1)

annotation_total=c(row.names(annotation_grf))
annotation_mRNA=c(row.names(annotation_grf)[which(annotation_grf$gene_type=='protein_coding')])
matrice_mRNA=matrice_count[annotation_mRNA,]
annotation_autre=setdiff(annotation_total,annotation_mRNA)
matrice_non_mRNA=matrice_count[annotation_autre,]

### read sample annotation
annotation_sample=read.csv("/Users/daiyuping/Desktop/gitstage/dataselection/creation_of_the_dataset/liver_final2.csv", sep = ";",row.names = 'SRR',header = TRUE)
annotation_sample=annotation_sample[c('Sex','age_class','series_id')]
annotation_sample['SRR17261652',]


### Deseq2
library("DESeq2")
cts <- as.matrix(matrice_mRNA)
cts
names(matrice_mRNA)
cts= cts[,order(names(matrice_mRNA))]

coldata=annotation_sample[order(row.names(annotation_sample)),]


all(rownames(coldata) == colnames(cts))

dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = coldata,
                              design = ~ age_class+Sex+series_id)

dds <- DESeq(dds)
resMF <- results(dds)

keep <- rowSums(counts(dds)) >= 15							# Parameter: keep only counts >= 10 ***
dds <- dds[keep,]

vsd_nobatch <- removeBatchEffect(assay(vsd), 
                                 design = model.matrix(~ set2@phenoData@data$genotype), 
                                 covariates = set2@phenoData@data$W_1)

library("pheatmap")
select <- order(rowMeans(counts(dds,normalized=TRUE)),
                decreasing=TRUE)[1:20]
df <- as.data.frame(colData(dds)[,c("age_class","Sex")])
pheatmap(assay(ntd)[select,], cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=FALSE, annotation_col=df)


### Heatmap of the sample-to-sample distances
vsd <- vst(dds, blind=FALSE)
sampleDists <- dist(t(assay(vsd)))


library("RColorBrewer")
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd$Sex, vsd$age_class, sep="-")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)


plotPCA(vsd, intgroup=c("age_class","Sex"))


pcaData <- plotPCA(vsd, intgroup=c("age_class","Sex"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
p=ggplot(pcaData, aes(PC1, PC2, color=age_class, shape=Sex,label = rownames(coldata))) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed()+ggtitle('Principal component plot of the sample-to-sample distances')

ggsave(p,filename = "ACP_male_sex.pdf",width = 24,height = 18,)


### only male

row_name_male=c(row.names(coldata)[which(coldata$Sex=='Male')])
cts=matrice_mRNA[,row_name_male]
coldata=coldata[row_name_male,]


all(rownames(coldata) == colnames(cts))

dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = coldata,
                              design = ~age_class)
dds <- DESeq(dds)
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

res <- results(dds)
summary(res)




'''
it can also be useful to examine the counts of reads for a single gene across the groups. 
A simple function for making this plot is plotCounts, which normalizes counts by the estimated size factors 
(or normalization factors if these were used) and adds a pseudocount of 1/2 to allow for log scale plotting. 
The counts are grouped by the variables in intgroup, where more than one variable can be specified. 
ere we specify the gene which had the smallest p value from the results table created above. 
You can select the gene to plot by rowname or by numeric index.
'''

plotCounts(dds, gene=which.min(res$padj), intgroup="age_class")

d <- plotCounts(dds, gene=which.min(res$padj), intgroup="age_class", 
                returnData=TRUE)
library("ggplot2")
ggplot(d, aes(x=age_class, y=count)) + 
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10(breaks=c(25,100,400))




head(dispersions(dds))

vsd <- vst(dds, blind=FALSE)
sampleDists <- dist(t(assay(vsd)))

sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd$age_class)
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)


plotPCA(vsd, intgroup=c("age_class"))
###seperation of the matrix

row_name_Female=c(row.names(annotation_sample)[which(annotation_sample$Sex=='Female')])
cts=matrice_mRNA[,row_name_Female]
coldata=annotation_sample[row_name_Female,]

cts_male_class=cts[,c(row.names(coldata)[which(coldata$age_class=='class4')])]
coldata_male_class=coldata[c(row.names(coldata)[which(coldata$age_class=='class4')]),]

write.csv(cts_male_class,"/Users/daiyuping/Desktop/gitstage/dataselection/creation_of_the_dataset/matrice_liver/liver_class4_Female.csv")
write.csv(coldata_male_class,"/Users/daiyuping/Desktop/gitstage/dataselection/creation_of_the_dataset/annotation_liver/annotation_liver_class4_Female.csv")


### annotation of ontology
###BiocManager::install('org.Mm.eg.db',force = TRUE)
###BiocManager::install("clusterProfiler")

library('org.Mm.eg.db')
library("clusterProfiler")

keytypes(org.Mm.eg.db)

a=c('ENSMUSG00000040660')



AnnotationDbi::select(org.Mm.eg.db, keys = a, columns = c("SYMBOL","GENENAME","ONTOLOGY","GO"), keytype = "ENSEMBL")





