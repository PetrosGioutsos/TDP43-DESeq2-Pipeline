###Please consider the README file for more information
###regarding the biological experiment and the  
###in silico transcriptomic analysis that followed

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("DESeq2")

library("DESeq2")

if (!require("rstudioapi"))
  install.packages("rstudioapi")

library("rstudioapi")

###set working directory
workingDir<-dirname(rstudioapi::getSourceEditorContext()$path)

outputDir=paste(workingDir,"deseqResults", sep = "/", collapse = "/")
dir.create(outputDir, showWarnings = FALSE)

setwd(workingDir)

###data calling 
#for the htseq counts file
directory="htseq_output/"
filename = "SRR10045016-17-18-19-20-21_counts_filtered.csv"
htseq_counts_filename=paste(directory, filename, sep = "/")

###reading the file
origData <- read.csv(htseq_counts_filename, sep ="\t", header = FALSE)

###columns and rows naming
rownames(origData) <-origData[,1]
countdata = origData[,2:7]
sampleName=c("rep1","rep2","rep3","_rep1","_rep2","_rep3")
colnames(countdata)<-sampleName

thresh <- 100
countdata <- countdata[rowSums(countdata[,-1]) > thresh, ]

condition=as.factor(c("H","H","H","D","D","D")) ## set the header
coldata=data.frame(sampleName,condition)

ddsTable <- DESeqDataSetFromMatrix(countData = countdata, colData = coldata, design = ~ condition )

## basic analysis
dds <- DESeq(ddsTable)
res <- results( dds )
summary(res)


##order by baseMean
results_Exp_Ordered<-res[ order( res$baseMean ,decreasing = TRUE), ] 
results_Exp_Ordered_DF<-data.frame(results_Exp_Ordered)
head(results_Exp_Ordered_DF)

##order by log2FoldChange
results_log2FoldChange_Ordered<-res[ order( res$log2FoldChange ,decreasing = TRUE), ] 
results_log2FoldChange_Ordered_DF<-data.frame(results_log2FoldChange_Ordered)
head(results_log2FoldChange_Ordered_DF)

##order by padj
results_padj_Ordered<-res[ order( res$padj ,decreasing = FALSE), ] 
results_padj_Ordered_DF<-data.frame(results_padj_Ordered)
head(results_padj_Ordered_DF)

##select padj < 0.01  and 
padj.cutoff <- 0.01
results_padj_0.01 <- res[ which(res$padj < padj.cutoff ), ]
results_padj_0.01


##select log2FoldChange>0.58  
log2FoldChange.cutoff <- 0.58  ##corresponds to 1.5 fold change
results_log2FoldChange_0.58 <- res[ which(res$log2FoldChange > log2FoldChange.cutoff ), ]
results_log2FoldChange_0.58


##select padj < 0.01  and  and select log2FoldChange>0.58  (absolute value)
res1<- res[ which(res$padj < padj.cutoff ), ]
res_final <- res1[ which(abs(res1$log2FoldChange) > log2FoldChange.cutoff ), ]
res_final

## write the results to a file
outputFileName="results_DESeq.csv"
outputFullFileName=paste(outputDir,outputFileName,sep="/")
write.table(res_final,outputFullFileName,sep="\t", row.names = TRUE,col.names=NA)##col.name needed to add blank header to cols



deseq2ResDF <- as.data.frame(res)
# Examine this data frame
head(deseq2ResDF)


DESeq2::plotMA(res,ylim=c(-2,2))  
##plot.default(deseq2ResDF$baseMean,deseq2ResDF$log2FoldChange,xlim=c(10,10000),ylim=c(-1,1),log = "x") ##same as plotMA(res)
DESeq2::plotMA(res_final)

#interactive gene identification
#idx <- identify(res$baseMean, res$log2FoldChange) #Uncomment to identify genes
#rownames(res)[idx] 
#idx

plotCounts(dds, gene="ENSG00000130312", intgroup="condition")


rld <- rlog(dds)
plotPCA(rld)

sessionInfo()
