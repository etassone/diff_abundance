#install DESeq2 only needed on first time using/running DESEq2
#source("http://bioconductor.org/biocLite.R"biobiocLite("DESeq2")

#update bioconductor packages - only needed once in a while
#source("http://bioconductor.org/biocLite.R")
#biocLite(character(), ask=FALSE)


#load required packages
library("DESeq2") #DE testing

library("cummeRbund") #has some clustering algorithms that I like to use
library("Hmisc") #graphics, clustering, data analysis, etc.

library("RColorBrewer") #color palettes
library("gplots") #plotting graphs

#color palette for heatmaps
hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)


#LET'S GO!
#http://dwheelerau.com/2014/02/17/how-to-use-deseq2-to-analyse-rnaseq-data/
###########IMPORTING COUNT DATA FILES (FROM HTSEQ)

directory <- "D:/set/directory/to/sorted/files/"
sampleFiles <- grep("sorted",list.files(directory),value=TRUE) #search for counts files for file input
sampleFiles

#Make sure to include condition for each rep (e.g. if you have 6 hot and 6 cold treatments, you need 6 "HOT", and then 6 "COLD")
sampleCondition <- c("TN","TN","TN","TN","TN","TN","TS","TS","TS","TS","TS","TS")

#This creates a dataframe of your files and the conditions to import into DESeq later
sampleTable <- data.frame(sampleName=sampleFiles, fileName=sampleFiles, condition=sampleCondition)
sampleTable

ddsHTSeq<-DESeqDataSetFromHTSeqCount(sampleTable=sampleTable, directory=directory, design=~condition)
ddsHTSeq

#Levels in colData are used in log calculation, set 'no treatment' first. 
colData(ddsHTSeq)$condition<-factor(colData(ddsHTSeq)$condition, levels=c("TN","TS"))
colData(ddsHTSeq)$condition


######

#ddsHTSeq$condition <- relevel(ddsHTSeq$condition, "Female") #changes base level - this doesn't change testing, just is nice for plots

##RUN DESeq on collapsed replicates using default WALD test
dds<-DESeq(ddsHTSeq, test="Wald")
res<-results(dds)
res<-res[order(res$padj),]
head(res)

res

#order by adjusted p-value
res.Ordered <- res[order(res$padj),] 
summary(res.Ordered) #will give values 

#pull out significant genes at 0.05 level - if too many, reduce padj #
res.sig <- subset(res.Ordered, padj < 0.05) #
summary(res.sig)
head(res.sig)




#separate out lists of upregulated and downregulated
res1UP <- subset(res.sig, log2FoldChange > 1) # genes 59
summary(res1UP)


#write sig results to a csv file: Change values to reflect experiment and date tested
write.csv(as.data.frame(res.sig), file = "File_DATE.csv")

#See if this is helpful!!!! WRITE only UP genes
write.csv(as.data.frame(resUP), file = "res3UP_sig_12AUG15.csv")


##Plot results - DE genes should be in red
#plotMA(dds,alpha=0.05,ylim=c(-2,2),main="DESeq2")
plotMA(dds,ylim=c(-2,2),main="DESeq2")
dev.copy(png,"deseq2_MAplot_alpha05_date.png")
dev.off()


#What are the columns
mcols(res,use.names=TRUE)
write.csv(as.data.frame(res),file="sim_condition_results_deseq2_27May16.csv")


#Transform the data to see if conditions cluster based purely on the datasets


rld <- rlogTransformation(dds, blind=TRUE)
vsd <- varianceStabilizingTransformation(dds, blind=TRUE)


##Lets make a heatmap 
library("RColorBrewer")
library("gplots")
select <- order(rowMeans(counts(dds,normalized=TRUE)),decreasing=TRUE)[1:60]
hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)
heatmap.2(counts(dds,normalized=TRUE)[select,], col = hmcol,
          Rowv = FALSE, Colv = FALSE, scale="none",
          dendrogram="none", trace="none", margin=c(10,6))
dev.copy(png,"DESeq2_heatmap1_27MAY16.png")
dev.off()
heatmap.2(assay(rld)[select,], col = hmcol,
          Rowv = FALSE, Colv = FALSE, scale="none",
          dendrogram="none", trace="none", margin=c(10, 6))
dev.copy(png,"DESeq2_rld_heatmap2_12AUG15")
dev.off()

heatmap.2(assay(vsd)[select,], col = hmcol,
          Rowv = FALSE, Colv = FALSE, scale="none",
          dendrogram="none", trace="none", margin=c(10, 6))
dev.copy(png,"DESeq2_heatmap3_27MAY16")
dev.off()

###Cluster to make a nice dendrogram

distsRL <- dist(t(assay(rld)))
mat <- as.matrix(distsRL)
rownames(mat) <- colnames(mat) <- with(colData(dds),
  paste(condition,sampleFiles , sep=" : "))

library("RColorBrewer")
library("gplots")
select <- order(rowMeans(counts(dds,normalized=TRUE)),decreasing=TRUE)[1:30]
hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)
heatmap.2(counts(dds,normalized=TRUE)[select,], col = hmcol,
          Rowv = FALSE, Colv = FALSE, scale="none",
          dendrogram="none", trace="none", margin=c(10,6))
dev.copy(png,"DESeq2_heatmap1.png")
dev.off()
heatmap.2(assay(rld)[select,], col = hmcol,
          Rowv = FALSE, Colv = FALSE, scale="none",
          dendrogram="none", trace="none", margin=c(10, 6))
dev.copy(png,"DESeq2_heatmap2.png")
dev.off()
heatmap.2(assay(vsd)[select,], col = hmcol,
          Rowv = FALSE, Colv = FALSE, scale="none",
          dendrogram="none", trace="none", margin=c(10, 6))
dev.copy(png,"DESeq2_heatmap3.png")
dev.off()


#updated in latest vignette (See comment by Michael Love)
#this line was incorrect
#heatmap.2(mat, trace="none", col = rev(hmcol), margin=c(16, 16))
#From the Apr 2015 vignette
hc <- hclust(distsRL)
heatmap.2(mat, Rowv=as.dendrogram(hc),
          symm=TRUE, trace="none",
          col = rev(hmcol), margin=c(13, 13))
dev.copy(png,"deseq2_heatmaps_samplebysample_12AUG15.png")
dev.off()

#Cluster based on PCA
print(plotPCA(rld, intgroup=c("condition")))
dev.copy(png,"deseq2_pca.png")
dev.off()



