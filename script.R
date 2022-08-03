#Script: analise comparativa de RNE vs NE das amostras V3

library(DESeq2)
getwd()

cts = read.delim("cutt_bw1-ReadCount.tab", sep="", header = T)
#rownames(cts)=cts[,1]
#cts=cts[,-1]
dim(cts)



Name=colnames(cts)
Time=c(rep("Control", 3), "Case", rep("Control",6), rep("Case", 4), rep("Control",12), "Case", rep("Control", 17), rep("Case",48))
coldata=data.frame(Name,Time)
head(coldata,14)
dim(coldata)

log_cts<- log(cts+1, 10) # This is log base 10 + 1 for "0"
table(rowSums(log_cts>log10(3))>=5)

keep.exprs<- rowSums(log_cts>log10(3))>=5
cts_filt<-cts[keep.exprs,]
dim(cts_filt)
head(cts_filt,10)

dds <- DESeqDataSetFromMatrix(countData = round(cts_filt), colData = coldata, design = ~ Time)
dds$Time <- factor(dds$Time, levels = c("Control","Case"))
head(dds)

dds <- DESeq(dds)
res <- results(dds)
head(res)
dim(res)

resOrdered_padj <- res[which(res$padj < 0.0001),]
dim(resOrdered_padj) #436 DEGs


write.csv(res,"DEGs_arvc.csv", row.names = T)