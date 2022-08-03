#Script: analise comparativa de RNE vs NE das amostras V3 e V2 sem amostras com muita expressão

library(DESeq2)
library(dplyr)
getwd()

cts = read.delim("cutt_bw1-ReadCount.tab", sep="", header = T)
dim(cts)
cts2 <- select(cts, -c (HFW29_S4_V3,HFW32_S3_V3, HFW01_V2, HFW01_V3, HFW03_V3, HFW03_V2, HFW10_V2, HFW10_V3, HFW24_V3, HFW24_V2, HFW32_S6_V3, HFW32_V2 ))
Name=colnames(cts2)
Time=c(rep("NRE", 2), rep("RE", 2), rep("NRE", 2), rep("RE", 2), rep("NRE", 4), rep("NRE", 5), rep("RE", 4),
       rep("NRE", 4), rep("RE", 4), rep("NRE", 2), rep("RE", 2), rep("NRE", 9), rep("RE", 2), rep("NRE", 10),
       rep("RE", 2), rep("NRE", 2),  rep("RE", 2), rep("NRE", 2), rep("RE", 2), rep("NRE", 2))
coldata=data.frame(Name,Time)
head(coldata,76)
dim(coldata)

log_cts<- log(cts2+1, 10) # This is log base 10 + 1 for "0"
table(rowSums(log_cts>log10(3))>=15)

keep.exprs<- rowSums(log_cts>log10(3))>=15
cts_filt<-cts2[keep.exprs,]
dim(cts_filt)
head(cts_filt,11)

dds <- DESeqDataSetFromMatrix(countData = round(cts_filt), colData = coldata, design = ~ Time)
dds$Time <- factor(dds$Time, levels = c("NRE","RE"))
head(dds)

dds <- DESeq(dds)
res <- results(dds)
head(res)
dim(res)


write.csv(res,"V2xV366.csv", row.names = T)
