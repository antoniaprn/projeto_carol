#Script: analise comparativa de RNE vs NE das amostras V3 sem as amostras de alta expressão

library(DESeq2)
library(dplyr)

getwd()

cts = read.delim("cutt_bw1-ReadCount.tab", sep="", header = T)
dim(cts)
cts2 <- select(cts, -c (HFW29_S4_V3,HFW32_S3_V3, HFW01_V2,HFW02_V2, HFW03_V2,
                        HFW04_V2, HFW05_V2, HFW07_V2, HFW08_V2, HFW09_V2,
                        HFW10_V2,HFW11_V2,HFW12_V2, HFW15_V2, HFW16_V2, HFW17_V2,
                        HFW18_V2, HFW21_V2, HFW23_V2,HFW23_V2, HFW25_V2, HFW30_V2,
                        HFW34_V2, HFW35_V2, HFW38_V2, HFW39_V2, HFW40_V2, HFW50_V2,
                        HFW18_V2, HFW24_V2, HFW26_V2, HFW27_V2, HFW144_V2, HFW31_V2,
                        HFW32_V2, HFW36_V2, HFW39_V2, HFW41_V2, HFW43_V2, HFW45_V2, 
                        HFW48_V2, HFW49_V2, HFW51_V2, HFW03_V3, HFW10_V3, HFW32_S6_V3 ))
Name=colnames(cts2)
Time=c(rep("NRE", 2), "RE", "NRE", "RE", rep("NRE", 2), rep("NRE", 2), rep("RE", 2),
       rep("NRE", 2), rep("RE", 2), rep("NRE", 2), "RE", rep("NRE", 3), "RE", rep("NRE", 2),
      rep("NRE", 5), "RE", "NRE", "RE", "NRE", "RE", "NRE")
coldata=data.frame(Name,Time)
head(coldata,38)
dim(coldata)

log_cts<- log(cts2+1, 10) # This is log base 10 + 1 for "0"
table(rowSums(log_cts>log10(3))>=10)

keep.exprs<- rowSums(log_cts>log10(3))>=10
cts_filt<-cts2[keep.exprs,]
dim(cts_filt)
head(cts_filt,12)

dds <- DESeqDataSetFromMatrix(countData = round(cts_filt), colData = coldata, design = ~ Time)
dds$Time <- factor(dds$Time, levels = c("NRE","RE"))
head(dds)

dds <- DESeq(dds)
res <- results(dds)
head(res)
dim(res)


write.csv(res,"v3_RExNRE35.csv", row.names = T)
