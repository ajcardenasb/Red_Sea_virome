library(phyloseq)
library(microbiome)
library(plyr)

setwd("~/Documents/GitHub/Red_Sea_virome")

#### Data in
dna=read.table("./Input_files/DNA_counts", header = TRUE, row.names = 1)
rna=read.table("./Input_files/RNA_counts", header = TRUE, row.names = 1)
met=read.table("./Input_files/redSea_virome_metadata", header = TRUE, row.names = 2, fill = TRUE)
met$Color=paste("#",met$Color, sep ="")
met$Color=as.factor(met$Color)

#### Identify and remove low-count samples
dna.o1<-dna[, colSums(dna) > 50] # 57/61
dna.o<-dna.o1[ , -which(names(dna.o1) %in% "MG_ACR5")]
rna.o<-rna[, colSums(rna) > 50] # 22/32, >50 leaves 29 samples

##ACR was left with no replicates, therefore ACR5 dna was removed too

## Removing low-count samples
colnames(dna.o)=gsub("MG_", "", colnames(dna.o))
colnames(rna.o)=gsub("MT_", "", colnames(rna.o))


## TRansform the data
require(Seurat)
clr.dna=NormalizeData((dna.o), normalization.method = "CLR")
clr.rna=NormalizeData((rna.o), normalization.method = "CLR")

##Spearman coorrelations
long_dna = melt(clr.dna, value.name = c("MG"))#### 
long_rna = melt(clr.rna, value.name = c("MT"))#### 
long_all=merge(long_dna, long_rna, by.x=c("Var1","Var2"), by.y=c("Var1","Var2"))
cor.test(long_all$MG, long_all$MT, method = "spearman")

SpearCC = cor(otus_log, method = "spearman")#### Used to make final SRCC plots
long_SpearCC = melt(SpearCC, value.name = c("SCC"))#### Melting Data into Long Format
