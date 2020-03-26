library(data.table)

setwd("~/Documents/GitHub/Red_Sea_virome/")

dna=read.table("./Input_files/DNA_counts", header = TRUE, row.names = 1)
rna=read.table("./Input_files/RNA_counts", header = TRUE, row.names = 1)
met=read.table("./Input_files/redSea_virome_metadata", header = TRUE)
all=merge(dna, rna, by = "row.names", all = TRUE)
all[is.na(all)] = 0
all.n=apply(all[2:ncol(all)], 2, as.numeric)
rownames(all.n)=all$Row.names

all.2=all.n[, colSums(all.n) > 50]
all.3=all.2[, -which(colnames(all.2) %in% "MG_ACR5")]
write.table(all.3, "./Input_files/MG_MT_filtered_counts", quote =FALSE)


dna.filtered=all.3[, which(colnames(all.3) %like% "MG")]
colnames(dna.filtered)=gsub("MG_", "", colnames(dna.filtered))
write.table(dna.filtered, "./Input_files/MG_filtered_counts", quote =FALSE)

rna.filtered=all.2[, which(colnames(all.2) %like% "MT")]
colnames(rna.filtered)=gsub("MT_", "", colnames(rna.filtered))
write.table(rna.filtered, "./Input_files/MT_filtered_counts", quote =FALSE)


### export log10 normalized data
dna=read.table("./Input_files/MG_filtered_counts", header = TRUE, row.names = 1)
rna=read.table("./Input_files/MT_filtered_counts", header = TRUE, row.names = 1)
met=read.table("./Input_files/redSea_virome_metadata", header = TRUE, row.names = 1, fill = TRUE)

# cretae phyloseq objects
otu.dna= otu_table(dna, taxa_are_rows=TRUE)
sam.dna= sample_data(met)
phy.dna= phyloseq(otu.dna, sam.dna)

otu.rna= otu_table(rna, taxa_are_rows=TRUE)
sam.rna= sample_data(met)
phy.rna= phyloseq(otu.rna, sam.rna)

####Transforming counts 
dna.t=microbiome::transform(phy.dna, transform = "log10p", target = "OTU", shift = 0, scale = 1)
rna.t=microbiome::transform(phy.rna, transform = "log10p", target = "OTU", shift = 0, scale = 1)

log10.dna=otu_table(dna.t)
log10.rna=otu_table(rna.t)

write.table(log10.dna, "./MG_log10_counts", quote =FALSE)
write.table(log10.rna, "./MT_log10_counts", quote =FALSE)

