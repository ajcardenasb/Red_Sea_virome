library(VennDiagram)
library(gridExtra)
library(data.table)
library(reshape2)

setwd("~/Documents/GitHub/Red_Sea_virome")


### venn
dna=read.table("./Input_files/MG_filtered_counts", header = TRUE)
rna=read.table("./Input_files/MT_filtered_counts", header = TRUE)
vir.info=read.table("./Input_files/all_viral_types.txt", header = FALSE)
dna.l=reshape2::melt(as.matrix(dna), keep.rownames = TRUE)
rna.l=reshape2::melt(as.matrix(rna), keep.rownames = TRUE)

dna[order(rowSums(dna),decreasing = TRUE),][1:20,]
rna[order(rowSums(rna),decreasing = TRUE),][1:20,]

### All samples 
mg=subset(dna.l, value > 0 )
mt=subset(rna.l, value > 0 )
A=venn.diagram(list(mg$Var1, mt$Var1), NULL, category.names = c("Metagenomes", "Metatranscriptomes"), rotation.degree = 90, fill = c("#88419d",  "#02818a"), col = "white", lty="blank", cat.pos = c(180,0), cat.fontfamily = "Arial",  fontfamily="Arial" , alpha = rep(0.7, 2))
grid.arrange(grobTree(A), ncol = 1, nrow = 1)

message("Only found in metagenomes: ")
only.mg=as.data.frame(unique(subset(mg, !mg$Var1 %in% mt$Var1)$Var1)) # Qinviridae     Nyamiviridae   Solinviviridae
only.mg$viralType=vir.info$V2[match(only.mg$`unique(subset(mg, !mg$Var1 %in% mt$Var1)$Var1)`, vir.info$V1)]

message("Only found in metatranscriptomes: ")
unique(subset(mt, !mt$Var1 %in% mg$Var1)$Var1)

###### MYC
mg.myc=subset(mg, value > 0 & Var2 %like% "MYC" )
mt.myc=subset(mt, value > 0 & Var2 %like% "MYC" )
B=venn.diagram(list(mg.myc$Var1, mt.myc$Var1), NULL, category.names = c("MYC MG", "MYC MT"), rotation.degree = 90, fill = c("#88419d",  "#02818a"), col = "white", lty="blank", cat.pos = c(180,0), cat.fontfamily = "Arial",  fontfamily="Arial" , alpha = rep(0.7, 2))
grid.arrange(grobTree(B), ncol = 1, nrow = 1)

message("Only found in metagenomes: ")
unique(subset(mg.myc, !mg.myc$Var1 %in% mt.myc$Var1)$Var1) # Qinviridae     Nyamiviridae   Solinviviridae

message("Only found in metatranscriptomes: ")
unique(subset(mt.myc, !mt.myc$Var1 %in% mg.myc$Var1)$Var1)


### per categories

dna$type=vir.info$V2[match(rownames(dna), vir.info$V1)]
dna.type=aggregate(dna[,1:57], by = list(Type=dna$type), FUN =  sum)
dna.type.perc=sweep(dna.type[,2:58],2,colSums(dna.type[,2:58]),"/")
