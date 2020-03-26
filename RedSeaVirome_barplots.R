library(scales)
library(reshape2)
library(ggplot2)
library(vegan)
library(gridExtra)
library(ggdendro)
library(Seurat)
library(plyr)
library(ggpubr)

setwd("~/Documents/GitHub/Red_Sea_virome")

all=read.table("./Input_files/MG_MT_filtered_counts", header = TRUE)
met=read.table("./Input_files/redSea_virome_metadata", header = TRUE)


## replicated barplots top 20 families 

test=as.data.frame(rowSums(all))
topFamilies=all[order(rowSums(all),decreasing = TRUE),][1:20,]
fam.top=subset(all, rownames(all) %in% rownames(topFamilies)) 
fam.bot=subset(all, !rownames(all) %in% rownames(topFamilies)) 
others=as.data.frame(t(colSums(fam.bot)))
rownames(others)="zOthers"
all.3 =rbind(fam.top, others)
all.l=melt(as.matrix(all.3), keep.rownames = TRUE)
colnames(all.l)=c("Family", "Host", "Abundance")
all.l$sample=ifelse(grepl("MG_",all.l$Host), "Metagenome", "Metatranscriptome")
all.l$sample=as.factor(all.l$sample)
all.l$Host=gsub(".*_", "",all.l$Host )
all.l1=subset(all.l, !Host == "ACR5")

P21=c("#D68D93","#ffb5e8", "#fcc2ff", "#ecd4ff", "#d5aaff", "#b28dff","#dcd3ff", "#b5b9ff", "#97a2ff", "#ace7ff", "#85E3FF", "#94E6F9", "#A3EAF3" ,"#B2EDED", "#C2F1E8" ,"#D1F4E2", "#E0F8DC" ,"#EFFBD6" ,"#FFFFD1","#EFC497", "#EFAD97", "#C4836D", "#B76652", "#A59089" )
# sup1=ggplot() +geom_bar(aes(y = Abundance, x = Host, fill = Family), data = all.l1, stat="identity", position = "fill") +  scale_y_continuous(labels = percent_format(), expand = c(0, 0))  + labs( y= "Percentage of sequences", x="Host species") + scale_fill_manual(values=P21) + facet_grid(sample~. ) + theme_bw()  + theme(axis.text.x=element_text(angle=90,hjust=1), plot.title = element_text(hjust = 0.5), strip.background =element_rect(fill="white"))
# svg("../outputs/RedSeaVirome_replicated_barplots.svg",  width = 10, height = 5, pointsize = 12) 
# plot(sup1)
# dev.off()


## mean barplots top 20 families 

all.l1$Host=gsub("[0-9]", "", all.l1$Host)
mg.l=subset(all.l1, sample == "Metagenome")
mg.l$Host=factor(mg.l$Host, levels = c("PAV", "PLE", "FUN", "POR", "STY", "POC", "PAC","GAL", "MIL",  "XEN", "ACA", "DIP", "MYC"))

mt.l=subset(all.l1, sample == "Metatranscriptome")
mt.l$Host=factor(mt.l$Host, levels = c("POC", "ACR", "STY", "GAL", "PAC", "POR", "DIP", "MYC"))


mg.bar=ggplot() +geom_bar(aes(y = Abundance, x = Host, fill = Family), data = mg.l, stat="identity", position = "fill") +  scale_y_continuous(labels = percent_format(), expand = c(0, 0))  + labs( y= "", x="") + scale_fill_manual(values=P21)  + theme_bw()  + theme(axis.text.x=element_text(angle=90,hjust=1), plot.title = element_text(hjust = 0.5), strip.background =element_rect(fill="white"), text = element_text(size=12),  legend.position="none") + coord_flip()
mt.bar=ggplot() +geom_bar(aes(y = Abundance, x = Host, fill = Family), data = mt.l, stat="identity", position = "fill") +  scale_y_continuous(labels = percent_format(), expand = c(0, 0))  + labs( y= "Percentage of sequences", x="") + scale_fill_manual(values=P21)  + theme_bw()  + theme(axis.text.x=element_text(angle=90,hjust=1), plot.title = element_text(hjust = 0.5), strip.background =element_rect(fill="white"), text = element_text(size=12),  legend.position="none") + coord_flip()


####clustering

dna=read.table("./Input_files/MG_filtered_counts", header = TRUE)
clr.dna=as.data.frame(NormalizeData((dna), normalization.method = "CLR"))
colnames(clr.dna)=gsub("MG_", "", colnames(clr.dna))
colnames(clr.dna)=gsub("[0-9]", "", colnames(clr.dna))
dna.mean= as.data.frame(sapply(unique(names(clr.dna)), function(col) rowMeans(clr.dna[names(clr.dna) == col])))
dna.dist= vegdist(t(dna.mean), method  = "bray") 
dna.hclu = as.dendrogram(hclust(dna.dist, method = "ward.D"))
mg.dendo=ggdendrogram(dna.hclu,  rotate = FALSE) + coord_flip() + scale_y_reverse(expand=c(0.2, 0))  + labs( y= "Bray Curtis dissimilarity", x="")

rna=read.table("./Input_files/MT_filtered_counts", header = TRUE)
clr.rna=as.data.frame(NormalizeData((rna), normalization.method = "CLR"))
colnames(clr.rna)=gsub("MT_", "", colnames(clr.rna))
colnames(clr.rna)=gsub("[0-9]", "", colnames(clr.rna))
rna.mean= as.data.frame(sapply(unique(names(clr.rna)), function(col) rowMeans(clr.rna[names(clr.rna) == col])))
rna.dist= vegdist(t(rna.mean), method  = "bray") 
rna.hclu = as.dendrogram(hclust(rna.dist, method = "ward.D"))
mt.dendo=ggdendrogram(rna.hclu,  rotate = FALSE) + coord_flip() + scale_y_reverse(expand=c(0.2, 0))  + labs(caption =  "Bray Curtis dissimilarity", x="")

## boxplots alpha
COLdna=c("#CCEBC5", "#FFFFB3", "#BEBADA", "#FB8072", "#9fb9bf", "#8DD3C7" ,"#B3DE69", "#80B1D3", "#D9D9D9", "#b88c8c", "#FDB462", "#BC80BD", "#FCCDE5") 
COLrna=c( "#8DD3C7", "#FFED6F", "#9fb9bf", "#80B1D3", "#B3DE69", "#FB8072", "#BC80BD", "#FCCDE5") 

#clr.dna.veg=as.data.frame(NormalizeData((dna), normalization.method = "CLR"))
#clr.dna.veg.t=as.data.frame(t(clr.dna.veg))
dna.t=t(dna)
dna.alpha=specpool(dna.t[,1:117], rownames(dna.t))
dna.alpha$host=gsub("[0-9]", "", rownames(dna.t))
dna.alpha$host=gsub("MG_", "", dna.alpha$host)
dna.alpha$host=factor(dna.alpha$host,  levels = c("PAV", "PLE", "FUN", "POR", "STY", "POC", "PAC","GAL", "MIL",  "XEN", "ACA", "DIP", "MYC"))
mg.alpha=ggplot(dna.alpha,aes(x=host, y=chao, fill=host)) +  geom_boxplot() + geom_jitter(width=0.25, alpha=0.5) + labs( y= "", x="" ) + theme(legend.title = element_blank(), axis.text.x = element_text(angle = 90, hjust = 1), panel.grid.major = element_line(colour = "transparent"), panel.grid.minor = element_line(colour = "transparent"), panel.background = element_blank(), axis.line = element_line(colour = "black"),  legend.position="none") + ggtitle("") + scale_fill_manual(values=COLdna) + coord_flip()  + ylim(1, 60) 
#clr.rna.veg=as.data.frame(NormalizeData((rna), normalization.method = "CLR"))
#clr.rna.veg.t=as.data.frame(t(clr.rna.veg))
rna.t=t(rna)
rna.alpha=specpool(rna.t[,1:117], rownames(rna.t))
rna.alpha$host=gsub("[0-9]", "", rownames(rna.t))
rna.alpha$host=gsub("MT_", "", rna.alpha$host)
rna.alpha$host=factor(rna.alpha$host, levels = c("POC", "ACR", "STY", "GAL", "PAC", "POR", "DIP", "MYC"))
mt.alpha=ggplot(rna.alpha,aes(x=host, y=chao, fill=host)) +  geom_boxplot() + geom_jitter(width=0.25, alpha=0.5) + labs( y= "", x="" ) + theme(legend.title = element_blank(), axis.text.x = element_text(angle = 90, hjust = 1), panel.grid.major = element_line(colour = "transparent"), panel.grid.minor = element_line(colour = "transparent"), panel.background = element_blank(), axis.line = element_line(colour = "black"),  legend.position="none") + ggtitle("") + scale_fill_manual(values=COLrna) + coord_flip() 

grid.all=grid.arrange( mg.alpha,  mt.alpha)

###comparison between MG and MT 
all.t=t(all)
all.alpha=specpool(all.t[,1:117], rownames(all.t))
all.alpha$groups=gsub("[0-9]","",rownames(all.alpha))
shapiro.test(all.alpha$chao)
pairwise.t.test(all.alpha$chao,all.alpha$groups, p.adj = "fdr")

t.test(subset(all.alpha, rownames(all.alpha) %like% "POC")$chao ~ subset(all.alpha, rownames(all.alpha) %like% "POC")$groups)
t.test(subset(all.alpha, rownames(all.alpha) %like% "POR")$chao ~ subset(all.alpha, rownames(all.alpha) %like% "POR")$groups)
t.test(subset(all.alpha, rownames(all.alpha) %like% "STY")$chao ~ subset(all.alpha, rownames(all.alpha) %like% "STY")$groups)

###leyend

names.bar=ggplot() +geom_bar(aes(y = Abundance, x = Host, fill = Family), data = mt.l, stat="identity", position = "fill") +  scale_y_continuous(labels = percent_format(), expand = c(0, 0))  + labs( y= "Percentage of sequences", x="") + scale_fill_manual(values=P21)  + theme_bw()  + theme(axis.text.x=element_text(angle=90,hjust=1), plot.title = element_text(hjust = 0.5), strip.background =element_rect(fill="white"), text = element_text(size=12),  legend.position="bottom") + coord_flip()
leg.bar=get_legend(names.bar)

names.box=ggplot(dna.alpha,aes(x=host, y=chao, fill=host)) +  geom_boxplot() + geom_jitter(width=0.25, alpha=0.5) + labs( y= "", x="" ) + theme(legend.title = element_blank(), axis.text.x = element_text(angle = 90, hjust = 1), panel.grid.major = element_line(colour = "transparent"), panel.grid.minor = element_line(colour = "transparent"), panel.background = element_blank(), axis.line = element_line(colour = "black"),  legend.position="bottom") + ggtitle("") + scale_fill_manual(values=COLdna) + coord_flip()  
leg.box=get_legend(names.box)

pdf("./outputs/RedSeaVirome_dend-barpl-alpha.pdf",  width = 20, height = 10, pointsize = 12) 
grid.all=grid.arrange(mg.dendo, mg.bar, mg.alpha, mt.dendo, mt.bar, mt.alpha, as_ggplot(leg.bar),as_ggplot(leg.box), ncol=3, nrow =3)
dev.off()
grid.all=grid.arrange( as_ggplot(leg.bar), ncol=1, nrow =1)


### viral types
vir=read.table("./Input_files/Viral_families_metadata.txt", header = TRUE, row.names = 1, sep = "\t")
all.l1$Type=vir$Presumed[match(all.l1$Family, rownames(vir))]
all.l1$Type=factor(all.l1$Type,  levels = c("Bacteriophage", "Unicellular eukaryotes", "Invertebrates", "Vertebrates"))
COLtype1=c("#BC4B51", "#8FB339","#5B8E7D", "#084C61") 
all.l1$Type2=vir$Type[match(all.l1$Family, rownames(vir))]
all.l1$Type2=factor(all.l1$Type2,  levels = c("dsDNA", "ssRNA", "ssDNA", "dsRNA"))
COLtype2=c("#7E4E60", "#73937E","#85BDBF", "#C6A15B" ) 
vyt.typ.l=subset(all.l1, !all.l1$Type == "NA")

type1=ggplot() +geom_bar(aes(y = Abundance, x = Host, fill = Type), data = vyt.typ.l, stat="identity", position = "fill") +  scale_y_continuous(labels = percent_format(), expand = c(0, 0))  + labs( y= "Percentage of sequences", x="") + scale_fill_manual(values=COLtype1)  + theme_bw()  + theme(axis.text.x=element_text(angle=90,hjust=1), plot.title = element_text(hjust = 0.5), strip.background =element_rect(fill="white"), text = element_text(size=12),  legend.position="bottom") + coord_flip() + facet_grid(~sample)
type2=ggplot() +geom_bar(aes(y = Abundance, x = Host, fill = Type2), data = vyt.typ.l, stat="identity", position = "fill") +  scale_y_continuous(labels = percent_format(), expand = c(0, 0))  + labs( y= "Percentage of sequences", x="") + scale_fill_manual(values=COLtype2)  + theme_bw()  + theme(axis.text.x=element_text(angle=90,hjust=1), plot.title = element_text(hjust = 0.5), strip.background =element_rect(fill="white"), text = element_text(size=12),  legend.position="bottom") + coord_flip() + facet_grid(~sample)
pdf("./outputs/RedSeaVirome_viral_types.pdf",  width = 7, height = 7, pointsize = 12) 
grid.all=grid.arrange(type1,type2, ncol=1, nrow =2)
dev.off()


##find out what ssRNA viruses are in metagenomes
# library(reshape2)
# RNAvirs=subset(vyt.typ.l, vyt.typ.l$Type2 == "ssRNA" | vyt.typ.l$Type2 == "dsRNA")
# RNAvirs.l=aggregate(RNAvirs$Abundance, by = list(RNAvirs$sample, RNAvirs$Family), FUN =  sum)



