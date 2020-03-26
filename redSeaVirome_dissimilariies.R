library(harrietr)
library(ggtree)
setwd("~/Documents/GitHub/Red_Sea_virome")


##intraspecific variation
met=read.table("./Input_files/redSea_virome_metadata", header = TRUE, row.names = 1, fill = TRUE)


dna=read.table("MG_log10_counts", header = TRUE)
dna.dist= as.matrix(vegdist(t(dna), method  = "bray", diag=TRUE))
dna.pairs=melt_dist(dna.dist, order = NULL, dist_name = "MG")
dna.pairs$pair1=gsub("[0-9]", "", dna.pairs$iso1)
dna.pairs$pair2=gsub("[0-9]", "", dna.pairs$iso2)
dna.pairs.filt=unique(subset(dna.pairs, pair1 == pair2 ))
dna.pairs.filt$clade1=met$Life_history_strategy[match(dna.pairs.filt$iso1, rownames(met))]
dna.pairs.filt$clade2=met$Life_history_strategy[match(dna.pairs.filt$iso2, rownames(met))]
a=ggplot(dna.pairs.filt,aes(x=pair1, y=MG)) +  geom_boxplot() + geom_jitter(width=0.25, alpha=0.5) + labs( y= "Bray-Curtis dissimilarities", title="Metagenomes" ,x="") + theme(legend.title = element_blank(), axis.text.x = element_text(angle = 90, hjust = 1), panel.grid.major = element_line(colour = "transparent"), panel.grid.minor = element_line(colour = "transparent"), panel.background = element_blank(), axis.line = element_line(colour = "black"),  legend.position="none") + ggtitle("") +  ylim(0, 0.8)#+ scale_fill_manual(values=COLrna) + coord_flip() 
c=ggplot(dna.pairs.filt,aes(x=clade1, y=MG)) +  geom_boxplot() + geom_jitter(width=0.25, alpha=0.5) + labs( y= "Bray-Curtis dissimilarities", title="" , x="") + theme(legend.title = element_blank(), axis.text.x = element_text(angle = 90, hjust = 1), panel.grid.major = element_line(colour = "transparent"), panel.grid.minor = element_line(colour = "transparent"), panel.background = element_blank(), axis.line = element_line(colour = "black"),  legend.position="none") + ggtitle("") +  ylim(0, 0.8)#+ scale_fill_manual(values=COLrna) + coord_flip() 

rna=read.table("./Input_files/MT_log10_counts", header = TRUE)
rna.dist= as.matrix(vegdist(t(rna), method  = "bray", diag=TRUE))
rna.pairs=melt_dist(rna.dist, order = NULL, dist_name = "MT")
rna.pairs$pair1=gsub("[0-9]", "", rna.pairs$iso1)
rna.pairs$pair2=gsub("[0-9]", "", rna.pairs$iso2)
rna.pairs.filt=subset(rna.pairs, pair1 == pair2 )
rna.pairs.filt$clade1=met$Life_history_strategy[match(rna.pairs.filt$iso1, rownames(met))]
rna.pairs.filt$clade2=met$Life_history_strategy[match(rna.pairs.filt$iso2, rownames(met))]

b=ggplot(rna.pairs.filt,aes(x=pair1, y=MT)) +  geom_boxplot() + geom_jitter(width=0.25, alpha=0.5) + labs( y= "", x="", title="Metatranscriptomes" ) + theme(legend.title = element_blank(), axis.text.x = element_text(angle = 90, hjust = 1), panel.grid.major = element_line(colour = "transparent"), panel.grid.minor = element_line(colour = "transparent"), panel.background = element_blank(), axis.line = element_line(colour = "black"),  legend.position="none") + ggtitle("")+  ylim(0, 0.8) #+ scale_fill_manual(values=COLrna) + coord_flip() 
d=ggplot(dna.pairs.filt,aes(x=clade1, y=MG)) +  geom_boxplot() + geom_jitter(width=0.25, alpha=0.5) + labs( y= "", x="" ) + theme(legend.title = element_blank(), axis.text.x = element_text(angle = 90, hjust = 1), panel.grid.major = element_line(colour = "transparent"), panel.grid.minor = element_line(colour = "transparent"), panel.background = element_blank(), axis.line = element_line(colour = "black"),  legend.position="none") + ggtitle("")+  ylim(0,0.8) #+ scale_fill_manual(values=COLrna) + coord_flip() 


pdf("./outputs/RedSeaVirome_brayDissimilarities.pdf",  width = 7, height = 7, pointsize = 10) 
grid.all=grid.arrange(a,b,c,d, ncol=2, nrow =2, top = grid::textGrob("A",x=1,hjust=0))
dev.off()

