#library(compositions)
library(vegan)
library(pairwiseAdonis)
library(reshape2)
library(phyloseq)
library(microbiome)
library(tidyr)
library(gridExtra)
library(ggpubr)


setwd("~/Documents/GitHub/Red_Sea_virome/")

#data in
dna=read.table("./Input_files/MG_filtered_counts", header = TRUE, row.names = 1)
rna=read.table("./Input_files/MT_filtered_counts", header = TRUE, row.names = 1)
met=read.table("./Input_files/redSea_virome_metadata", header = TRUE, row.names = 1, fill = TRUE)
met$Life_history_strategy=factor(met$Life_history_strategy, levels = c("Competitive", "Generalist", "Stress-tolerant", "Weedy", "outgroup"))

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

####Ordination plots
COLdna=c("#FDB462","#BC80BD","#BEBADA","#80B1D3","#AFAFAF","#FCCDE5","#B3DE69","#CCEBC5","#FFED6F","#8DD3C7","#FB8072","#9FB9BF","#b88c8c") 
COLrna=c("#518C4B","#BC80BD","#80B1D3","#FCCDE5","#B3DE69","#8DD3C7","#FB8072","#9FB9BF") ##FFED6F


## RDA
dna.rda=ordinate(dna.t , "RDA", "bray" , ~Species )  # 18.5 + 11.2.3%
a=plot_ordination(dna.t,dna.rda , color = "Species", shape = "Life_history_strategy")  + geom_point(size = 3) + theme_bw()  + ggtitle("Metagenomes") + theme(plot.title = element_text(hjust = 0.5)) + theme(legend.title = element_blank(), legend.position="none")   + scale_shape_manual(values=c(15,19,17,18,8,19 )) + scale_colour_manual(values=COLdna) 
rna.rda=ordinate(rna.t, "RDA" , "bray", ~Species ) # 18.5 + 11.2.3%
b=plot_ordination(rna.t,rna.rda , color = "Species", shape = "Life_history_strategy")  + geom_point(size = 3) + theme_bw()  + ggtitle("Metatranscriptomes") + theme(plot.title = element_text(hjust = 0.5)) + theme(legend.title = element_blank(), legend.position="none")   + scale_shape_manual(values=c(15,19,17,18,8 )) + scale_colour_manual(values=COLrna) 

names.bar=plot_ordination(dna.t,dna.rda , color = "Species", shape = "Life_history_strategy")  + geom_point(size = 3) + theme_bw()  + ggtitle("Metagenomes") + theme(plot.title = element_text(hjust = 0.5)) + theme(legend.title = element_blank(), legend.position="right")   + scale_shape_manual(values=c(15,19,17,18,8,19 )) + scale_colour_manual(values=COLdna) 
leg.bar=get_legend(names.bar)

#plot

pdf("./outputs/RedSeaVirome_ordination_rda.pdf",  width = 7, height =5, pointsize = 12) 
grid.arrange(a,b, ncol=2, nrow=1)
dev.off()

pdf("./outputs/RedSeaVirome_ordination_rda_legend.pdf",  width = 20, height = 5, pointsize = 12) 
grid.arrange(leg.bar, ncol=1, nrow=1)
dev.off()


permutest(dna.rda)
#######################
### inertia heatmaps ##
#######################

trait.rda=read.table("./Input_files/RDA_inertia.txt", header = TRUE, row.names = 1)
colnames(trait.rda)=gsub("_", " ", colnames(trait.rda))

colheat=c( "#f7fcf0", "#e5f3e1", "#d0ead5", "#b9e1cc", "#a0d8c7", "#81c9c0", "#61babc", "#3dabb9", "#0091b1", "#0077a6", "#005c96", "#084081") # 12 colors

library(gplots)
pdf("./outputs/RedSeaVirome_rda_inertia.pdf",  width = 10, height = 5, pointsize = 12) 
heatmap.2(t(as.matrix(trait.rda)), scale = "col",  col =  colheat, dendrogram = "none",   notecol="black", trace="none", density.info= "none",  key.title = "", key.xlab = "Z-score for adjusted adonis R2", srtCol = 0, cexRow=1,cexCol=1,margins=c(5,15),  Rowv=FALSE,Colv=FALSE)
dev.off()

#######################
### Beta dispersion ###
#######################

rownames(mg.log.met)=mg.log.met[,1]
bray.dna=vegdist(mg.log.met[2:118], method = "bray")
beta_dna=lapply(mg.log.met[, 120:ncol(mg.log.met)], function(x) betadisper(bray.dna, x))
lapply(beta_dna, function(x) permutest(x, permutations = 99, pairwise = TRUE)) # eveyrthing is significant except for Major Clade, life history strategy and sexual system


########################
### Full model adonis###
########################

mg.log10=as.data.frame(t(otu_table(dna.t)))
mg.log.met=merge(mg.log10, met, by = "row.names")


adonis(mg.log.met[2:118] ~ mg.log.met$Abbreviation+mg.log.met$Life_history_strategy)

outgrup=ifelse(mg.log.met$Major_Clade ==  "outgroup", "outgrup", "coral" )
test2=mg.log.met$Life_history_strategy

adonis(mg.log10 ~ mg.log.met$Life_history_strategy + mg.log.met$Sexual_system, method = "bray", perm = 999)

full_ado_dna=lapply(mg.log.met[, 120:ncol(mg.log.met)], function(x) adonis(mg.log.met[2:118] ~ x))
lapply(names(full_ado_dna), function(x) write.table(full_ado_dna[[x]][1], paste("./outputs/Adonis_Results/MG_", x, "_fullAdo.txt", sep = ""),col.names=FALSE, row.names=FALSE, sep="\t", quote=FALSE))

mt.log10=as.data.frame(t(otu_table(rna.t)))
mt.log.met=merge(mt.log10, met, by = "row.names")

full_ado_rna=lapply(mt.log.met[, c(120:121,123:128,130)], function(x) adonis(mt.log.met[2:118] ~ x)) #removed substrate, to have only factors with 2 or more levels
lapply(names(full_ado_rna), function(x) write.table(full_ado_rna[[x]][1], paste("./outputs/Adonis_Results/MT_", x, "_fullAdo.txt", sep = ""),col.names=FALSE, row.names=FALSE, sep="\t", quote=FALSE))

trait.cor=read.table("./Input_files/adonis_R2.txt", header = TRUE, row.names = 1)
colnames(trait.cor)=gsub("_", " ", colnames(trait.cor))

#pval=read.table("./Input_files/adonis_pval.txt", header = TRUE, row.names = 1)
#colnames(trait.cor)=gsub("_", " ", colnames(trait.cor))
#colheat=c('#f7fcf0','#e0f3db','#ccebc5','#a8ddb5','#7bccc4','#4eb3d3','#2b8cbe','#0868ac','#084081') # 9 colors
colheat=c( "#f7fcf0", "#e5f3e1", "#d0ead5", "#b9e1cc", "#a0d8c7", "#81c9c0", "#61babc", "#3dabb9", "#0091b1", "#0077a6", "#005c96", "#084081") # 12 colors

library(gplots)
pdf("./outputs/RedSeaVirome_adonisR2.pdf",  width = 10, height = 5, pointsize = 12) 
heatmap.2(t(as.matrix(trait.cor)), scale = "col",  col =  colheat, dendrogram = "none",   notecol="black", trace="none", density.info= "none",  key.title = "", key.xlab = "Z-score for adjusted adonis R2", srtCol = 0, cexRow=1,cexCol=1,margins=c(5,15),  Rowv=FALSE,Colv=FALSE)
dev.off()


############################
### Pairwise permanovas ###
###########################


traits_pair.ado_dna=lapply(mg.log.met[, 119:ncol(mg.log.met)], function(x) pairwise.adonis(mg.log.met[, 2:118], x, p.adjust.m ='fdr', sim.method = 'bray'))
all.pair.ado.dna=rbind(traits_pair.ado_dna$Species, traits_pair.ado_dna$Family,  traits_pair.ado_dna$Order,  traits_pair.ado_dna$Major_Clade,  traits_pair.ado_dna$Minor_Clade,  traits_pair.ado_dna$Life_history_strategy,  traits_pair.ado_dna$Growth_form_typical,  traits_pair.ado_dna$Mode_of_larval_development,  traits_pair.ado_dna$Sexual_system,  traits_pair.ado_dna$Substrate)
write.table(all.pair.ado.dna, "./outputs/Adonis_Results/MG_all_adoPairWise.txt", row.names=FALSE, sep="\t", quote=FALSE)

traits_pair.ado_rna=lapply(mt.log.met[, c(120:121,123:128,130)], function(x) pairwise.adonis(mt.log.met[, 2:118], x, p.adjust.m ='fdr', sim.method = 'bray'))
all.pair.ado.dna=rbind(traits_pair.ado_rna$Species, traits_pair.ado_rna$Family,  traits_pair.ado_rna$Order,  traits_pair.ado_rna$Major_Clade,  traits_pair.ado_rna$Minor_Clade,  traits_pair.ado_rna$Life_history_strategy,  traits_pair.ado_rna$Growth_form_typical,  traits_pair.ado_rna$Mode_of_larval_development,  traits_pair.ado_rna$Sexual_system,  traits_pair.ado_rna$Substrate)
write.table(all.pair.ado.dna, "./outputs/Adonis_Results/MT_all_adoPairWise.txt", row.names=FALSE, sep="\t", quote=FALSE)


##############################
### Constrained ordination ###
##############################

dna.rda=ordinate(dna.t, "RDA", "bray")  # 18.5 + 11.2.3%
plot_ordination(dna.t,dna.rda , color = "Species", shape = "Life_history_strategy", type="split", label="Abbreviation")  + geom_point(size = 3) + theme_bw()  + ggtitle("Metagenomes") + theme(plot.title = element_text(hjust = 0.5)) + theme(legend.title = element_blank(), legend.position="bottom")   + scale_shape_manual(values=c(15,19,17,18,8,19 )) + scale_colour_manual(values=COLdna) 

rna.rda=ordinate(rna.t, "RDA" , "bray") # 18.5 + 11.2.3%
b=plot_ordination(rna.t,rna.rda , color = "Species", shape = "Life_history_strategy")  + geom_point(size = 3) + theme_bw()  + ggtitle("Metagenomes") + theme(plot.title = element_text(hjust = 0.5)) + theme(legend.title = element_blank(), legend.position="bottom")   + scale_shape_manual(values=c(15,19,17,18,8 )) + scale_colour_manual(values=COLrna) 
grid.arrange(a,b,  ncol=2, nrow=1)

dna.log=otu_table(dna.t)
dna.rda2=rda(t(dna.log) ,  distance = "bray")
dna.log.met=merge(t(dna.log), met , by = "row.names")
rda(dna.log.met[,2:118] ~ dna.log.met$Species + dna.log.met$Life_history_strategy,  distance = "bray")

rna.log=otu_table(rna.t)
rna.rda2=rda(t(rna.log),  distance = "bray")
rna.log.met=merge(t(rna.log), met , by = "row.names")


dna.sp=c("#8DD3C7", "#518C4B", "#9fb9bf","#FB8072","#80B1D3","#FDB462","#B3DE69", "#CCEBC5","#D9D9D9","#BC80BD","#FCCDE5","#BEBADA","#b88c8c")
rna.sp=c("#8DD3C7", "#518C4B", "#9fb9bf","#FB8072","#80B1D3","#FDB462","#B3DE69", "#CCEBC5","#D9D9D9","#BC80BD","#FCCDE5","#BEBADA","#b88c8c")
life=c("#A54657", "#C8963E", "#058ED9", "#4CB275", "#476A6F")


##Life strategies
par(mfrow=c(2,1)) 
biplot(dna.rda, display =  "sites", type = c("text","points"), main="Metagenomes", scaling = "species", xlim=c(-0.7, 0.7), ylim=c(-1, 1), expand=10)
plot(dna.rda,  main="Metagenomes")
ordihull(dna.rda, group = dna.log.met$Life_history_strategy,  draw = "polygon" , col = life)

biplot(rna.rda, display =  "sites", type = c("text","points"), main="Metatranscriptomes", scaling = "species")
ordihull(rna.rda, group = rna.log.met$Life_history_strategy,  draw = "polygon" , col = life)


##species
par(mfrow=c(2,1)) 
biplot(dna.rda, display =  "sites", type = c("text","points"), main="Metagenomes", scaling = "species")
ordihull(dna.rda, group = dna.log.met$Species,  draw = "polygon" , col = dna.sp)

biplot(rna.rda2, display =  "sites", type = c("text","points"), main="Metatranscriptomes", scaling = "species")
ordihull(rna.rda2, group = rna.log.met$Species,  draw = "polygon" , col = rna.sp)
biplot(rna.rda, display =  "sites", type = c("text","points"), main="Metatranscriptomes", scaling = "species")
ordihull(rna.rda, group = rna.log.met$Species,  draw = "polygon" , col = rna.sp)


