setwd("~/Documents/GitHub/Red_Sea_virome/")


dna=read.table("./Input_files/MG_log10_counts", header = TRUE, row.names = 1)
rna=read.table("./Input_files/MT_log10_counts", header = TRUE, row.names = 1)
met=read.table("./Input_files/redSea_virome_metadata", header = TRUE, row.names = 1, fill = TRUE)


### Biplots
dna.pca=rda(t(dna),  distance = "bray")
biplot(dna.pca, display = c('species', "sites"), type = c("text","points"), main="MG")
mg.scores=as.data.frame(dna.pca[1])
mg.scores$Familes=rownames(mg.scores)
topFamilies.mg=mg.scores[order(mg.scores$colsum, decreasing = TRUE),][1:10,]
dna.filtered=subset(dna, rownames(dna) %in% rownames(topFamilies.mg))
dna.met=merge(t(dna.filtered), met, by="row.names")
dna.pca2=rda(dna.met[,2:11],  distance = "bray")

rna.pca=rda(t(rna),  distance = "bray")
biplot(rna.pca, display = c('species', "sites"), type = c("text","points"), main="mt")
mt.scores=as.data.frame(rna.pca[1])
mt.scores$Familes=rownames(mt.scores)
topFamilies.mt=mt.scores[order(mt.scores$colsum, decreasing = TRUE),][1:10,]
rna.filtered=subset(rna, rownames(rna) %in% rownames(topFamilies.mt))
rna.met=merge(t(rna.filtered), met, by="row.names")
rna.pca2=rda(rna.met[,2:11],  distance = "bray")

dna.sp=c("#8DD3C7", "#518C4B", "#9fb9bf","#FB8072","#80B1D3","#FDB462","#B3DE69", "#CCEBC5","#D9D9D9","#BC80BD","#FCCDE5","#BEBADA","#b88c8c")
rna.sp=c("#8DD3C7", "#518C4B", "#9fb9bf","#FB8072","#80B1D3","#FDB462","#B3DE69", "#CCEBC5","#D9D9D9","#BC80BD","#FCCDE5","#BEBADA","#b88c8c")
life=c("#A54657", "#C8963E", "#058ED9", "#4CB275", "#476A6F")

pdf("./outputs/RedSeaVirome_biplots.pdf",  width = 10, height = 10, pointsize = 10) 
par(mfrow=c(2,2)) 
biplot(dna.pca2, display = c('species', "sites"), type = c("text","points"), main="Host species", scaling = "species")
ordihull(dna.pca2, group = dna.met$Species,  draw = "polygon" , col = dna.sp)#,  label = TRUE to see groups
spp.names <- levels(dna.met$Species)
legend("bottomleft",   legend = spp.names, fill = dna.sp, cex=0.5 )

biplot(dna.pca2, display = c('species', "sites"), type = c("text","points"), main="Life history strategies", scaling = "species")
ordihull(dna.pca2, group = dna.met$Life_history_strategy,  draw = "polygon" , col = life)
life.names <- levels(dna.met$Life_history_strategy)
legend("bottomleft",   legend = life.names, fill = life, cex=0.5 )

# biplot(dna.pca2, display = c('species', "sites"), type = c("text","points"), main="Growth form")
# ordihull(dna.pca2, group = dna.met$Growth_form_typical,  draw = "polygon" , col = c("#8DD3C7", "#518C4B", "#9fb9bf","#FB8072","#80B1D3","#FDB462","#B3DE69", "#CCEBC5","#D9D9D9","#BC80BD","#FCCDE5","#BEBADA","#b88c8c"))
# growth.names <- levels(dna.met$Growth_form_typical)
# legend("bottomleft",   legend = growth.names, fill = c("#8DD3C7", "#518C4B", "#9fb9bf","#FB8072","#80B1D3","#FDB462","#B3DE69", "#CCEBC5","#D9D9D9","#BC80BD","#FCCDE5","#BEBADA","#b88c8c"), cex=0.5 )


biplot(rna.pca2, display = c('species', "sites"), type = c("text","points"), main="Host species", scaling = "species")
ordihull(rna.pca2, group = rna.met$Species,  draw = "polygon" , col = rna.sp) # ,  label = TRUE to see groups
legend("bottomleft",   legend = spp.names, fill = rna.sp, cex=0.5 )

biplot(rna.pca2, display = c('species', "sites"), type = c("text","points"), main="Life history strategies", scaling = "species")
ordihull(rna.pca2, group = rna.met$Life_history_strategy,  draw = "polygon" , col =life)
legend("bottomleft",   legend = life.names, fill = life, cex=0.5 )
# 
# biplot(rna.pca2, display = c('species', "sites"), type = c("text","points"), main="Growth form")
# ordihull(rna.pca2, group = rna.met$Growth_form_typical,  draw = "polygon" , col = c("#8DD3C7", "#518C4B", "#9fb9bf","#FB8072","#80B1D3","#FDB462","#B3DE69", "#CCEBC5","#D9D9D9","#BC80BD","#FCCDE5","#BEBADA","#b88c8c"))
# legend("bottomleft",   legend = growth.names, fill = c("#8DD3C7", "#518C4B", "#9fb9bf","#FB8072","#80B1D3","#FDB462","#B3DE69", "#CCEBC5","#D9D9D9","#BC80BD","#FCCDE5","#BEBADA","#b88c8c"), cex=0.5 )

dev.off()



