
setwd("~/Documents/Bioinformatics_scripts/R_scripts/Red_Sea_viromes/kaiju2viralreads/")

library(stringr)
library(dplyr)
library(data.table)
library(reshape2)
library(ggplot2)
library(scales)
library(gridExtra)
library(ggpubr)

################
## CCmetagene ##
################

ccmet=read.csv("ccmetagene_class", header = TRUE)
#met<-read.table("./Input_files/redSea_virome_metadata", header = TRUE)

ccmet.spk=aggregate(ccmet[,1:101], by = list(ccmet[, 105]), FUN =  sum) 
ccmet.l=melt(ccmet.spk, keep.rownames = TRUE)
ccmet.l$variable=gsub("[0-9]_results.res", "", ccmet.l$variable)
ccmet.l$variable=gsub("[0-9]_ccmetagenOut.res", "", ccmet.l$variable)
ccmet.l$sample=ifelse(ccmet.l$variable %like% "MG", "Metagenomes", "Metatranscriptomes")
ccmet.l$variable=gsub("M[TG]_", "", ccmet.l$variable)
topClass=aggregate(value ~ Group.1, ccmet.l , sum)[order(aggregate(value ~ Group.1, ccmet.l , sum)$value,decreasing = TRUE),][c(1:2,4:6,8:21),1]
ccmet.l$Group.1=ifelse(ccmet.l$Group.1 %in% topClass,  as.character(ccmet.l$Group.1), gsub(".*", "zOthers", ccmet.l$Group.1))

P21=c("#771155",  "#CC99BB", "#114477", "#4477AA", "#77AADD", "#117777", "#44AAAA", "#77CCCC", "#117744", "#44AA77", "#88CCAA", "#777711", "#AAAA44", "#DDDD77", "#774411", "#AA7744", "#DDAA77", "#771122", "#AA4455", "#C0C0C0")
a=ggplot() +geom_bar(aes(y = value, x = variable, fill = Group.1), data = ccmet.l, stat="identity", position = "fill", width=1) +  scale_y_continuous(labels = percent_format(), expand = c(0, 0))  + labs( y= "", x="") + scale_fill_manual(values=P21)  + theme_bw()  + theme(axis.text.x=element_text(angle=90,hjust=1), plot.title = element_text(hjust = 0.5), strip.background =element_rect(fill="white"), text = element_text(size=12),  legend.position="none") + facet_grid(.~sample)
names.a=ggplot() +geom_bar(aes(y = value, x = variable, fill = Group.1), data = ccmet.l, stat="identity", position = "fill", width=1) +  scale_y_continuous(labels = percent_format(), expand = c(0, 0))  + labs( y= "", x="") + scale_fill_manual(values=P21)  + theme_bw()  + theme(axis.text.x=element_text(angle=90,hjust=1), plot.title = element_text(hjust = 0.5), strip.background =element_rect(fill="white"), text = element_text(size=12),  legend.position="right") + facet_grid(.~sample)
leg.a=get_legend(names.a)

test=ccmet.spk[c(8,26),50:ncol(ccmet.spk)]

############
## KAIJU ##
###########

kai=read.table("all_kaiju_fullReport3", sep = "\t", header = T)
kai.s=subset(kai, file %like% "MG" | file %like% "MT")
kai.s[c('superkingdom','phylum','class','order','family','genus','species')] <- colsplit(kai.s$taxon_name,";",c('superkingdom','phylum','class','order','family','genus','species'))
kai.s$reads=as.numeric(as.character(kai.s$reads))
kai.s$superkingdom=gsub("cannot be assigned to.*", "unclassified", kai.s$superkingdom)
kai.s$superkingdom=gsub("Eukaryota", "Microbial Eukaryotes", kai.s$superkingdom)
supk=aggregate(reads ~ file+superkingdom, kai.s , sum)
supk$file=gsub("[0-9].family.out2", "", supk$file)
supk$sample=ifelse(supk$file %like% "MG", "Metagenomes", "Metatranscriptomes")
supk$file=gsub("M[TG]_", "", supk$file)


# P5=c("#AA4488","#F9E176","#117777","#C2C2C2","#ff7f0e")
# a=ggplot() +geom_bar(aes(y = reads, x = file, fill = superkingdom), data = supk, stat="identity", position = "fill", width=1) +  scale_y_continuous(labels = percent_format(), expand = c(0, 0))  + labs( y= "", x="") + scale_fill_manual(values=P5)  + theme_bw()  + theme(axis.text.x=element_text(angle=90,hjust=1), plot.title = element_text(hjust = 0.5), strip.background =element_rect(fill="white"), text = element_text(size=12),  legend.position="none")  + facet_grid(.~sample)
# names.a=ggplot() +geom_bar(aes(y = reads, x = file, fill = superkingdom), data = supk, stat="identity", position = "fill", width=1) +  scale_y_continuous(labels = percent_format(), expand = c(0, 0))  + labs( y= "", x="") + scale_fill_manual(values=P5)  + theme_bw()  + theme(axis.text.x=element_text(angle=90,hjust=1), plot.title = element_text(hjust = 0.5), strip.background =element_rect(fill="white"), text = element_text(size=12),  legend.position="right")  + facet_grid(.~sample)
# leg.a=get_legend(names.a)


#euks
kai.euk=subset(kai.s, superkingdom == "Microbial Eukaryotes")
euk=aggregate(reads ~ file+class, kai.euk , sum)
euk$file=gsub("[0-9].family.out2", "", euk$file)
euk$sample=ifelse(euk$file %like% "MG", "Metagenomes", "Metatranscriptomes")
euk$file=gsub("M[TG]_", "", euk$file)
topClass=aggregate(reads ~ class, euk , sum)[order(aggregate(reads ~ class, euk , sum)$reads,decreasing = TRUE),][1:20,1]

euk$class=ifelse(euk$class %in% topClass,  as.character(euk$class), gsub(".*", "zOthers", euk$class))

P21=c("#771155", "#AA4488", "#CC99BB", "#114477", "#4477AA", "#77AADD", "#117777", "#44AAAA", "#77CCCC", "#117744", "#44AA77", "#88CCAA", "#777711", "#AAAA44", "#DDDD77", "#774411", "#AA7744", "#DDAA77", "#771122", "#AA4455", "#C0C0C0")
b=ggplot() +geom_bar(aes(y = reads, x = file, fill = class), data = euk, stat="identity", position = "fill", width=1) +  scale_y_continuous(labels = percent_format(), expand = c(0, 0))  + labs( y= "", x="") + scale_fill_manual(values=P21)  + theme_bw()  + theme(axis.text.x=element_text(angle=90,hjust=1), plot.title = element_text(hjust = 0.5), strip.background =element_rect(fill="white"), text = element_text(size=12),  legend.position="none") + facet_grid(.~sample)
names.b=ggplot() +geom_bar(aes(y = reads, x = file, fill = class), data = euk, stat="identity", position = "fill", width=1) +  scale_y_continuous(labels = percent_format(), expand = c(0, 0))  + labs( y= "", x="") + scale_fill_manual(values=P21)  + theme_bw()  + theme(axis.text.x=element_text(angle=90,hjust=1), plot.title = element_text(hjust = 0.5), strip.background =element_rect(fill="white"), text = element_text(size=12),  legend.position="right") + facet_grid(.~sample)
leg.b=get_legend(names.b)

#bac
kai.bac=subset(kai.s, superkingdom == "Bacteria")
bac=aggregate(reads ~ file+family, kai.bac , sum)
bac$file=gsub("[0-9].family.out2", "", bac$file)
bac$sample=ifelse(bac$file %like% "MG", "Metagenomes", "Metatranscriptomes")
bac$file=gsub("M[TG]_", "", bac$file)
topClass=aggregate(reads ~ family, bac , sum)[order(aggregate(reads ~ family, bac , sum)$reads,decreasing = TRUE),][1:20,1]
bac$family=ifelse(bac$family %in% topClass,  as.character(bac$family), gsub(".*", "zOthers", bac$family))

P21=c("#771155", "#AA4488", "#CC99BB", "#114477", "#4477AA", "#77AADD", "#117777", "#44AAAA", "#77CCCC", "#117744", "#44AA77", "#88CCAA", "#777711", "#AAAA44", "#DDDD77", "#774411", "#AA7744", "#DDAA77", "#771122", "#AA4455", "#C0C0C0")
c=ggplot() +geom_bar(aes(y = reads, x = file, fill = family), data = bac, stat="identity", position = "fill", width=1) +  scale_y_continuous(labels = percent_format(), expand = c(0, 0))  + labs( y= "", x="") + scale_fill_manual(values=P21)  + theme_bw()  + theme(axis.text.x=element_text(angle=90,hjust=1), plot.title = element_text(hjust = 0.5), strip.background =element_rect(fill="white"), text = element_text(size=12),  legend.position="none")  + facet_grid(.~sample) 
names.c=ggplot() +geom_bar(aes(y = reads, x = file, fill = family), data = bac, stat="identity", position = "fill", width=1) +  scale_y_continuous(labels = percent_format(), expand = c(0, 0))  + labs( y= "", x="") + scale_fill_manual(values=P21)  + theme_bw()  + theme(axis.text.x=element_text(angle=90,hjust=1), plot.title = element_text(hjust = 0.5), strip.background =element_rect(fill="white"), text = element_text(size=12),  legend.position="right")  + facet_grid(.~sample)
leg.c=get_legend(names.c)

pdf(file = "../outputs/otherHosts_relAbunda.pdf",  width = 5, height = 7, pointsize = 8) 
grid.arrange(a,b,c,ncol=1,nrow=3)
dev.off()

pdf(file = "../outputs/otherHosts_relAbunda_legends.pdf",  width = 5, height = 10, pointsize = 8) 
grid.arrange(leg.a,leg.b,leg.c,ncol=1,nrow=3)
dev.off()



#stats

stats=aggregate(reads ~ file+superkingdom, kai.s , sum)
stats$file=gsub(".family.out2", "", stats$file)
stats.w=dcast(stats, file ~ superkingdom, value.var= "reads", fun.aggregate = sum)
stats.w$Sum=rowSums(stats.w[,2:6])
write.table(stats.w, "stats_taxonomic_assigments_kaiju", row.names = F, quote = F)

bac.fam=kai.bac[,c(1:4,10)]
bac.fam$file=gsub(".family.out2", "", bac.fam$file)
write.table(stats.w, "../outputs/bacterial_abundances", row.names = F, quote = F)

