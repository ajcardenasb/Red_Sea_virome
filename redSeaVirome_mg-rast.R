library(reshape2)
library(data.table)
library(dplyr)
library(scales)

#############
## new
############

#setwd("~/Documents/GitHub/Red_Sea_virome/")

seed=read.csv("~/Input_files/Subsystems_hierarchy.csv")
my_ann=read.table("~/Input_files/my_annotations", fill = TRUE, sep = "\t", header = TRUE)
temp1=colsplit(my_ann$semicolon.separated.list.of.annotations,',',c('SSaccession.number', 'function'))
temp1$SSaccession.number=gsub("accession=|\\[|\\]","", temp1$SSaccession.number)
temp2=colsplit(my_ann$sequence.id,'_',c('Sample', 'Host', 'Header'))
temp2$Sample=gsub("mgm4881329.3\\|","", temp2$Sample)
full_ann=cbind(temp2, temp1$SSaccession.number)
colnames(full_ann)[4]="SSaccession.number"
anno_seed=merge(full_ann,seed, by.x="SSaccession.number", by.y="accession")
lv1.long=anno_seed %>% count(level1, Sample, Host)

## remove low Q samples
goodSamples=read.table("./Input_files/MG_MT_filtered_counts", header = TRUE, row.names = 1)
lv1.long.filt=subset(lv1.long, paste(lv1.long$Sample, lv1.long$Host, sep = "_") %in% colnames(goodSamples))
lv1.long.filt$Sample=gsub("MG", "Metagenomes", lv1.long.filt$Sample)
lv1.long.filt$Sample=gsub("MT", "Metatranscriptomes", lv1.long.filt$Sample)

#col20=c('#d73027','#f46d43','#fdae61','#fee08b','#ffffbf','#d9ef8b','#a6d96a','#66bd63','#1a9850', '#c51b7d','#de77ae','#f1b6da','#fde0ef', "#d1e5f0", "#92c5de", "#4393c3", "#2166ac" , "#053061" , "#023858", "#023858")
#ggplot() +geom_bar(aes(x=Host, y=n, fill=level1), data = lv1.long, stat="identity", position = "fill") +  scale_y_continuous(labels = percent_format(), expand = c(0, 0)) + facet_grid(~Sample)  + labs( y= "", x="") + scale_fill_manual(values=col)  + theme_bw()  + theme(axis.text.x=element_text(angle=90,hjust=1), plot.title = element_text(hjust = 0.5), strip.background =element_rect(fill="white"), text = element_text(size=12)) + coord_flip()

#top 10 categories

lv1.frec=lv1.long.filt %>% count(level1)
lv1.top10=lv1.frec[order(lv1.frec$n,decreasing = TRUE),][1:10,]
lv1.top10.long=subset(lv1.long.filt, lv1.long.filt$level1 %in%lv1.top10$level1)

col10=c('#8e0152','#c51b7d','#de77ae','#f1b6da','#fde0ef','#e6f5d0','#b8e186','#7fbc41','#4d9221','#276419')


pdf("./outputs/RedSeaVirome_level1Subsystems.pdf",  width = 10, height = 5, pointsize = 12) 
ggplot() +geom_bar(aes(x=Host, y=n, fill=level1), data = lv1.top10.long, stat="identity", position = "fill") +  scale_y_continuous(labels = percent_format(), expand = c(0, 0)) + facet_grid(~Sample)  + labs( y= "", x="") + scale_fill_manual(values=col10)  + theme_bw()  + theme(axis.text.x=element_text(angle=90,hjust=1), plot.title = element_text(hjust = 0.5), strip.background =element_rect(fill="white"), text = element_text(size=12), legend.position="bottom") + coord_flip()
dev.off()


plot=ggplot() +geom_bar(aes(x=Host, y=n, fill=level1), data = lv1.top10.long, stat="identity", position = "fill") +  scale_y_continuous(labels = percent_format(), expand = c(0, 0)) + facet_grid(~Sample)  + labs( y= "", x="") + scale_fill_manual(values=col10)  + theme_bw()  + theme(axis.text.x=element_text(angle=90,hjust=1), plot.title = element_text(hjust = 0.5), strip.background =element_rect(fill="white"), text = element_text(size=12), legend.position="bottom") + coord_flip()
leg.bar=get_legend(plot)

pdf("./outputs/RedSeaVirome_legend_Figure3.pdf",  width = 20, height = 3, pointsize = 8) 
grid.arrange(leg.bar)
dev.off()

#stats
# QCstats=anno_seed %>% count(Sample, Host)
# 
# viral=anno_seed %>% count(Sample, level1, level3)
# 
# phages=viral %>% filter(level1 == "Phages, Prophages, Transposable elements, Plasmids")
# lvl3.p=dcast(phages,level3~Sample  , value.var="n")
# 
# viru=viral %>% filter(level1 == "Virulence, Disease and Defense")
# lvl3.v=dcast(viru,level3~Sample  , value.var="n")
# 
# viral1=anno_seed %>% count(Sample, level1)
# lvl1=dcast(viral1,level1~Sample  , value.var="n")
# 
# viral2=anno_seed %>% count(Sample, level2)
# lvl1=dcast(viral1,level1~Sample  , value.var="n")
# 
# test=viral %>% filter(level1 == "Phages, Prophages, Transposable elements, Plasmids")
# 
