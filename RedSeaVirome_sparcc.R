setwd("~/Documents/Bioinformatics_scripts/R_scripts/Red_Sea_viromes/")
library(Hmisc)


kai=read.table("./kaiju2viralreads/all_kaiju_fullReport3", sep = "\t", header = T,)
vir=as.data.frame(t(read.table("./Input_files/MG_MT_filtered_counts", header = TRUE,)))
kai.s=subset(kai, file %like% "MG" | file %like% "MT")
kai.s[c('superkingdom','phylum','class','order','family','genus','species')] <- colsplit(kai.s$taxon_name,";",c('superkingdom','phylum','class','order','family','genus','species'))
kai.s$reads=as.numeric(as.character(kai.s$reads))
kai.s$superkingdom=gsub("cannot be assigned to.*", "unclassified", kai.s$superkingdom)
kai.s$superkingdom=gsub("Eukaryota", "Microbial Eukaryotes", kai.s$superkingdom)
fam=aggregate(reads ~ file+family, kai.s , sum)
fam$family=gsub("^$|^ $", 'Unclassified', fam$family)
fam$file=gsub(".family.out2", "", fam$file)
badSamples=c("MG_ACA5","MG_ACR1","MG_ACR2","MG_ACR3","MG_ACR4","MG_ACR5","MG_GAL2","MG_PAC1","MG_POR1","MG_XEN2","MT_ACR1","MT_PAC5","MT_POR1","MT_POR5", "MT_STY5")
fam.s=subset(fam, !fam$file %in% badSamples)
fam.m=dcast(fam.s, file ~family, value.var= "reads", sum )
cor=fam.m[,-1]
rownames(cor)=fam.m[,1]

merged=merge(cor,vir, by= "row.names",all.x = T)
all=merged[,-1]
rownames(all)=merged[,1]

all.n=sweep(all,1,rowSums(all),`/`)
all.n=log10(all+1)

SpearCC2 = rcorr(as.matrix(all.n), type = "spearman")
long_SpearCC = melt(SpearCC2$r, value.name = c("SCC"))
long_SpearCC$pval=melt(SpearCC2$P, value.name = c("pval"))$pval
long_SpearCC$Superkingdom1=kai.s$superkingdom[match(long_SpearCC$Var1,kai.s$family)]

### Phages
phages=c("Myoviridae", "Picobirnaviridae", "Podoviridae", "Siphoviridae")
bac_SCC=subset(long_SpearCC, Superkingdom1 == "Bacteria" & Var2 %in% phages & SCC > 0.8 )
bac_SCC$Class=paste(kai.s$class, kai.s$phylum, sep = " - ")[match(bac_SCC$Var1,kai.s$family)]

inverte=c("Baculoviridae", "Caulimoviridae", "Flaviviridae","Iridoviridae", "Marseilleviridae","Mimiviridae", "Parvoviridae", "Phycodnaviridae", "Pithoviridae", "Polydnaviridae", "Poxviridae")
inv_SCC=subset(long_SpearCC, Superkingdom1 == "Microbial Eukaryotes" & Var2 %in% inverte & SCC > 0.8)
inv_SCC$Class=paste(kai.s$class, kai.s$phylum, sep = " - ")[match(inv_SCC$Var1,kai.s$family)]


write.table(inv_SCC, "./outputs/Spearman_corr_coh_Euka",quote = F, row.names = F)

test=inv_SCC %>% group_by(Var2, Class) %>% tally()
test2=test %>% group_by(Class) %>% tally()

