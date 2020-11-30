#sprot
sprot_metadata=read.table("sprot_accession_complementary_information", sep = "\t",  quote = "", header = T)
sprot_annotations=read.table("sprot_annotations", sep = "\t",  quote = "", header = F)
sprot_annotations_f=sprot_annotations[,c(1,2,3,12,13)]
colnames(sprot_annotations_f) =c("Sample", "Query_sequence_identifier", "UniProt_accession","E-val","Bit_score")
sprot_annotations_f$Sample=gsub("_viral_contigs_mmseqsDB_viral_sprot_outfmt6","",sprot_annotations_f$Sample)

sprot_annotations_f$Protein_name=sprot_metadata$Protein.names[match(sprot_annotations_f$UniProt_accession,sprot_metadata$Entry)]
sprot_annotations_f$Taxonomic_lineage=sprot_metadata$Taxonomic.lineage..ALL.[match(sprot_annotations_f$UniProt_accession,sprot_metadata$Entry)]
sprot_annotations_f$GO_term=sprot_metadata$Gene.ontology.IDs[match(sprot_annotations_f$UniProt_accession,sprot_metadata$Entry)]
sprot_annotations_f$GO_BP=sprot_metadata$Gene.ontology..biological.process.[match(sprot_annotations_f$UniProt_accession,sprot_metadata$Entry)]
sprot_annotations_f$GO_MF=sprot_metadata$Gene.ontology..molecular.function.[match(sprot_annotations_f$UniProt_accession,sprot_metadata$Entry)]

#trembl
trembl_metadata=read.table("trembl_accession_complementary_information", sep = "\t",  quote = "", header = T)
trembl_annotations=read.table("trembl_annotations", sep = "\t",  quote = "", header = F)
trembl_annotations_f=trembl_annotations[,c(1,2,3,12,13)]
colnames(trembl_annotations_f) =c("Sample", "Query_sequence_identifier", "UniProt_accession","E-val","Bit_score")
trembl_annotations_f$Sample=gsub("_viral_contigs_mmseqsDB_viral_trembl_outfmt6","",trembl_annotations_f$Sample)

trembl_annotations_f$Protein_name=trembl_metadata$Protein.names[match(trembl_annotations_f$UniProt_accession,trembl_metadata$Entry)]
trembl_annotations_f$Taxonomic_lineage=trembl_metadata$Taxonomic.lineage..ALL.[match(trembl_annotations_f$UniProt_accession,trembl_metadata$Entry)]
trembl_annotations_f$GO_term=trembl_metadata$Gene.ontology.IDs[match(trembl_annotations_f$UniProt_accession,trembl_metadata$Entry)]
trembl_annotations_f$GO_BP=trembl_metadata$Gene.ontology..biological.process.[match(trembl_annotations_f$UniProt_accession,trembl_metadata$Entry)]
trembl_annotations_f$GO_MF=trembl_metadata$Gene.ontology..molecular.function.[match(trembl_annotations_f$UniProt_accession,trembl_metadata$Entry)]

#compilation
supple_trembl=subset(trembl_annotations_f, !trembl_annotations_f$Query_sequence_identifier %in% sprot_annotations_f$Query_sequence_identifier)
final=rbind(sprot_annotations_f, supple_trembl)
write.table(final, "final_viral_annotations_Uniprot", quote = F, row.names = F, sep="\t")

## Go terms
split_go=final %>% mutate(GO_term = strsplit(as.character(GO_term), ";")) %>% unnest(GO_term)
go_summary=split_go %>% group_by(GO_term) %>% tally()
go_summary$GO_term=sub("^$", "", go_summary$GO_term)
write.table(go_summary, "annotated_GO_terms", quote = F, row.names = F, sep="\t", col.names = F)

split_go=final %>% mutate(GO_term = strsplit(as.character(GO_term), ";")) %>% unnest(GO_term)
split_go$GO_term=gsub("\t", "", split_go$GO_term)
split_go$GO_term=gsub(" ", "", split_go$GO_term)
#split_go$GO_term=gsub(":", "", split_go$GO_term)
#split_go$GO_term=sub(":", "", split_go$GO_term)
write.table(split_go$GO_term, "annotated_GO_terms", quote = F, row.names = F,  col.names = F)


#### Barplots of most abundant GO BPs
final2=final
final2$GO_BP=gsub(".*;", "", final2$GO_BP)
final2$GO_BP=gsub("^ ", "", final2$GO_BP)
top_BP=final2 %>% group_by( GO_BP) %>% tally()
top_BP_s=top_BP[order(top_BP$n,decreasing = TRUE),][2:21,1]

barp=final2 %>% group_by(Sample, GO_BP) %>% tally()
barp$GO_BP=sub("^$", "Unclassified", barp$GO_BP)
barp_2=barp[!barp$GO_BP == "Unclassified",]
barp_2$GO_BP=ifelse(barp_2$GO_BP %in% top_BP_s$GO_BP, as.character(barp_2$GO_BP), gsub(".*","Less abundant categories", barp_2$GO_BP))
barp_2$Sample=gsub("[0-9]", "", barp_2$Sample)
barp_2$Library=ifelse(barp_2$Sample %like% "MG_", "Metagenomes", "Metatranscriptomes")
barp_2$Sample=gsub("M[TG]_", "", barp_2$Sample)

P21=c("#771155", "#AA4488", "#CC99BB", "#114477", "#4477AA", "#77AADD", "#117777", "#44AAAA", "#77CCCC", "#117744", "#44AA77", "#88CCAA", "#777711", "#AAAA44", "#DDDD77", "#774411", "#AA7744", "#DDAA77", "#771122", "#AA4455", "#C0C0C0")
pdf("Barplots_UniProt-GO.pdf",  width = 10, height =5, pointsize = 10) 
ggplot() +geom_bar(aes(y = n, x = Sample, fill = GO_BP), data = barp_2, stat="identity", position = "fill")  + labs( y= "Gene Onltology (GO) biological processes", x="", fill = "GO biological process") + scale_fill_manual(values=P21) + facet_grid(Library~.) +  theme_classic() + theme( legend.key.size = unit(0.2, "cm"),legend.key.width = unit(0.2,"cm"), legend.position = 'bottom', legend.key = element_blank(), strip.background = element_blank()) + guides(fill=guide_legend(ncol=3)) + scale_y_continuous( expand = c(0, 0)) 
dev.off()

## table
tab=final
tab$Sample=ifelse(tab$Sample %like% "MG_", "Metagenomes", "Metatranscriptomes")
tab$GO_BP=gsub(";.*", "", tab$GO_BP)
tab_stat=tab %>% group_by(Sample, GO_BP) %>% tally()
tab_fin=reshape2::dcast(tab_stat, GO_BP~Sample,value.var= "n", fun.aggregate = sum)
write.table(tab_fin, "TableS7", quote = F, row.names = F,  col.names = T, sep = "\t")


