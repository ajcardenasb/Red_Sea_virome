
setwd("~/Documents/GitHub/Red_Sea_virome/kaiju2viralreads/")
kai=read.table("/Input_files/class_ACA1.family.out2", fill = TRUE,  nrows=10000)
taxid=read.table("viral_taxID.txt",  nrows=10000)

viral=subset(kai, kai$V3 %in% taxid$V1)$V2


### pulling viral reads out
## Rscript kaiju2viral_reads.R kaiju_output fasta_out. THe viral_taxID file is a single column with viral taxID numbers
args=commandArgs(TRUE)
message("Reading Kaiju output file : ", args[1])
kai=read.table(args[1])
message("Reading virus TaxIDs")
taxid=read.table("viral_taxID")
message("Done!")
viral=subset(kai, kai$V3 %in% taxid$V1)$V2
write.table(viral, args[2], row.names = FALSE, col.names = FALSE)


### command line
args=commandArgs(TRUE)
message("Reading Kaiju output file : ", args[1])
kai=read.table(args[1])
message("Reading virus TaxIDs")
taxid=read.table("viral_taxID")
message("Done!")
viral=subset(kai, kai$V3 %in% taxid$V1)$V2
write.table(viral, args[2], row.names = FALSE, col.names = FALSE)