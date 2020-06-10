#!/bin/bash

###########
### QC  ###
###########

fastqc *HT1_kBR2_[12]P -o QC_trimmed/ -t 25

########################
### Trimming reads   ###
########################

#metatranscriptomes
for i in *R1.fastq;do java -jar /home/yej/software/Trimmomatic-0.36/trimmomatic-0.36.jar PE $i ${i//R1/R2} -baseout  ~/Ana_2/1_Trimmed_MG/${i//_R1.fastq/_HT1_kBR2} ILLUMINACLIP:/home/yej/software/Trimmomatic-0.36/adapters/Trueseq_HT.fa:2:30:10:8:TRUE LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:40;done
#metagenomes
for i in *R1.fastq;do java -jar /home/yej/software/Trimmomatic-0.36/trimmomatic-0.36.jar PE $i ${i//R1/R2} -baseout  ~/Ana_2/1_Trimmed_MG/${i//_R1.fastq/_HT1_kBR2} ILLUMINACLIP:/home/yej/software/Trimmomatic-0.36/adapters/Trueseq_HT.fa:2:30:10:8:TRUE LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:40;done


#############################
### Removing host reads   ###
#############################

#metatranscriptomes
for i in *HT1_kBR2_1P; do /home/cardena/software/bbmap/bbsplit.sh in1=$i in2=${i//1P/2P} basename=${i//HT1_kBR2_1P/cnidarian}_%_#.fq outu=../2_bbsplit/${i//HT1_kBR2_1P/clean}_%_#.fq ref=~/software/Databases/cnidarian.fa threads=26 ; done
#metageomes
for i in *HT1_kBR2_1P; do /home/cardena/software/bbmap/bbsplit.sh in1=$i in2=${i//1P/2P} basename=${i//HT1_kBR2_1P/_cnidarian}_%_#.fq outu=../2_bbsplit/${i//HT1_kBR2_1P/_clean}_%_#.fq  ref=/home/cardena/database/contaminants_anny/cnidarian.fa threads=16 ; done

#############################
### Removing rRNA reads  ###
#############################

for i in *_nonrRNA.fastq.fastq;do /home/yej/software/sortmerna-2.1b/scripts/unmerge-paired-reads.sh $i ${i//.fastq.fastq/_1P} ${i//.fastq.fastq/_2P};done
