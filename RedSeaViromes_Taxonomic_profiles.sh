###########################################################
### CCmetagen https://github.com/vrmarcelino/CCMetagen  ###
###########################################################
export PATH="/home/cardena/miniconda3/bin:$PATH"
source activate ccmetagen
export PATH="/home/cardena/software/CCMetagen:$PATH"

#metagenomes on raw reads @/home/cardena/projects/redSeaViruses/ccmetagene/raw_mg
for i in /home/cardena/projects/redSeaViruses/raw_mg/quality_checked/*_R1_paired.fastq  ; do kma -ipe $i ${i//_R1_paired.fastq/_R2_paired.fastq} -o ${i//_R1_paired.fastq/_ccmetagenOut} -t_db /home/cardena/databases/nt_kma_indexed/compress_ncbi_nt/ncbi_nt -t 24 -1t1 -mem_mode -and -apm f ; done
##running
for i in *_ccmetagenOut.res ; do CCMetagen.py -i $i -o ${i//_ccmetagenOut/_results} ; done
CCMetagen_merge.py -i . -l Class -o class_table
CCMetagen_merge.py -i . -l Phylum -o phylum_table

#metatrasncriptomes on raw reads
for i in /home/cardena/projects/redSeaViruses/raw_mt/*_kBR2_1P  ; do kma -ipe $i ${i//_kBR2_1P/_kBR2_2P} -o ${i//_HT1_kBR2_1P/_ccmetagenOut} -t_db /home/cardena/databases/nt_kma_indexed/compress_ncbi_nt/ncbi_nt -t 200 -1t1 -mem_mode -and -apm f ; done
for i in *_ccmetagenOut.res ; do CCMetagen.py -i $i -o ${i//__ccmetagenOut/_results} ; done
CCMetagen_merge.py -i . -l Class -o class_table
CCMetagen_merge.py -i . -l Phylum -o phylum_table


########################
### Kaiju extended   ###
########################

#metagenomes
for i in /home/yej/Ana_2/2_bbsplit/*_clean_%_1.fq; do /home/yej/software/kaiju/bin/kaiju -z 16 -t nodes.dmp -f ./nr_euk/kaiju_db_nr_euk.fmi -i $i -j ${i//1.fq/2.fq} -o ${i//_clean_%_1.fq/.family.out2} -v; done
#metatrasncriptomes
for i in /home/yej/Ana_2_MT/2_5_sortmeRNA/NonrRNA/*merged_nonrRNA_1P; do /home/yej/software/kaiju/bin/kaiju -z 20 -t nodes.dmp -f ./nr_euk/kaiju_db_nr_euk.fmi -i $i -j ${i//1P/2P} -o ${i//merged_nonrRNA_1P/.family.out2} -v; done

#using  kaiju outputs to create tables using an old nodes.dmp
for i in *.family.out2; do /home/cardena/software/kaiju/bin/kaiju2table -t /home/cardena/databases/kaiju_old_Jin/nodes.dmp -n  /home/cardena/databases/kaiju_old_Jin/names.dmp  -o ${i//.family.out2/_report_full2} -r class -l superkingdom,phylum,class,order,family,genus,species  $i; done
for i in *.family.out2; do /home/cardena/software/kaiju/bin/kaiju2table -t /home/cardena/databases/kaiju_old_Jin/nodes.dmp -n  /home/cardena/databases/kaiju_old_Jin/names.dmp  -o ${i//.family.out2/_report_full3} -r family -l superkingdom,phylum,class,order,family,genus,species  $i; done

cat *_report_full3 > all_kaiju_fullReport #(phylum)
cat *_report_full3 > all_kaiju_fullReport2 #(class)
cat *_report_full3 > all_kaiju_fullReport3 #(family)

### Trying taxonkit
export PATH="/home/cardena/miniconda3/bin:$PATH"
source activate taxonkit
cat *.classified | cut -f3 | sort | uniq > all_taxIDs_found
#42701 were found in total, of which 4445 are Euks
taxonkit lineage --data-dir /home/cardena/databases/NCBI_tax/ all_taxIDs_found -j 56 -o lineage_annnotations
for i in *.classified ; do awk '{print FILENAME"\t"$3}' $i > ${i}_2 ; done
cat *classified_2 > kaiju_out_all_taxID
