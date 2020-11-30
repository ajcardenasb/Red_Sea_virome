#### annotating viral reads
export PATH="/home/cardena/miniconda3/bin:$PATH"
source activate mmseqs2

for i in ~/databases/uniprot_viruses/*.fasta ; do mmseqs createdb $i $(basename $i | sed 's/.fasta/_mmseqsDB/')  ; done

for i in *_mmseqsDB ; do mmseqs search $i /share/databases/uniprot_sprot/mmseqs_uniprot_sprot/mmseqs_uniprot_sprot.targetDB ${i/_mmseqsDB/_sprot} tmp --max-accept 1  --threads 56 ; done
for i in *_mmseqsDB ; do mmseqs search $i /share/databases/uniprot_trembl/mmseqs_uniprot_trembl/mmseqs_uniprot_trembl.targetDB ${i/_mmseqsDB/_trembl} tmp --max-accept 1  --threads 56 ; done
#Running

for i in *_mmseqsDB ; do mmseqs convertalis  ${i} /share/databases/uniprot_sprot/mmseqs_uniprot_sprot/mmseqs_uniprot_sprot.targetDB  ${i/_mmseqsDB/_sprot}  ${i}_outfmt6 ; done
for i in *_mmseqsDB; do mmseqs convertalis  ${i} /share/databases/uniprot_trembl/mmseqs_uniprot_trembl/mmseqs_uniprot_trembl.targetDB  ${i/_mmseqsDB/_trembl} ${i}_trembl_outfmt6 ; done



## viral sprot
for i in *_mmseqsDB ; do mmseqs search $i /home/cardena/databases/uniprot_viruses/sequences/viral_uniprot_swissprot_mmseqsDB ./viral_databases/${i/_mmseqsDB/_viral_sprot} tmp --max-accept 1  --threads 64 ; done
for i in *_mmseqsDB ; do mmseqs convertalis  ${i} /home/cardena/databases/uniprot_viruses/sequences/viral_uniprot_swissprot_mmseqsDB  ./viral_databases/${i/_mmseqsDB/_viral_sprot}  ./viral_databases/${i}_viral_sprot_outfmt6 ; done

cat *_sprot_outfmt6 | cut -f2 | sort | uniq > sprot_accessions ## scp and then using the Uniprot mapping tool, got gene name, tax and GO terms
for i in *_viral_sprot_outfmt6 ; do awk '{print FILENAME"\t"$0}' $i > ${i//viral_contigs_mmseqsDB_viral_sprot_outfmt6/_final_output} ; done
cat *_final_output > sprot_annotations



## viral trembl
mmseqs createdb viral_uniprot_trembl.fasta viral_uniprot_trembl_mmseqsDB
for i in *_mmseqsDB ; do mmseqs search $i /home/cardena/databases/uniprot_viruses/sequences/viral_uniprot_trembl_mmseqsDB ./viral_databases/${i/_mmseqsDB/_viral_trembl} tmp --max-accept 1  --threads 64 ; done
for i in *_mmseqsDB ; do mmseqs convertalis  ${i} /home/cardena/databases/uniprot_viruses/sequences/viral_uniprot_trembl_mmseqsDB  ./viral_databases/${i/_mmseqsDB/_viral_trembl}  ./viral_databases/${i}_viral_trembl_outfmt6 ; done
cat *_trembl_outfmt6 | cut -f2 | sort | uniq > trembl_accessions
for i in *_viral_trembl_outfmt6 ; do awk '{print FILENAME"\t"$0}' $i > ${i//viral_contigs_mmseqsDB_viral_trembl_outfmt6/_final_trembl_output} ; done
cat *_final_trembl_output > trembl_annotations


#CirGO ## from ~/software/CirGO/Py3+/CirGO_Wind_Unix. #Input from REVIGO (treemap csv)
export PATH="/home/cardena/miniconda3/bin:$PATH"
source activate  cirGO
python CirGO.py -inputFile /home/cardena/projects/redSeaViruses/raw_mg/test_trimming_reviewers/viral_data_cirGO/REVIGO_treemap.csv -outputFile /home/cardena/projects/redSeaViruses/raw_mg/test_trimming_reviewers/viral_data_cirGO/redsea_virmome_GO.svg -fontSize 12 -numCat 40 -leg "GO biological processes"

python CirGO.py -inputFile /home/cardena/projects/redSeaViruses/raw_mg/test_trimming_reviewers/viral_data_cirGO/REVIGO_treemap.csv -outputFile /home/cardena/projects/redSeaViruses/raw_mg/test_trimming_reviewers/viral_data_cirGO/redsea_virmome_GO_50cat.svg -fontSize 10 -numCat 50 -leg "GO biological processes"
