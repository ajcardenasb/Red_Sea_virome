# Red Sea coral viromes

This repository contains the scripts used to analyze data and create figures for the manuscript “Coral-associated viral assemblages from the central Red Sea align with host species and contribute to the holobiont’s bacterial virulence”.  

##Summary
Viruses are the most abundant biological entities in marine ecosystems and their prevailing association with marine invertebrates, including corals, is assumed to have consequences for host health. However, comparatively little is known about function and diversity of coral-associated viruses. In this study we provide a first insight into the viral assemblages of Red Sea corals that also represents the first comprehensive survey of DNA and RNA viruses across 14 cnidarian species. We confirm the prevalence of several viral families associated with corals locally and globally. Functional annotation of viral reads suggests the contribution of phage-mediated transduction of bacterial virulence genes to the coral holobiont as a mechanism to increase genetic diversity. 

##Workflow

1. QC and Trimming were done using the script `` 
2. Overall taxonomic profiles were obtained using [link](https://github.com/vrmarcelino/CCMetagen)
3. Host (Cnidarian) reads and rRNA were removed usign the script ``
4. KAIJU was run using the script ``
5. Host relative abundances were analyzed using the scrip `RedSeaVirome_host_abundances.R`
6. Viral abundace tables were obtained usign the script `redSeaVirome_QC.R`
7. Viral diversity analysis and plots were done using the scripts `RedSeaVirome_barplots.R`, `redSeaVirome_biplots.R` and `RedSeaVirome_Ordination_permanovas.R`
8. Reads annotated to viruses were retrieve from KAIJU output using the script `RedSeaVirome_kaiju2viral_reads.R`
9. MG-RAST annotations were processed using the script `redSeaVirome_mg-rast.R`
