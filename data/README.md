Description of data sets:  

1. `condensed_species_NCBI.txt`: This file contains the latest database release from [Madin et al.](). This file is organized as a taxa-by-trait matrix with NCBI identifiers. Current timestamp:  
2. `gevers_biom.biom`: This is the biom file for the Gevers et al. microbiome data set and can be used as inputs to `picrust2`.    
3. `gevers_dada2.rds`: This is the `seqtabnochim` data export after running the Gevers et al. data set through `dada2`.  
4. `gevers_seq.fa`: This file contains the unique amplicon sequence variants for the Gevers et al. data set that will be used as inputs to `picrust2`.  
5. `madin_proc.rds`: This file contains pre-processed database from the Madin et al. comprehensive database.  
6. `pathabundances_3.tsv.gz`, `taxonomic_profiles_3.tsv.gz`: This contains all the latest pathway abundances and taxonomic profiles for the integrative human microbiome project. These are results obtained from the latest biobakery version (v3).  
7. `weissman.csv`: This is 