---
title: " Study II - Importing functional data from gene identified in the co-assembly to phyloseq"
author: "Florentin Constancias"
date: "January 10, 2025"
output: 
  html_document: 
    toc: yes
    keep_md: yes
---








After running SQM pipeline: 


``` bash
(SqueezeMeta) [ljc444@esrumhead01fl scratch]$ ~/.conda/envs/SqueezeMeta/SqueezeMeta/utils/sqm2tables.py 010_coassembly/Trial2-co-assembly Trial2-co-assembly_SQM_tables
```

You would need to downlaod the files locally from /maps/projects/hansen_ol-AUDIT//data/Trial2-co-assembly_SQM/Trial2-co-assembly/results/tables

These are the files:

├── Trial2-co-assembly.COG.abund.tsv
├── Trial2-co-assembly.COG.bases.tsv
├── Trial2-co-assembly.COG.copyNumber.tsv
├── Trial2-co-assembly.COG.cov.tsv
├── Trial2-co-assembly.COG.names.tsv
├── Trial2-co-assembly.COG.tpm.tsv
├── Trial2-co-assembly.KO.abund.tsv
├── Trial2-co-assembly.KO.bases.tsv
├── Trial2-co-assembly.KO.copyNumber.tsv
├── Trial2-co-assembly.KO.cov.tsv
├── Trial2-co-assembly.KO.names.tsv
├── Trial2-co-assembly.KO.tpm.tsv
├── Trial2-co-assembly.PFAM.abund.tsv
├── Trial2-co-assembly.PFAM.bases.tsv
├── Trial2-co-assembly.PFAM.copyNumber.tsv
├── Trial2-co-assembly.PFAM.cov.tsv
├── Trial2-co-assembly.PFAM.tpm.tsv
├── Trial2-co-assembly.RecA.tsv
├── Trial2-co-assembly.bin.tax.tsv
├── Trial2-co-assembly.class.allfilter.abund.tsv
├── Trial2-co-assembly.class.nofilter.abund.tsv
├── Trial2-co-assembly.class.prokfilter.abund.tsv
├── Trial2-co-assembly.contig.sequences.tsv
├── Trial2-co-assembly.contig.tax.allfilter.tsv
├── Trial2-co-assembly.contig.tax.nofilter.tsv
├── Trial2-co-assembly.contig.tax.prokfilter.tsv
├── Trial2-co-assembly.family.allfilter.abund.tsv
├── Trial2-co-assembly.family.nofilter.abund.tsv
├── Trial2-co-assembly.family.prokfilter.abund.tsv
├── Trial2-co-assembly.genus.allfilter.abund.tsv
├── Trial2-co-assembly.genus.nofilter.abund.tsv
├── Trial2-co-assembly.genus.prokfilter.abund.tsv
├── Trial2-co-assembly.order.allfilter.abund.tsv
├── Trial2-co-assembly.order.nofilter.abund.tsv
├── Trial2-co-assembly.order.prokfilter.abund.tsv
├── Trial2-co-assembly.orf.sequences.tsv
├── Trial2-co-assembly.orf.tax.allfilter.tsv
├── Trial2-co-assembly.orf.tax.nofilter.tsv
├── Trial2-co-assembly.orf.tax.prokfilter.tsv
├── Trial2-co-assembly.phylum.allfilter.abund.tsv
├── Trial2-co-assembly.phylum.nofilter.abund.tsv
├── Trial2-co-assembly.phylum.prokfilter.abund.tsv
├── Trial2-co-assembly.species.allfilter.abund.tsv
├── Trial2-co-assembly.species.nofilter.abund.tsv
├── Trial2-co-assembly.species.prokfilter.abund.tsv
├── Trial2-co-assembly.superkingdom.allfilter.abund.tsv
├── Trial2-co-assembly.superkingdom.nofilter.abund.tsv
└── Trial2-co-assembly.superkingdom.prokfilter.abund.tsv


# Import and get ready:

## metaphlan data - so we have metadata and colors:


``` r
load(here::here("../../data/processed_data/metaphlan/01_data.Rdata") )
```

## SQM tables:


``` r
#' Convert SQM Table to Phyloseq Object
#'
#' This function reads files from a specified directory and constructs a Phyloseq object
#' using metric and taxonomy information derived from SQM tables.
#'
#' @param dir Character string specifying the directory containing the SQM table files. Default: NULL.
#' @param db Character string specifying the database prefix used in file naming (e.g., "COG"). Default: "COG".
#' @param metric Character string specifying the metric type (e.g., "tpm"). Default: "tpm".
#'
#' @return A `phyloseq` object constructed from the specified SQM tables.
#'
#' @details
#' The function searches the directory for two files:
#' - A metric file matching the pattern `db.metric.tsv`.
#' - A names file matching the pattern `db.names.tsv`.
#'
#' It reads the files, processes them into matrices, and constructs a `phyloseq` object.
#'
#' @examples
#' # Example usage (ensure the directory and files exist):
#' # ps <- SQM_table_2phyloseq(dir = "path/to/dir", db = "COG", metric = "tpm")
#'
```


``` r
SQM_table_2phyloseq(dir =  "/Volumes/Elements/Caroline_KU/full_coassembly/SQM/full_coassembly_SQM_tables/tables/", 
                                db = "COG", metric = "tpm") -> ps

ps
```

```
## phyloseq-class experiment-level object
## otu_table()   OTU Table:         [ 35240 taxa and 239 samples ]
## tax_table()   Taxonomy Table:    [ 35240 taxa by 2 taxonomic ranks ]
```



``` r
ps %>% 
  sample_names() %>% 
  str_extract(.,"[^_]+") %>% 
  str_remove_all(., "[A-Z]") %>% 
  paste0("S_",.) -> sample_names(ps)
```


``` r
ps %>% 
  physeq_add_metadata(metadata = ps_up %>% 
                        sample_data() %>% 
                        data.frame() %>% 
                        rownames_to_column('sample_name')) -> ps
```

Phyloseq object with functional gene info from the gene identified in the coassembly together with the metadata of the metaphlan phyloseq obejct. 



``` r
ps
```

```
## phyloseq-class experiment-level object
## otu_table()   OTU Table:         [ 35240 taxa and 239 samples ]
## sample_data() Sample Data:       [ 239 samples by 25 sample variables ]
## tax_table()   Taxonomy Table:    [ 35240 taxa by 2 taxonomic ranks ]
```


``` r
# save(period_pal,
#      sample_pal,
#      sex_pal,
#      sub_pal,
#      time_pal,
#      p_meta_delta, p_meta,p_meta_ratio, 
#      ps_up, clin_xy, clin_xy_d, pclust_3, pclust_2, file = here::here("../../data/processed_data/metaphlan/01_data.Rdata"))
```



``` r
# saveRDS(object = ps_up, file = here::here("../../data/processed_data/metaphlan/metaphlan4.1.1_phyloseq.RDS"))
```




``` r
sessionInfo()
```

```
## R version 4.4.0 (2024-04-24)
## Platform: aarch64-apple-darwin20
## Running under: macOS 15.2
## 
## Matrix products: default
## BLAS:   /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/lib/libRblas.0.dylib 
## LAPACK: /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/lib/libRlapack.dylib;  LAPACK version 3.12.0
## 
## locale:
## [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
## 
## time zone: Europe/Zurich
## tzcode source: internal
## 
## attached base packages:
## [1] stats     graphics  grDevices utils     datasets  methods   base     
## 
## other attached packages:
##  [1] phyloseq_1.50.0 readxl_1.4.3    lubridate_1.9.4 forcats_1.0.0  
##  [5] stringr_1.5.1   dplyr_1.1.4     purrr_1.0.2     readr_2.1.5    
##  [9] tidyr_1.3.1     tibble_3.2.1    ggplot2_3.5.1   tidyverse_2.0.0
## 
## loaded via a namespace (and not attached):
##  [1] ade4_1.7-22             tidyselect_1.2.1        Biostrings_2.74.1      
##  [4] fastmap_1.2.0           digest_0.6.37           timechange_0.3.0       
##  [7] lifecycle_1.0.4         cluster_2.1.8           survival_3.8-3         
## [10] magrittr_2.0.3          compiler_4.4.0          rlang_1.1.4            
## [13] sass_0.4.9              tools_4.4.0             igraph_2.1.2           
## [16] yaml_2.3.10             data.table_1.16.4       knitr_1.49             
## [19] bit_4.5.0.1             here_1.0.1              plyr_1.8.9             
## [22] withr_3.0.2             BiocGenerics_0.52.0     grid_4.4.0             
## [25] stats4_4.4.0            multtest_2.62.0         biomformat_1.34.0      
## [28] colorspace_2.1-1        Rhdf5lib_1.28.0         scales_1.3.0           
## [31] iterators_1.0.14        MASS_7.3-63             cli_3.6.3              
## [34] rmarkdown_2.29          vegan_2.6-8             crayon_1.5.3           
## [37] generics_0.1.3          rstudioapi_0.17.1       httr_1.4.7             
## [40] reshape2_1.4.4          tzdb_0.4.0              ape_5.8-1              
## [43] cachem_1.1.0            rhdf5_2.50.1            zlibbioc_1.52.0        
## [46] splines_4.4.0           parallel_4.4.0          cellranger_1.1.0       
## [49] XVector_0.46.0          vctrs_0.6.5             Matrix_1.7-1           
## [52] jsonlite_1.8.9          IRanges_2.40.1          hms_1.1.3              
## [55] S4Vectors_0.44.0        bit64_4.5.2             foreach_1.5.2          
## [58] jquerylib_0.1.4         glue_1.8.0              codetools_0.2-20       
## [61] stringi_1.8.4           gtable_0.3.6            GenomeInfoDb_1.42.1    
## [64] UCSC.utils_1.2.0        munsell_0.5.1           pillar_1.10.0          
## [67] htmltools_0.5.8.1       rhdf5filters_1.18.0     GenomeInfoDbData_1.2.13
## [70] R6_2.5.1                rprojroot_2.0.4         vroom_1.6.5            
## [73] evaluate_1.0.1          Biobase_2.66.0          lattice_0.22-6         
## [76] bslib_0.8.0             Rcpp_1.0.13-1           nlme_3.1-166           
## [79] permute_0.9-7           mgcv_1.9-1              xfun_0.49              
## [82] pkgconfig_2.0.3
```
