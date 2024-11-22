---
title: " Study II - Taxonomic visualisation of metaphlan profiles"
author: "Florentin Constancias"
date: "November 20, 2024"
output: 
  html_document: 
    toc: yes
    keep_md: yes
---











``` r
load(here::here("../../data/processed_data/metaphlan/01_data.Rdata"))
```


``` r
pd <- position_dodge(0.3)
```

Source the function `phyloseq_top_heatmap_barplot`:


``` r
source("https://raw.githubusercontent.com/fconstancias/KU_Caroline/refs/heads/main/code/functions/phyloseq_functions.R")
```

below the Arguments of the function phyloseq_top_heatmap_barplot() + default values:


``` r
# parameters:
  
#' @param ps_up A phyloseq object containing microbiome data.
#' @param group_var Character string for the grouping variable (default: "Sample_Time").
#' @param tax_levels Character vector of taxonomic levels to include in plots (default: c("Phylum", "Family", "Genus", "Species")).
#' @param ntax Integer specifying the number of taxa to display at each level (default: 5).
#' @param ntax_species Integer specifying the number of species-level taxa to display (default: 14).
#' @param plot_heights Numeric vector for the heights of the plot panels (default: c(1.4, 1.5, 4)).
#' @param plot_x Character string for the x-axis label in bar plots (default: "Subject").
#' @param facet_by Character vector specifying variables to use for faceting in heatmaps (default: c("Sample_Type", "Time")).
#' @param group_by Character vector for grouping in heatmaps (default: c("Sample_Type", "Time")).
#' @param facet_heat Formula for faceting heatmap plots (default: "~ Sample_Type + Time").
#' @param facet_formula Formula for faceting bar plots (default: "Sample_Type ~ Time").
#' @param rm_unclassified Logical indicating if unclassified taxa should be removed (default: TRUE).
#' @param barplot_level Character string for taxonomic level in bar plots (default: "Species").
#' @param boxplot_main_group Character string for the main grouping variable in box plots (default: "Class").
```

Run the function:


``` r
ps_up %>% 
  phyloseq_top_heatmap_barplot(facet_formula = "Sample_Type ~ Time" , 
                               ntax = 5, ntax_species = 10, plot_heights = c(1.4, 1.4, 4),
                               boxplot_main_group = "Family") -> out
```

Output is a list which containts different objects:


``` r
ls(out)
```

```
## [1] "bar_plot"      "heat"          "heat_all"      "heat_legend"  
## [5] "most_ab_treat" "nested_legend" "p"             "ps_sub"
```



``` r
out$heat_all
```

<img src="04_metaphlan_taxa_plot_files/figure-html/unnamed-chunk-7-1.png" style="display: block; margin: auto;" />


``` r
out$bar_plot + 
  ggpubr::rotate_x_text(60)
```

<img src="04_metaphlan_taxa_plot_files/figure-html/unnamed-chunk-8-1.png" style="display: block; margin: auto;" />


``` r
out$p
```

<img src="04_metaphlan_taxa_plot_files/figure-html/unnamed-chunk-9-1.png" style="display: block; margin: auto;" />


``` r
out$nested_legend
```

<img src="04_metaphlan_taxa_plot_files/figure-html/unnamed-chunk-10-1.png" style="display: block; margin: auto;" />


``` r
out$ps_sub
```

<img src="04_metaphlan_taxa_plot_files/figure-html/unnamed-chunk-11-1.png" style="display: block; margin: auto;" />


``` r
out$ps_sub + facet_null() +
  facet_grid(rows = vars(Sample_Type), 
             cols = vars(cluster_Dtp2, Subject )) + theme(legend.position = "none") 
```

<img src="04_metaphlan_taxa_plot_files/figure-html/unnamed-chunk-12-1.png" style="display: block; margin: auto;" />


``` r
ps_up %>% 
  subset_samples(Sample_Type == "Saliva") %>% 
  phyloseq_top_heatmap_barplot(facet_formula = "Sample_Type + cluster_Dtp2 ~ Time" , 
                               ntax = 5, ntax_species = 10, plot_heights = c(1.4, 1.4, 4),
                               boxplot_main_group = "Family") -> out2

out2$p
```

<img src="04_metaphlan_taxa_plot_files/figure-html/unnamed-chunk-13-1.png" style="display: block; margin: auto;" />


``` r
save(out, 
     ps_up,
     sample_pal,
     time_pal, sub_pal, sex_pal, sample_pal, period_pal,
     file = here::here("../../data/processed_data/metaphlan/04_data.Rdata"))
```



``` r
sessionInfo()
```

```
## R version 4.4.0 (2024-04-24)
## Platform: aarch64-apple-darwin20
## Running under: macOS Sonoma 14.6.1
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
##  [1] ggnested_0.1.1       rstatix_0.7.2        ampvis2_2.8.9       
##  [4] ggpubr_0.6.0         reshape2_1.4.4       scales_1.3.0        
##  [7] compositions_2.0-8   nlme_3.1-166         microViz_0.12.4     
## [10] speedyseq_0.5.3.9021 phyloseq_1.50.0      readxl_1.4.3        
## [13] lubridate_1.9.3      forcats_1.0.0        stringr_1.5.1       
## [16] dplyr_1.1.4          purrr_1.0.2          readr_2.1.5         
## [19] tidyr_1.3.1          tibble_3.2.1         ggplot2_3.5.1       
## [22] tidyverse_2.0.0     
## 
## loaded via a namespace (and not attached):
##   [1] RColorBrewer_1.1-3      tensorA_0.36.2.1        rstudioapi_0.17.1      
##   [4] jsonlite_1.8.9          magrittr_2.0.3          farver_2.1.2           
##   [7] rmarkdown_2.29          zlibbioc_1.51.2         vctrs_0.6.5            
##  [10] multtest_2.61.0         htmltools_0.5.8.1       broom_1.0.7            
##  [13] cellranger_1.1.0        Rhdf5lib_1.27.0         Formula_1.2-5          
##  [16] rhdf5_2.49.0            sass_0.4.9              bslib_0.8.0            
##  [19] htmlwidgets_1.6.4       plyr_1.8.9              plotly_4.10.4          
##  [22] cachem_1.1.0            igraph_2.1.1            lifecycle_1.0.4        
##  [25] iterators_1.0.14        pkgconfig_2.0.3         Matrix_1.7-1           
##  [28] R6_2.5.1                fastmap_1.2.0           GenomeInfoDbData_1.2.13
##  [31] digest_0.6.37           colorspace_2.1-1        S4Vectors_0.43.2       
##  [34] rprojroot_2.0.4         seriation_1.5.6         vegan_2.6-8            
##  [37] labeling_0.4.3          fansi_1.0.6             timechange_0.3.0       
##  [40] httr_1.4.7              abind_1.4-8             mgcv_1.9-1             
##  [43] compiler_4.4.0          here_1.0.1              withr_3.0.2            
##  [46] backports_1.5.0         carData_3.0-5           highr_0.11             
##  [49] ggsignif_0.6.4          MASS_7.3-61             bayesm_3.1-6           
##  [52] biomformat_1.34.0       permute_0.9-7           tools_4.4.0            
##  [55] ape_5.8                 glue_1.8.0              rhdf5filters_1.17.0    
##  [58] fantaxtic_0.2.1         gridtext_0.1.5          grid_4.4.0             
##  [61] Rtsne_0.17              cluster_2.1.6           ade4_1.7-22            
##  [64] generics_0.1.3          gtable_0.3.6            microbiome_1.28.0      
##  [67] tzdb_0.4.0              ca_0.71.1               data.table_1.16.2      
##  [70] hms_1.1.3               xml2_1.3.6              car_3.1-3              
##  [73] utf8_1.2.4              XVector_0.45.0          BiocGenerics_0.52.0    
##  [76] ggrepel_0.9.6           foreach_1.5.2           pillar_1.9.0           
##  [79] robustbase_0.99-4-1     splines_4.4.0           ggtext_0.1.2           
##  [82] lattice_0.22-6          survival_3.7-0          tidyselect_1.2.1       
##  [85] registry_0.5-1          Biostrings_2.73.2       knitr_1.48             
##  [88] IRanges_2.39.2          stats4_4.4.0            xfun_0.49              
##  [91] Biobase_2.66.0          DEoptimR_1.1-3          stringi_1.8.4          
##  [94] UCSC.utils_1.2.0        lazyeval_0.2.2          yaml_2.3.10            
##  [97] evaluate_1.0.1          codetools_0.2-20        cli_3.6.3              
## [100] munsell_0.5.1           jquerylib_0.1.4         Rcpp_1.0.13-1          
## [103] GenomeInfoDb_1.42.0     parallel_4.4.0          viridisLite_0.4.2      
## [106] plotwidgets_0.5.1       crayon_1.5.3            rlang_1.1.4            
## [109] cowplot_1.1.3           TSP_1.2-4
```
