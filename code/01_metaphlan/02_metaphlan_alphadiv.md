---
title: " Study II - Alpha-diversity analyses of metaphlan profiles"
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
source("https://raw.githubusercontent.com/fconstancias/KU_Caroline/refs/heads/main/code/functions/phyloseq_functions.R")

# compute_plot_alpha()
# details:

#' @title Compute Alpha Diversity Plot and Statistical Summaries
#' @description
#' This function generates plots and statistical summaries for alpha diversity measures 
#' over time across different sample groups. The function returns a list containing summary
#' statistics, plots, and results from statistical tests.
#' @author Florentin Constancias
#' @param alpha_long_df Data frame in long format containing alpha diversity values. 
#'        Should have columns specified by `alpha_measures`, `x`, `y`, `color_point`, etc.
#' @param alpha_measures Character vector of alpha diversity measures to analyze and plot.
#' @param pd Position dodge to adjust plot element spacing (default is 0.3).
#' @param x Character string for x-axis variable (e.g., "Time").
#' @param y Character string for y-axis variable (e.g., "value").
#' @param color_point Column name for coloring points (e.g., "Subject").
#' @param shape_point Column name for point shape (e.g., "Sample").
#' @param group_point Column name for grouping points (e.g., "interaction(Sample,Subject)").
#' @param group_boxplot Column name for grouping boxplot elements (e.g., "interaction(Sample,Time)").
#' @param group_line Column name for grouping line elements (e.g., "Subject").
#' @param facet_formula Formula for faceting plots (default is "alphadiversiy ~ Sample + cluster_Dtp2").
#' @param col_pal Color palette for points and lines.
#' @param fill_pal Fill palette for boxplots.
#' @param stat_formula Formula for statistical testing (default is "value ~ Time").
#' @param anova_test_formula Formula for ANOVA testing (default is "value ~ Time*Sample + Error(Subject/(Time*Sample))").
#' @param group_by_stats Variables for grouping in statistical comparisons.
#' @param padjust_method Method for p-value adjustment (default is "fdr").
#' @param stat_paired Logical; whether statistical tests are paired (default is FALSE).
#' @param ref_group_stat Reference group for pairwise statistical tests.
#' @return A list with summary data frames, ggplot objects, and statistical test results.
```

# Alpha div exploration:


``` r
# alpha_measures = c("observed", "diversity_shannon", "diversity_inverse_simpson", "diversity_coverage", "evenness_pielou")
alpha_measures = c("observed", "diversity_shannon", "diversity_inverse_simpson", "evenness_pielou")


ps_up %>% 
  subset_taxa(Class != "UNCLASSIFIED") %>% 
  transform_sample_counts(function(x) x/sum(x) * 1) %>%
  phyloseq_alphas() -> alphas


# alphas
```



``` r
alphas %>% 
  pivot_longer(cols = all_of(alpha_measures), values_to = "value", names_to = 'alphadiversiy', values_drop_na  = TRUE) %>%
  mutate(alphadiversiy = fct_relevel(alphadiversiy, alpha_measures)) -> alpha_long_df


# alpha_long_df
```


``` r
alpha_long_df %>% 
  phyloseq_explore_alpha(facet_formula = "alphadiversiy ~ Sample",
                     group_by_stats = c("alphadiversiy", "Sample"),
                     stat_formula = "value ~ Time",
                     padjust_method = "fdr",
                     # stat_paired = FALSE,
                     ref_group_stat = "TP1") -> out_sample_time

ls(out_sample_time)
```

```
## [1] "alpha_legend"  "alpha_plot"    "alpha_summary" "anova_result" 
## [5] "anova_table"   "kruskal_test"  "t_test"        "wilcox_test"
```



``` r
out_sample_time$alpha_plot
```

<img src="02_metaphlan_alphadiv_files/figure-html/unnamed-chunk-6-1.png" style="display: block; margin: auto;" />


``` r
out_sample_time$anova %>% 
  DT::datatable()
```

```{=html}
<div class="datatables html-widget html-fill-item" id="htmlwidget-c20bf29149c5e9f458a9" style="width:100%;height:auto;"></div>
<script type="application/json" data-for="htmlwidget-c20bf29149c5e9f458a9">{"x":null,"evals":[],"jsHooks":[]}</script>
```


``` r
out_sample_time$wilcox.test_stat %>% 
  DT::datatable()
```

```{=html}
<div class="datatables html-widget html-fill-item" id="htmlwidget-5c6da86725651da650c6" style="width:100%;height:auto;"></div>
<script type="application/json" data-for="htmlwidget-5c6da86725651da650c6">{"x":null,"evals":[],"jsHooks":[]}</script>
```


``` r
alpha_long_df %>% 
  compute_plot_alpha(facet_formula = "alphadiversiy ~ Sample + cluster_Dtp2",
                     group_by_stats = c("alphadiversiy", "Sample", "cluster_Dtp2"),
                     stat_formula = "value ~ Time",
                     padjust_method = "fdr",
                     # stat_paired = FALSE,
                     ref_group_stat = "TP1") -> out_sample_time_clust


out_sample_time_clust$alpha_plot
```

<img src="02_metaphlan_alphadiv_files/figure-html/unnamed-chunk-9-1.png" style="display: block; margin: auto;" />



``` r
out_sample_time_clust$wilcox.test_stat %>% 
  DT::datatable()
```

```{=html}
<div class="datatables html-widget html-fill-item" id="htmlwidget-5babc318a4b8fb88e497" style="width:100%;height:auto;"></div>
<script type="application/json" data-for="htmlwidget-5babc318a4b8fb88e497">{"x":{"filter":"none","vertical":false,"data":[["1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24","25","26","27","28","29","30","31","32","33","34","35","36","37","38","39","40","41","42","43","44","45","46","47","48"],["diversity_shannon","observed","evenness_pielou","diversity_inverse_simpson","observed","evenness_pielou","diversity_shannon","diversity_inverse_simpson","observed","diversity_shannon","diversity_inverse_simpson","evenness_pielou","diversity_shannon","evenness_pielou","observed","diversity_inverse_simpson","evenness_pielou","observed","diversity_inverse_simpson","evenness_pielou","diversity_inverse_simpson","observed","observed","observed","evenness_pielou","diversity_inverse_simpson","diversity_shannon","evenness_pielou","observed","evenness_pielou","diversity_shannon","diversity_shannon","diversity_inverse_simpson","evenness_pielou","observed","diversity_inverse_simpson","evenness_pielou","diversity_shannon","diversity_shannon","diversity_shannon","diversity_inverse_simpson","diversity_inverse_simpson","diversity_inverse_simpson","diversity_shannon","observed","diversity_shannon","observed","evenness_pielou"],["Saliva","Saliva","Saliva","Saliva","Saliva","Plaque","Plaque","Plaque","Saliva","Saliva","Saliva","Saliva","Saliva","Plaque","Plaque","Saliva","Saliva","Plaque","Saliva","Saliva","Plaque","Plaque","Plaque","Saliva","Plaque","Saliva","Plaque","Saliva","Plaque","Plaque","Plaque","Saliva","Plaque","Saliva","Saliva","Plaque","Plaque","Saliva","Plaque","Plaque","Plaque","Saliva","Plaque","Plaque","Saliva","Saliva","Plaque","Plaque"],["3","3","3","3","2","1","1","1","2","2","2","2","2","2","1","2","1","3","1","1","2","1","2","1","2","3","2","2","2","3","3","1","3","3","1","2","3","3","3","2","3","1","1","1","3","1","3","1"],["TP1","TP1","TP1","TP1","TP1","TP1","TP1","TP1","TP1","TP1","TP1","TP1","TP1","TP1","TP1","TP1","TP1","TP1","TP1","TP1","TP1","TP1","TP1","TP1","TP1","TP1","TP1","TP1","TP1","TP1","TP1","TP1","TP1","TP1","TP1","TP1","TP1","TP1","TP1","TP1","TP1","TP1","TP1","TP1","TP1","TP1","TP1","TP1"],["TP2","TP2","TP2","TP2","TP2","TP2","TP2","TP2","TP3","TP2","TP3","TP2","TP3","TP3","TP2","TP2","TP2","TP2","TP3","TP3","TP3","TP3","TP2","TP2","TP2","TP3","TP3","TP3","TP3","TP2","TP2","TP3","TP3","TP3","TP3","TP2","TP3","TP3","TP3","TP2","TP2","TP2","TP3","TP3","TP3","TP2","TP3","TP3"],[0.003762871355390886,0.009298386763165967,0.009494498176303285,0.0314225799415293,0.03504499384294781,0.03998354586589881,0.06252570958453312,0.07700534759358289,0.08907510441393168,0.09386812441732351,0.09386812441732351,0.1370684953360953,0.1635523947692477,0.2649101369507884,0.2890658693972311,0.3064255196103358,0.3401069518716577,0.3701370109445391,0.3865076100370218,0.3865076100370218,0.4012630681563243,0.4234471410941999,0.4481827511114188,0.4526791810629212,0.4544380131915466,0.4744583445616631,0.4824007518084783,0.4824007518084783,0.4905691312298233,0.4958462923534718,0.5861480768593427,0.6048128342245989,0.6098285373956537,0.6098285373956537,0.6664747017688193,0.667331129756759,0.6832116705619683,0.6832116705619683,0.7596317505530132,0.7688116639994416,0.7856645814048268,0.7961744138214726,0.8148087206910737,0.8883587001234061,0.8903939681448466,0.9314273961332784,0.9862565482025321,1],[0.15,0.15,0.15,0.32,0.32,0.32,0.41,0.41,0.41,0.41,0.41,0.55,0.6,0.79,0.79,0.79,0.79,0.79,0.79,0.79,0.79,0.79,0.79,0.79,0.79,0.79,0.79,0.79,0.79,0.79,0.86,0.86,0.86,0.86,0.86,0.86,0.86,0.86,0.91,0.91,0.91,0.91,0.91,0.95,0.95,0.97,1,1],["Wilcoxon","Wilcoxon","Wilcoxon","Wilcoxon","Wilcoxon","Wilcoxon","Wilcoxon","Wilcoxon","Wilcoxon","Wilcoxon","Wilcoxon","Wilcoxon","Wilcoxon","Wilcoxon","Wilcoxon","Wilcoxon","Wilcoxon","Wilcoxon","Wilcoxon","Wilcoxon","Wilcoxon","Wilcoxon","Wilcoxon","Wilcoxon","Wilcoxon","Wilcoxon","Wilcoxon","Wilcoxon","Wilcoxon","Wilcoxon","Wilcoxon","Wilcoxon","Wilcoxon","Wilcoxon","Wilcoxon","Wilcoxon","Wilcoxon","Wilcoxon","Wilcoxon","Wilcoxon","Wilcoxon","Wilcoxon","Wilcoxon","Wilcoxon","Wilcoxon","Wilcoxon","Wilcoxon","Wilcoxon"],["ns","ns","ns","ns","ns","ns","ns","ns","ns","ns","ns","ns","ns","ns","ns","ns","ns","ns","ns","ns","ns","ns","ns","ns","ns","ns","ns","ns","ns","ns","ns","ns","ns","ns","ns","ns","ns","ns","ns","ns","ns","ns","ns","ns","ns","ns","ns","ns"]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>alphadiversiy<\/th>\n      <th>Sample<\/th>\n      <th>cluster_Dtp2<\/th>\n      <th>group1<\/th>\n      <th>group2<\/th>\n      <th>p<\/th>\n      <th>p.adj<\/th>\n      <th>method<\/th>\n      <th>p.adj.signif<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":[6,7]},{"orderable":false,"targets":0},{"name":" ","targets":0},{"name":"alphadiversiy","targets":1},{"name":"Sample","targets":2},{"name":"cluster_Dtp2","targets":3},{"name":"group1","targets":4},{"name":"group2","targets":5},{"name":"p","targets":6},{"name":"p.adj","targets":7},{"name":"method","targets":8},{"name":"p.adj.signif","targets":9}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
```

# Correlate alpha with metadata:

Then alpha diversity at baseline associated/predicts with bleeding + plaque at TP2 TP3:


``` r
alphas %>%
  correlate_alpha(colnums_to_plot = c("diversity_shannon", "diversity_inverse_simpson", "evenness_pielou",
                                      "mean_bleeding", "mean_plaque", "Sample", "Time"), colour =  "Time",
                  method = "spearman") -> alpha_corr

alpha_corr
```

<img src="02_metaphlan_alphadiv_files/figure-html/unnamed-chunk-11-1.png" style="display: block; margin: auto;" />



``` r
save(out_sample_time_clust,
     out_sample_time,
     alpha_corr,
     ps_up,
     sample_pal,
     time_pal, sub_pal, sex_pal, sample_pal,
     file = here::here("../../data/processed_data/metaphlan/02_data.Rdata"))
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
##  [1] rstatix_0.7.2        ampvis2_2.8.9        ggpubr_0.6.0        
##  [4] microbiome_1.28.0    metagMisc_0.5.0      reshape2_1.4.4      
##  [7] scales_1.3.0         microViz_0.12.4      speedyseq_0.5.3.9021
## [10] phyloseq_1.50.0      readxl_1.4.3         lubridate_1.9.3     
## [13] forcats_1.0.0        stringr_1.5.1        dplyr_1.1.4         
## [16] purrr_1.0.2          readr_2.1.5          tidyr_1.3.1         
## [19] tibble_3.2.1         ggplot2_3.5.1        tidyverse_2.0.0     
## 
## loaded via a namespace (and not attached):
##   [1] permute_0.9-7           rlang_1.1.4             magrittr_2.0.3         
##   [4] ade4_1.7-22             compiler_4.4.0          mgcv_1.9-1             
##   [7] vctrs_0.6.5             pkgconfig_2.0.3         crayon_1.5.3           
##  [10] fastmap_1.2.0           backports_1.5.0         XVector_0.45.0         
##  [13] labeling_0.4.3          utf8_1.2.4              rmarkdown_2.29         
##  [16] tzdb_0.4.0              UCSC.utils_1.2.0        xfun_0.49              
##  [19] zlibbioc_1.51.2         cachem_1.1.0            GenomeInfoDb_1.42.0    
##  [22] jsonlite_1.8.9          biomformat_1.34.0       highr_0.11             
##  [25] rhdf5filters_1.17.0     Rhdf5lib_1.27.0         broom_1.0.7            
##  [28] parallel_4.4.0          cluster_2.1.6           R6_2.5.1               
##  [31] RColorBrewer_1.1-3      bslib_0.8.0             stringi_1.8.4          
##  [34] GGally_2.2.1            car_3.1-3               jquerylib_0.1.4        
##  [37] cellranger_1.1.0        Rcpp_1.0.13-1           iterators_1.0.14       
##  [40] knitr_1.48              IRanges_2.39.2          Matrix_1.7-1           
##  [43] splines_4.4.0           igraph_2.1.1            timechange_0.3.0       
##  [46] tidyselect_1.2.1        rstudioapi_0.17.1       abind_1.4-8            
##  [49] yaml_2.3.10             vegan_2.6-8             codetools_0.2-20       
##  [52] lattice_0.22-6          plyr_1.8.9              Biobase_2.66.0         
##  [55] withr_3.0.2             evaluate_1.0.1          Rtsne_0.17             
##  [58] survival_3.7-0          ggstats_0.7.0           Biostrings_2.73.2      
##  [61] pillar_1.9.0            carData_3.0-5           DT_0.33                
##  [64] foreach_1.5.2           stats4_4.4.0            plotly_4.10.4          
##  [67] generics_0.1.3          rprojroot_2.0.4         S4Vectors_0.43.2       
##  [70] hms_1.1.3               munsell_0.5.1           glue_1.8.0             
##  [73] lazyeval_0.2.2          tools_4.4.0             data.table_1.16.2      
##  [76] ggsignif_0.6.4          cowplot_1.1.3           rhdf5_2.49.0           
##  [79] grid_4.4.0              ape_5.8                 crosstalk_1.2.1        
##  [82] colorspace_2.1-1        nlme_3.1-166            GenomeInfoDbData_1.2.13
##  [85] Formula_1.2-5           cli_3.6.3               fansi_1.0.6            
##  [88] viridisLite_0.4.2       gtable_0.3.6            sass_0.4.9             
##  [91] digest_0.6.37           BiocGenerics_0.52.0     ggrepel_0.9.6          
##  [94] farver_2.1.2            htmlwidgets_1.6.4       htmltools_0.5.8.1      
##  [97] multtest_2.61.0         lifecycle_1.0.4         httr_1.4.7             
## [100] here_1.0.1              MASS_7.3-61
```
