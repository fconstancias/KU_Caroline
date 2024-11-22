---
title: " Study II - Beta-diversity analyses of metaphlan profiles"
author: "Florentin Constancias"
date: "November 20, 2024"
output: 
  html_document: 
    toc: yes
    keep_md: yes
---













``` r
load(here::here("../../data/processed_data/metaphlan/01_data.Rdata"))
load(here::here("../../data/processed_data/metaphlan/02_data.Rdata"))
```


``` r
pd <- position_dodge(0.3)
```

Let's source the function:


``` r
source("https://raw.githubusercontent.com/fconstancias/KU_Caroline/refs/heads/main/code/functions/phyloseq_functions.R")

# details:

#' @title This function calculates beta diversity, generates ordination plots (e.g., PCoA),  performs PERMANOVA, and optionally creates boxplots of distances between metadata groups.
#' @author Florentin Constancias
#' @param ps_up A phyloseq object containing microbiome data.
#' @param beta A list of distance matrices for beta diversity (e.g., Bray-Curtis, Aitchison).
#' @param m Character string specifying the ordination method (default: "PCoA").
#' @param color_group Character string for the variable used to color points in ordination plots (default: "Time").
#' @param shape_group Character string for the variable used to shape points in ordination plots (default: "Sample").
#' @param alpha Optional variable for transparency adjustment in ordination plots (default: NULL).
#' @param col_pal Color palette for point colors (default: `time_pal`).
#' @param fill_pal Fill color palette for point fills (default: `time_pal`).
#' @param path_group Character string for the variable grouping paths in ordination plots (default: "interaction(Sample,Subject)").
#' @param facet_formula Formula for faceting ordination plots (default: "Sample ~ .").
#' @param axis1 Integer specifying the ordination axis to plot on the x-axis (default: 1).
#' @param axis2 Integer specifying the ordination axis to plot on the y-axis (default: 2).
#' @param seed Integer for random seed, ensuring reproducible results in ordination (default: 123).
#' @param permanova_terms Character vector of variables for PERMANOVA tests (default: c("Time", "cluster_Dtp2")).
#' @param metadata_dist_boxplot Optional variable for generating distance boxplots (default: NULL).
#' @param strata Variable for stratified PERMANOVA testing (default: "none").
#' @param perm Number of permutations for PERMANOVA (default: 999).
#' 
#' @return A list of results, including:
#'   - `pcoas`: Ordination plots with faceting and custom color/shape groups.
#'   - `expl_var`: Explained variance for each axis in ordination plots.
#'   - `PCOA`: Detailed ordination plot for the rAitchison distance with paths.
#'   - `PCOA_leg`: Legend for the ordination plot.
#'   - `perm`: PERMANOVA results.
#'   - `pw_perm`: Pairwise PERMANOVA results for each variable in `permanova_terms`.
#'   - `tw_perm`: Two-way PERMANOVA results for each variable in `permanova_terms`.
#'   - Optional: `dist_boxplot` if `metadata_dist_boxplot` is provided, showing boxplots of distances.
```


# Compute the distance for all samples:




# Run fonction for full dataset:




Explore PCoA : distances show similar patterns


``` r
out$pcoas
```

<img src="03_metaphlan_betaadiv_files/figure-html/unnamed-chunk-6-1.png" style="display: block; margin: auto;" />

Explained variance of the plots above.


``` r
out$expl_var %>% 
  mutate_if(is.numeric, ~round(., 3)) %>% 
  DT::datatable()
```

```{=html}
<div class="datatables html-widget html-fill-item" id="htmlwidget-9eceb20dc2a931752ae4" style="width:100%;height:auto;"></div>
<script type="application/json" data-for="htmlwidget-9eceb20dc2a931752ae4">{"x":{"filter":"none","vertical":false,"data":[["1","2","3"],["bjaccard","wjaccard","rAitchison"],["Axis.1   [21%]","Axis.1   [19.9%]","Axis.1   [20%]"],["bjaccard","wjaccard","rAitchison"],["Axis.2   [9.9%]","Axis.2   [5.9%]","Axis.2   [3.9%]"]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>.id...1<\/th>\n      <th>V1...2<\/th>\n      <th>.id...3<\/th>\n      <th>V1...4<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"orderable":false,"targets":0},{"name":" ","targets":0},{"name":".id...1","targets":1},{"name":"V1...2","targets":2},{"name":".id...3","targets":3},{"name":"V1...4","targets":4}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
```

Focus on Robust Aitchinson distance, clear discrimination between site, display trajectories.


``` r
out$PCOA + facet_null()
```

<img src="03_metaphlan_betaadiv_files/figure-html/unnamed-chunk-8-1.png" style="display: block; margin: auto;" />
Corresponding legend:


``` r
out$PCOA_leg
```

<img src="03_metaphlan_betaadiv_files/figure-html/unnamed-chunk-9-1.png" style="display: block; margin: auto;" />

Statistical evaluation of the beta diversity: association with Sample_Type, Time and Cluster of inflamation at TP2 -PERMANOVA:

Sample Type explain most of the variation (R2) > Time > Cluster

More info: https://gist.github.com/claczny/3415270a6c6919969bff79d2c246a527



``` r
out$perm %>% 
  mutate_if(is.double, ~round(., 3)) %>% 
  dplyr::filter(Distance == "rAitchison", terms_margin == "terms") %>% 
  DT::datatable()
```

```{=html}
<div class="datatables html-widget html-fill-item" id="htmlwidget-f30d7598595affc4f25a" style="width:100%;height:auto;"></div>
<script type="application/json" data-for="htmlwidget-f30d7598595affc4f25a">{"x":{"filter":"none","vertical":false,"data":[["1","2","3","4","5","6","7","8","9","10","11","12","13","14"],["rAitchison","rAitchison","rAitchison","rAitchison","rAitchison","rAitchison","rAitchison","rAitchison","rAitchison","rAitchison","rAitchison","rAitchison","rAitchison","rAitchison"],["Sample","Time","cluster_Dtp2","Sample:Time","Sample:cluster_Dtp2","Time:cluster_Dtp2","Sample:Time:cluster_Dtp2","Residual","Total","Sample","Time","cluster_Dtp2","Residual","Total"],[1,2,2,2,2,4,4,221,238,1,2,2,233,238],[3.11,0.318,0.175,0.263,0.142,0.214,0.165,12.652,17.039,3.11,0.318,0.175,13.436,17.039],[0.183,0.019,0.01,0.015,0.008,0.013,0.01,0.743,1,0.183,0.019,0.01,0.789,1],[54.326,2.778,1.527,2.298,1.236,0.9330000000000001,0.722,null,null,53.935,2.758,1.516,null,null],[0.001,0.001,0.027,0.001,0.129,0.605,0.998,null,null,0.001,0.001,0.028,null,null],["Sample * Time * cluster_Dtp2","Sample * Time * cluster_Dtp2","Sample * Time * cluster_Dtp2","Sample * Time * cluster_Dtp2","Sample * Time * cluster_Dtp2","Sample * Time * cluster_Dtp2","Sample * Time * cluster_Dtp2","Sample * Time * cluster_Dtp2","Sample * Time * cluster_Dtp2","Sample + Time + cluster_Dtp2","Sample + Time + cluster_Dtp2","Sample + Time + cluster_Dtp2","Sample + Time + cluster_Dtp2","Sample + Time + cluster_Dtp2"],["terms","terms","terms","terms","terms","terms","terms","terms","terms","terms","terms","terms","terms","terms"]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>Distance<\/th>\n      <th>terms<\/th>\n      <th>Df<\/th>\n      <th>SumOfSqs<\/th>\n      <th>R2<\/th>\n      <th>F<\/th>\n      <th>Pr(&gt;F)<\/th>\n      <th>Group<\/th>\n      <th>terms_margin<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":[3,4,5,6,7]},{"orderable":false,"targets":0},{"name":" ","targets":0},{"name":"Distance","targets":1},{"name":"terms","targets":2},{"name":"Df","targets":3},{"name":"SumOfSqs","targets":4},{"name":"R2","targets":5},{"name":"F","targets":6},{"name":"Pr(>F)","targets":7},{"name":"Group","targets":8},{"name":"terms_margin","targets":9}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
```

Statistical evaluation of the beta diversity - Pairwise PERMANOVA tests between Time, Sample Type, Inflamation response clusters:



``` r
out$pw_perm %>% 
  dplyr::filter(Distance == "rAitchison") %>% 
  mutate_if(is.numeric, ~round(., 3)) %>% 
  DT::datatable()
```

```{=html}
<div class="datatables html-widget html-fill-item" id="htmlwidget-b85bf888dc62471906a9" style="width:100%;height:auto;"></div>
<script type="application/json" data-for="htmlwidget-b85bf888dc62471906a9">{"x":{"filter":"none","vertical":false,"data":[["1","2","3","4","5","6","7"],["rAitchison","rAitchison","rAitchison","rAitchison","rAitchison","rAitchison","rAitchison"],["1","1","2","TP1","TP1","TP2","Saliva"],["2","3","3","TP2","TP3","TP3","Plaque"],[0.008999999999999999,0.007,0.007,0.018,0.008,0.017,0.183],[0.148,0.214,0.134,0.001,0.155,0.002,0.001],[0.444,0.642,0.402,0.003,0.465,0.006,0.001],[0.222,0.214,0.402,0.003,0.155,0.003,0.001],["cluster_Dtp2","cluster_Dtp2","cluster_Dtp2","Time","Time","Time","Sample"]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>Distance<\/th>\n      <th>X1<\/th>\n      <th>X2<\/th>\n      <th>R2<\/th>\n      <th>pval<\/th>\n      <th>pvalBon<\/th>\n      <th>pvalFDR<\/th>\n      <th>Group<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":[4,5,6,7]},{"orderable":false,"targets":0},{"name":" ","targets":0},{"name":"Distance","targets":1},{"name":"X1","targets":2},{"name":"X2","targets":3},{"name":"R2","targets":4},{"name":"pval","targets":5},{"name":"pvalBon","targets":6},{"name":"pvalFDR","targets":7},{"name":"Group","targets":8}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
```

We could also play the intra-subject distance accross time | maybe with within group variation for TP1: within TP1 / Sample:


``` r
out$dist_box
```

<img src="03_metaphlan_betaadiv_files/figure-html/unnamed-chunk-12-1.png" style="display: block; margin: auto;" />


Statistical evaluation of the beta diversity - Pairwise Multivariate T.tests between Time, Sample Type, Inflamation response clusters:


``` r
# out$tw_perm %>% 
#   DT::datatable()
```

We can also explore the groups of the different variables:

Sample


``` r
out$Sample$Saliva
```

<img src="03_metaphlan_betaadiv_files/figure-html/unnamed-chunk-14-1.png" style="display: block; margin: auto;" />


``` r
out$Sample$Plaque
```

<img src="03_metaphlan_betaadiv_files/figure-html/unnamed-chunk-15-1.png" style="display: block; margin: auto;" />

``` r
out$Time$TP1
```

<img src="03_metaphlan_betaadiv_files/figure-html/unnamed-chunk-16-1.png" style="display: block; margin: auto;" />


``` r
out$Time$TP2
```

<img src="03_metaphlan_betaadiv_files/figure-html/unnamed-chunk-17-1.png" style="display: block; margin: auto;" />


``` r
out$Time$TP3
```

<img src="03_metaphlan_betaadiv_files/figure-html/unnamed-chunk-18-1.png" style="display: block; margin: auto;" />


# Now we can explore more Saliva and Plaque samples:


## Saliva:


``` r
ps_up %>%
  subset_samples(Sample == "Saliva") %>% 
  subset_taxa(Class != "UNCLASSIFIED") %>% 
  transform_sample_counts(function(x) x/sum(x) * 1) %>% 
  microViz::ps_mutate(cluster_Dtp2 = as.factor(cluster_Dtp2)) %>% 
  compute_plot_beta(ps_up = .,
                    beta = beta,
                    color_group = "Time",
                    shape_group = "Sample",
                    alpha = NULL,
                    col_pal = time_pal,
                    fill_pal = time_pal,
                    path_group = "interaction(Subject)",
                    facet_formula = "Sample ~ .",
                    axis1 = 1,
                    axis2 = 2,
                    seed = 123,
                    permanova_terms = c("Time","cluster_Dtp2"),
                    # metadata_dist_boxplot = NULL ,
                    strata = "none") -> saliva
```



``` r
saliva$pcoas
```

<img src="03_metaphlan_betaadiv_files/figure-html/unnamed-chunk-20-1.png" style="display: block; margin: auto;" />



``` r
saliva$PCOA
```

<img src="03_metaphlan_betaadiv_files/figure-html/unnamed-chunk-21-1.png" style="display: block; margin: auto;" />


``` r
saliva$perm %>% 
  mutate_if(is.double, ~round(., 3)) %>% 
  dplyr::filter(Distance == "rAitchison", terms_margin == "terms") %>%
  DT::datatable()
```

```{=html}
<div class="datatables html-widget html-fill-item" id="htmlwidget-88fdd2a34d502a66744b" style="width:100%;height:auto;"></div>
<script type="application/json" data-for="htmlwidget-88fdd2a34d502a66744b">{"x":{"filter":"none","vertical":false,"data":[["1","2","3","4","5","6","7","8","9"],["rAitchison","rAitchison","rAitchison","rAitchison","rAitchison","rAitchison","rAitchison","rAitchison","rAitchison"],["Time","cluster_Dtp2","Time:cluster_Dtp2","Residual","Total","Time","cluster_Dtp2","Residual","Total"],[2,2,4,111,119,2,2,115,119],[0.335,0.15,0.161,5.942,6.589,0.335,0.15,6.103,6.589],[0.051,0.023,0.024,0.902,1,0.051,0.023,0.926,1],[3.133,1.402,0.754,null,null,3.16,1.414,null,null],[0.001,0.002,1,null,null,0.001,0.003,null,null],["Time * cluster_Dtp2","Time * cluster_Dtp2","Time * cluster_Dtp2","Time * cluster_Dtp2","Time * cluster_Dtp2","Time + cluster_Dtp2","Time + cluster_Dtp2","Time + cluster_Dtp2","Time + cluster_Dtp2"],["terms","terms","terms","terms","terms","terms","terms","terms","terms"]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>Distance<\/th>\n      <th>terms<\/th>\n      <th>Df<\/th>\n      <th>SumOfSqs<\/th>\n      <th>R2<\/th>\n      <th>F<\/th>\n      <th>Pr(&gt;F)<\/th>\n      <th>Group<\/th>\n      <th>terms_margin<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":[3,4,5,6,7]},{"orderable":false,"targets":0},{"name":" ","targets":0},{"name":"Distance","targets":1},{"name":"terms","targets":2},{"name":"Df","targets":3},{"name":"SumOfSqs","targets":4},{"name":"R2","targets":5},{"name":"F","targets":6},{"name":"Pr(>F)","targets":7},{"name":"Group","targets":8},{"name":"terms_margin","targets":9}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
```


Pairwise comparisons between groups: Time, ... 


Interesting to see that when considering only presence and absence of taxa (i.e., binary bJaccard), no significant differences between TP1 and TP3, but still when considering proportion of taxa. 



``` r
saliva$pw_perm %>% 
  # dplyr::filter(Distance == "rAitchison") %>% 
  mutate_if(is.numeric, ~round(., 3)) %>% 
  DT::datatable()
```

```{=html}
<div class="datatables html-widget html-fill-item" id="htmlwidget-8e1fe96c0688a849169a" style="width:100%;height:auto;"></div>
<script type="application/json" data-for="htmlwidget-8e1fe96c0688a849169a">{"x":{"filter":"none","vertical":false,"data":[["1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18"],["bjaccard","bjaccard","bjaccard","wjaccard","wjaccard","wjaccard","rAitchison","rAitchison","rAitchison","bjaccard","bjaccard","bjaccard","wjaccard","wjaccard","wjaccard","rAitchison","rAitchison","rAitchison"],["1","1","2","1","1","2","1","1","2","TP1","TP1","TP2","TP1","TP1","TP2","TP1","TP1","TP2"],["2","3","3","2","3","3","2","3","3","TP2","TP3","TP3","TP2","TP3","TP3","TP2","TP3","TP3"],[0.034,0.025,0.02,0.023,0.019,0.013,0.023,0.018,0.013,0.073,0.018,0.074,0.04,0.043,0.042,0.05,0.025,0.04],[0.003,0.016,0.005,0.032,0.041,0.238,0.003,0.022,0.095,0.001,0.064,0.001,0.001,0.001,0.001,0.001,0.001,0.001],[0.008999999999999999,0.048,0.015,0.096,0.123,0.714,0.008999999999999999,0.066,0.285,0.003,0.192,0.003,0.003,0.003,0.003,0.003,0.003,0.003],[0.008999999999999999,0.016,0.007,0.096,0.062,0.238,0.008999999999999999,0.033,0.095,0.002,0.064,0.002,0.002,0.002,0.002,0.002,0.002,0.002],["cluster_Dtp2","cluster_Dtp2","cluster_Dtp2","cluster_Dtp2","cluster_Dtp2","cluster_Dtp2","cluster_Dtp2","cluster_Dtp2","cluster_Dtp2","Time","Time","Time","Time","Time","Time","Time","Time","Time"]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>Distance<\/th>\n      <th>X1<\/th>\n      <th>X2<\/th>\n      <th>R2<\/th>\n      <th>pval<\/th>\n      <th>pvalBon<\/th>\n      <th>pvalFDR<\/th>\n      <th>Group<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":[4,5,6,7]},{"orderable":false,"targets":0},{"name":" ","targets":0},{"name":"Distance","targets":1},{"name":"X1","targets":2},{"name":"X2","targets":3},{"name":"R2","targets":4},{"name":"pval","targets":5},{"name":"pvalBon","targets":6},{"name":"pvalFDR","targets":7},{"name":"Group","targets":8}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
```


``` r
# saliva$tw_perm %>% 
#   DT::datatable()
```

We could test if the distance between saliva and plaque per individual changes over time:

We have this distance to the baseline boxplot - which is the same as prev. but only Saliva here. But for this we need to use the ouputs generated with the full dataset.


``` r
saliva$dist_box
```

```
## NULL
```

During perturbation : TP2 : saliva and Plaque microbiome are similar:


``` r
out$dist_df %>% 
  dplyr::filter(Sample_1 != Sample_2,
                Time_1 == Time_2,
                Subject_1 == Subject_2) %>%
  arrange(Distance, Subject_2, Time_2) %>% 
  ggplot(data = ., aes_string(x="Time_2", y="value")) +
  geom_boxplot(outlier.colour = NA, alpha=0.7, aes_string(fill = "Time_2")) +
  geom_jitter(size=1, position = pd, aes_string(color = "Subject_1")) +
  geom_line(aes_string(group="Subject_1"), position = pd, linetype = "dashed", color = "grey50", linewidth = 0.08) +
  facet_grid(as.formula(paste0("Distance ~ Sample_1")), scales = "free_y", space = "fixed", switch = "y") +
  scale_color_manual(name = "", values = sub_pal,
                     na.value = "black") +
  scale_fill_manual(name = "", values = time_pal,
                    na.value = "black") +
  theme_light() + ylab("Between Site Intra-individual distance ") + xlab(NULL) + theme(
    axis.text.x = element_blank()) +  theme_linedraw() + theme(strip.placement = "outside") -> d_box

d_box %>%
  ggpubr::get_legend(.) %>%
  ggpubr::as_ggplot(.) -> d_box_leg

d_box + theme(legend.position = "none") -> d_box

d_box
```

<img src="03_metaphlan_betaadiv_files/figure-html/unnamed-chunk-26-1.png" style="display: block; margin: auto;" />

## Plaque:


``` r
ps_up %>%
  subset_samples(Sample == "Plaque") %>% 
  subset_taxa(Class != "UNCLASSIFIED") %>% 
  transform_sample_counts(function(x) x/sum(x) * 1) %>% 
  microViz::ps_mutate(cluster_Dtp2 = as.factor(cluster_Dtp2)) %>% 
  compute_plot_beta(ps_up = .,
                    beta = beta,
                    color_group = "Time",
                    shape_group = "Sample",
                    alpha = NULL,
                    col_pal = time_pal,
                    fill_pal = time_pal,
                    path_group = "interaction(Subject)",
                    facet_formula = "Sample ~ .",
                    axis1 = 1,
                    axis2 = 2,
                    seed = 123,
                    permanova_terms = c("Time","cluster_Dtp2"),
                    # metadata_dist_boxplot = NULL ,
                    strata = "none") -> plaque
```



``` r
plaque$pcoas
```

<img src="03_metaphlan_betaadiv_files/figure-html/unnamed-chunk-28-1.png" style="display: block; margin: auto;" />



``` r
plaque$PCOA
```

<img src="03_metaphlan_betaadiv_files/figure-html/unnamed-chunk-29-1.png" style="display: block; margin: auto;" />


``` r
plaque$perm %>% 
  mutate_if(is.double, ~round(., 3)) %>% 
  dplyr::filter(Distance == "rAitchison", terms_margin == "terms") %>%
  DT::datatable()
```

```{=html}
<div class="datatables html-widget html-fill-item" id="htmlwidget-bb3486349adfa5d821dc" style="width:100%;height:auto;"></div>
<script type="application/json" data-for="htmlwidget-bb3486349adfa5d821dc">{"x":{"filter":"none","vertical":false,"data":[["1","2","3","4","5","6","7","8","9"],["rAitchison","rAitchison","rAitchison","rAitchison","rAitchison","rAitchison","rAitchison","rAitchison","rAitchison"],["Time","cluster_Dtp2","Time:cluster_Dtp2","Residual","Total","Time","cluster_Dtp2","Residual","Total"],[2,2,4,110,118,2,2,114,118],[0.246,0.166,0.218,6.71,7.34,0.246,0.166,6.928,7.34],[0.034,0.023,0.03,0.914,1,0.034,0.023,0.944,1],[2.016,1.362,0.892,null,null,2.024,1.367,null,null],[0.001,0.01,0.888,null,null,0.001,0.01,null,null],["Time * cluster_Dtp2","Time * cluster_Dtp2","Time * cluster_Dtp2","Time * cluster_Dtp2","Time * cluster_Dtp2","Time + cluster_Dtp2","Time + cluster_Dtp2","Time + cluster_Dtp2","Time + cluster_Dtp2"],["terms","terms","terms","terms","terms","terms","terms","terms","terms"]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>Distance<\/th>\n      <th>terms<\/th>\n      <th>Df<\/th>\n      <th>SumOfSqs<\/th>\n      <th>R2<\/th>\n      <th>F<\/th>\n      <th>Pr(&gt;F)<\/th>\n      <th>Group<\/th>\n      <th>terms_margin<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":[3,4,5,6,7]},{"orderable":false,"targets":0},{"name":" ","targets":0},{"name":"Distance","targets":1},{"name":"terms","targets":2},{"name":"Df","targets":3},{"name":"SumOfSqs","targets":4},{"name":"R2","targets":5},{"name":"F","targets":6},{"name":"Pr(>F)","targets":7},{"name":"Group","targets":8},{"name":"terms_margin","targets":9}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
```

Pairwise comparisons between groups: Time, ... 

Interesting to compare between saliva and plaque... Maybe a plot displaying R2 and adjpavalue?



``` r
plaque$pw_perm %>% 
  # dplyr::filter(Distance == "rAitchison") %>% 
  mutate_if(is.numeric, ~round(., 3)) %>% 
  DT::datatable()
```

```{=html}
<div class="datatables html-widget html-fill-item" id="htmlwidget-44474b3e8cd160f85a98" style="width:100%;height:auto;"></div>
<script type="application/json" data-for="htmlwidget-44474b3e8cd160f85a98">{"x":{"filter":"none","vertical":false,"data":[["1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18"],["bjaccard","bjaccard","bjaccard","wjaccard","wjaccard","wjaccard","rAitchison","rAitchison","rAitchison","bjaccard","bjaccard","bjaccard","wjaccard","wjaccard","wjaccard","rAitchison","rAitchison","rAitchison"],["1","1","2","1","1","2","1","1","2","TP1","TP1","TP2","TP1","TP1","TP2","TP1","TP1","TP2"],["2","3","3","2","3","3","2","3","3","TP2","TP3","TP3","TP2","TP3","TP3","TP2","TP3","TP3"],[0.027,0.022,0.025,0.023,0.021,0.018,0.018,0.016,0.017,0.041,0.008,0.036,0.038,0.017,0.029,0.033,0.014,0.028],[0.019,0.03,0.001,0.041,0.033,0.019,0.092,0.133,0.008,0.001,0.981,0.001,0.001,0.08500000000000001,0.003,0.001,0.287,0.001],[0.057,0.09,0.003,0.123,0.099,0.057,0.276,0.399,0.024,0.003,2.943,0.003,0.003,0.255,0.008999999999999999,0.003,0.861,0.003],[0.028,0.03,0.003,0.041,0.05,0.057,0.138,0.133,0.024,0.002,0.981,0.002,0.003,0.08500000000000001,0.005,0.002,0.287,0.002],["cluster_Dtp2","cluster_Dtp2","cluster_Dtp2","cluster_Dtp2","cluster_Dtp2","cluster_Dtp2","cluster_Dtp2","cluster_Dtp2","cluster_Dtp2","Time","Time","Time","Time","Time","Time","Time","Time","Time"]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>Distance<\/th>\n      <th>X1<\/th>\n      <th>X2<\/th>\n      <th>R2<\/th>\n      <th>pval<\/th>\n      <th>pvalBon<\/th>\n      <th>pvalFDR<\/th>\n      <th>Group<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":[4,5,6,7]},{"orderable":false,"targets":0},{"name":" ","targets":0},{"name":"Distance","targets":1},{"name":"X1","targets":2},{"name":"X2","targets":3},{"name":"R2","targets":4},{"name":"pval","targets":5},{"name":"pvalBon","targets":6},{"name":"pvalFDR","targets":7},{"name":"Group","targets":8}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
```


``` r
plaque$tw_perm %>% 
  DT::datatable()
```

```{=html}
<div class="datatables html-widget html-fill-item" id="htmlwidget-2514ce28ae593bde0f73" style="width:100%;height:auto;"></div>
<script type="application/json" data-for="htmlwidget-2514ce28ae593bde0f73">{"x":{"filter":"none","vertical":false,"data":[["1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18"],["bjaccard","bjaccard","bjaccard","wjaccard","wjaccard","wjaccard","rAitchison","rAitchison","rAitchison","bjaccard","bjaccard","bjaccard","wjaccard","wjaccard","wjaccard","rAitchison","rAitchison","rAitchison"],[["1"],["1"],["2"],["1"],["1"],["2"],["1"],["1"],["2"],["TP1"],["TP1"],["TP2"],["TP1"],["TP1"],["TP2"],["TP1"],["TP1"],["TP2"]],[["2"],["3"],["3"],["2"],["3"],["3"],["2"],["3"],["3"],["TP2"],["TP3"],["TP3"],["TP2"],["TP3"],["TP3"],["TP2"],["TP3"],["TP3"]],[[26],[26],[42],[26],[26],[42],[26],[26],[42],[40],[40],[40],[40],[40],[40],[40],[40],[40]],[[42],[51],[51],[42],[51],[51],[42],[51],[51],[40],[39],[39],[40],[39],[39],[40],[39],[39]],[[0.032],[0.04],[0.003],[0.047],[0.037],[0.016],[0.097],[0.128],[0.003],[0.001],[0.983],[0.001],[0.001],[0.08500000000000001],[0.001],[0.002],[0.275],[0.002]],[[1.73461418947618],[1.634795735512252],[2.365034924791286],[1.495570207723217],[1.548623281171168],[1.671255386679231],[1.235564436552785],[1.19799955626987],[1.555500854485347],[3.35603097867184],[0.5892193271707984],[2.880293118537155],[3.094735538500077],[1.31930578439638],[2.28195584873116],[2.697926297564111],[1.064331802150514],[2.240378589349513]],[[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999],[999]],[0.04,0.04,0.009000000000000001,0.047,0.047,0.047,0.128,0.128,0.009000000000000001,0.0015,0.983,0.0015,0.0015,0.08500000000000001,0.0015,0.003,0.275,0.003],["cluster_Dtp2","cluster_Dtp2","cluster_Dtp2","cluster_Dtp2","cluster_Dtp2","cluster_Dtp2","cluster_Dtp2","cluster_Dtp2","cluster_Dtp2","Time","Time","Time","Time","Time","Time","Time","Time","Time"]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>Distance<\/th>\n      <th>Level1<\/th>\n      <th>Level2<\/th>\n      <th>N1<\/th>\n      <th>N2<\/th>\n      <th>p.value<\/th>\n      <th>tw2.stat<\/th>\n      <th>nrep<\/th>\n      <th>pvalFDR<\/th>\n      <th>Group<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":9},{"orderable":false,"targets":0},{"name":" ","targets":0},{"name":"Distance","targets":1},{"name":"Level1","targets":2},{"name":"Level2","targets":3},{"name":"N1","targets":4},{"name":"N2","targets":5},{"name":"p.value","targets":6},{"name":"tw2.stat","targets":7},{"name":"nrep","targets":8},{"name":"pvalFDR","targets":9},{"name":"Group","targets":10}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
```

Display PERMANOVA results as a plot - in progress:


``` r
permanova_terms = c("Time","cluster_Dtp2")


saliva$perm %>% 
  mutate(Data = "Saliva") %>% 
  rbind(.,
        plaque$perm %>% 
          mutate(Data = "Plaque")    
  ) %>% 
  filter(!terms %in% c("Residual","Total")) %>%
  mutate(terms = fct_relevel(terms, c(permanova_terms, "Time:cluster_Dtp2"))) %>%
  # arrange(Distance, Group2) %>%
  # filter(p.adj <= 0.05) %>%
  # mutate(p.adj = replace(p.adj, p.adj > 0.05, NA)) %>% 
  ggplot(mapping = aes(x = terms,
                       y = R2, #`Pr(>F)` ,
                       size = R2,
                       # fill = `Pr(>F)`,
                       color = `Pr(>F)`)) + #,
  # size = R2)) +
  facet_grid(Distance + Data ~ terms_margin+ Group,  scales = "fixed", space = "fixed", drop = FALSE) +
  geom_point(alpha=0.8) +
  # geom_tile() +
  # geom_bar(stat="identity") +
  theme_linedraw() + 
  xlab(NULL) + ylab(NULL) +
  scale_size(range = c(1, 8), 
             breaks = c(0.05, 0.1,0.2,0.6),
             labels = c(0.05, 0.1,0.2,0.6),
             limits = c(0,1),
             # trans = "sqrt",
             name="PERMANOVA R2 value") + 
  scale_color_viridis_c(name = "adj p.value",
                        na.value = "grey10",
                        trans = scales::pseudo_log_trans(sigma = 0.001),
                        limits = c(0,0.05),
                        breaks = c(1,0.05, 0.01, 0.001),
                        labels = c(1,0.05, 0.01, 0.001)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) #+
```

<img src="03_metaphlan_betaadiv_files/figure-html/unnamed-chunk-33-1.png" style="display: block; margin: auto;" />

``` r
# coord_equal()
```

Site across time?


``` r
tmp = ps_up
group = "Time"
Time_out=NULL
# Time_out <- vector("list", length(tmp %>%
#                                get_variable(group) %>%
#                                levels()))
# names(Time_out) <- tmp %>%
#   get_variable(group) %>% levels()

for(tp in tmp %>%
    get_variable(group) %>%
    unique()){
  
  prune_samples(get_variable(tmp, group) == tp,
                tmp) %>%
    subset_taxa(Class != "UNCLASSIFIED") %>% 
    transform_sample_counts(function(x) x/sum(x) * 1) %>% 
    microViz::ps_mutate(cluster_Dtp2 = as.factor(cluster_Dtp2)) %>% 
    compute_plot_beta(ps_up = .,
                      beta = beta,
                      color_group = "Sample",
                      shape_group = "Sample",
                      alpha = NULL,
                      col_pal = sample_pal,
                      fill_pal = sample_pal,
                      path_group = "interaction(Subject)",
                      facet_formula = "Sample ~ .",
                      axis1 = 1,
                      axis2 = 2,
                      seed = 123,
                      permanova_terms = c("Sample"),
                      metadata_dist_boxplot = NULL ,
                      strata = "none") -> Time_out[[tp]]
}
```



``` r
Time_out %>% 
  plyr::ldply(., function(x) x$perm) %>%
  rename(Time = '.id') %>% 
  DT::datatable()
```

```{=html}
<div class="datatables html-widget html-fill-item" id="htmlwidget-5d05eac75555269eefd0" style="width:100%;height:auto;"></div>
<script type="application/json" data-for="htmlwidget-5d05eac75555269eefd0">{"x":{"filter":"none","vertical":false,"data":[["1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24","25","26","27","28","29","30","31","32","33","34","35","36","37","38","39","40","41","42","43","44","45","46","47","48","49","50","51","52","53","54","55","56","57","58","59","60","61","62","63","64","65","66","67","68","69","70","71","72","73","74","75","76","77","78","79","80","81","82","83","84","85","86","87","88","89","90","91","92","93","94","95","96","97","98","99","100","101","102","103","104","105","106","107","108"],["TP1","TP1","TP1","TP1","TP1","TP1","TP1","TP1","TP1","TP1","TP1","TP1","TP1","TP1","TP1","TP1","TP1","TP1","TP1","TP1","TP1","TP1","TP1","TP1","TP1","TP1","TP1","TP1","TP1","TP1","TP1","TP1","TP1","TP1","TP1","TP1","TP2","TP2","TP2","TP2","TP2","TP2","TP2","TP2","TP2","TP2","TP2","TP2","TP2","TP2","TP2","TP2","TP2","TP2","TP2","TP2","TP2","TP2","TP2","TP2","TP2","TP2","TP2","TP2","TP2","TP2","TP2","TP2","TP2","TP2","TP2","TP2","TP3","TP3","TP3","TP3","TP3","TP3","TP3","TP3","TP3","TP3","TP3","TP3","TP3","TP3","TP3","TP3","TP3","TP3","TP3","TP3","TP3","TP3","TP3","TP3","TP3","TP3","TP3","TP3","TP3","TP3","TP3","TP3","TP3","TP3","TP3","TP3"],["bjaccard","bjaccard","bjaccard","wjaccard","wjaccard","wjaccard","rAitchison","rAitchison","rAitchison","bjaccard","bjaccard","bjaccard","wjaccard","wjaccard","wjaccard","rAitchison","rAitchison","rAitchison","bjaccard","bjaccard","bjaccard","wjaccard","wjaccard","wjaccard","rAitchison","rAitchison","rAitchison","bjaccard","bjaccard","bjaccard","wjaccard","wjaccard","wjaccard","rAitchison","rAitchison","rAitchison","bjaccard","bjaccard","bjaccard","wjaccard","wjaccard","wjaccard","rAitchison","rAitchison","rAitchison","bjaccard","bjaccard","bjaccard","wjaccard","wjaccard","wjaccard","rAitchison","rAitchison","rAitchison","bjaccard","bjaccard","bjaccard","wjaccard","wjaccard","wjaccard","rAitchison","rAitchison","rAitchison","bjaccard","bjaccard","bjaccard","wjaccard","wjaccard","wjaccard","rAitchison","rAitchison","rAitchison","bjaccard","bjaccard","bjaccard","wjaccard","wjaccard","wjaccard","rAitchison","rAitchison","rAitchison","bjaccard","bjaccard","bjaccard","wjaccard","wjaccard","wjaccard","rAitchison","rAitchison","rAitchison","bjaccard","bjaccard","bjaccard","wjaccard","wjaccard","wjaccard","rAitchison","rAitchison","rAitchison","bjaccard","bjaccard","bjaccard","wjaccard","wjaccard","wjaccard","rAitchison","rAitchison","rAitchison"],["Sample","Residual","Total","Sample","Residual","Total","Sample","Residual","Total","Sample","Residual","Total","Sample","Residual","Total","Sample","Residual","Total","Sample","Residual","Total","Sample","Residual","Total","Sample","Residual","Total","Sample","Residual","Total","Sample","Residual","Total","Sample","Residual","Total","Sample","Residual","Total","Sample","Residual","Total","Sample","Residual","Total","Sample","Residual","Total","Sample","Residual","Total","Sample","Residual","Total","Sample","Residual","Total","Sample","Residual","Total","Sample","Residual","Total","Sample","Residual","Total","Sample","Residual","Total","Sample","Residual","Total","Sample","Residual","Total","Sample","Residual","Total","Sample","Residual","Total","Sample","Residual","Total","Sample","Residual","Total","Sample","Residual","Total","Sample","Residual","Total","Sample","Residual","Total","Sample","Residual","Total","Sample","Residual","Total","Sample","Residual","Total","Sample","Residual","Total"],[1,78,79,1,78,79,1,78,79,1,78,79,1,78,79,1,78,79,1,78,79,1,78,79,1,78,79,1,78,79,1,78,79,1,78,79,1,78,79,1,78,79,1,78,79,1,78,79,1,78,79,1,78,79,1,78,79,1,78,79,1,78,79,1,78,79,1,78,79,1,78,79,1,77,78,1,77,78,1,77,78,1,77,78,1,77,78,1,77,78,1,77,78,1,77,78,1,77,78,1,77,78,1,77,78,1,77,78],[3.248473994802223,12.35562375156915,15.60409774637138,6.85446847203233,22.65882067574316,29.5132891477755,1.318247494652351,4.429811132201448,5.748058626853797,3.248473994802233,12.35562375156915,15.60409774637138,6.854468472032324,22.65882067574316,29.5132891477755,1.318247494652351,4.429811132201448,5.748058626853797,3.248473994802223,12.35562375156915,15.60409774637138,6.85446847203233,22.65882067574316,29.5132891477755,1.318247494652351,4.429811132201448,5.748058626853797,3.248473994802233,12.35562375156915,15.60409774637138,6.854468472032324,22.65882067574316,29.5132891477755,1.318247494652351,4.429811132201448,5.748058626853797,2.980978549295354,11.8476691708499,14.82864772014526,5.276211156400361,24.92519924920596,30.20141040560633,0.9444654786299971,4.753420341439515,5.697885820069514,2.980978549295353,11.8476691708499,14.82864772014526,5.276211156400367,24.92519924920596,30.20141040560633,0.9444654786299993,4.753420341439515,5.697885820069514,2.980978549295354,11.8476691708499,14.82864772014526,5.276211156400361,24.92519924920596,30.20141040560633,0.9444654786299971,4.753420341439515,5.697885820069514,2.980978549295353,11.8476691708499,14.82864772014526,5.276211156400367,24.92519924920596,30.20141040560633,0.9444654786299993,4.753420341439515,5.697885820069514,3.193385928827109,12.57697341542997,15.77035934425709,6.441306935044484,22.74657072629838,29.18787766134289,1.109826122086349,4.164200420848565,5.274026542934912,3.19338592882711,12.57697341542997,15.77035934425709,6.441306935044498,22.74657072629838,29.18787766134289,1.10982612208635,4.164200420848565,5.274026542934912,3.193385928827109,12.57697341542997,15.77035934425709,6.441306935044484,22.74657072629838,29.18787766134289,1.109826122086349,4.164200420848565,5.274026542934912,3.19338592882711,12.57697341542997,15.77035934425709,6.441306935044498,22.74657072629838,29.18787766134289,1.10982612208635,4.164200420848565,5.274026542934912],[0.2081808283697552,0.7918191716302444,1,0.2322502394670901,0.7677497605329098,1,0.2293378652217914,0.770662134778209,1,0.2081808283697558,0.7918191716302444,1,0.2322502394670899,0.7677497605329098,1,0.2293378652217914,0.770662134778209,1,0.2081808283697552,0.7918191716302444,1,0.2322502394670901,0.7677497605329098,1,0.2293378652217914,0.770662134778209,1,0.2081808283697558,0.7918191716302444,1,0.2322502394670899,0.7677497605329098,1,0.2293378652217914,0.770662134778209,1,0.2010283476655519,0.7989716523344478,1,0.1747008197809507,0.8252991802190489,1,0.1657571787948665,0.8342428212051332,1,0.2010283476655518,0.7989716523344478,1,0.1747008197809509,0.8252991802190489,1,0.1657571787948669,0.8342428212051332,1,0.2010283476655519,0.7989716523344478,1,0.1747008197809507,0.8252991802190489,1,0.1657571787948665,0.8342428212051332,1,0.2010283476655518,0.7989716523344478,1,0.1747008197809509,0.8252991802190489,1,0.1657571787948669,0.8342428212051332,1,0.2024929083172735,0.7975070916827258,1,0.2206843200379554,0.7793156799620438,1,0.2104324111855433,0.789567588814457,1,0.2024929083172736,0.7975070916827258,1,0.2206843200379559,0.7793156799620438,1,0.2104324111855436,0.789567588814457,1,0.2024929083172735,0.7975070916827258,1,0.2206843200379554,0.7793156799620438,1,0.2104324111855433,0.789567588814457,1,0.2024929083172736,0.7975070916827258,1,0.2206843200379559,0.7793156799620438,1,0.2104324111855436,0.789567588814457,1],[20.50733954749913,null,null,23.5956031635343,null,null,23.21166783735634,null,null,20.50733954749919,null,null,23.59560316353427,null,null,23.21166783735634,null,null,20.50733954749913,null,null,23.5956031635343,null,null,23.21166783735634,null,null,20.50733954749919,null,null,23.59560316353427,null,null,23.21166783735634,null,null,19.62549118244479,null,null,16.51118075665288,null,null,15.49795768973173,null,null,19.62549118244478,null,null,16.5111807566529,null,null,15.49795768973177,null,null,19.62549118244479,null,null,16.51118075665288,null,null,15.49795768973173,null,null,19.62549118244478,null,null,16.5111807566529,null,null,15.49795768973177,null,null,19.55086556977357,null,null,21.80463332105699,null,null,20.52173352963517,null,null,19.55086556977358,null,null,21.80463332105704,null,null,20.5217335296352,null,null,19.55086556977357,null,null,21.80463332105699,null,null,20.52173352963517,null,null,19.55086556977358,null,null,21.80463332105704,null,null,20.5217335296352,null,null],[0.001,null,null,0.001,null,null,0.001,null,null,0.001,null,null,0.001,null,null,0.001,null,null,0.001,null,null,0.001,null,null,0.001,null,null,0.001,null,null,0.001,null,null,0.001,null,null,0.001,null,null,0.001,null,null,0.001,null,null,0.001,null,null,0.001,null,null,0.001,null,null,0.001,null,null,0.001,null,null,0.001,null,null,0.001,null,null,0.001,null,null,0.001,null,null,0.001,null,null,0.001,null,null,0.001,null,null,0.001,null,null,0.001,null,null,0.001,null,null,0.001,null,null,0.001,null,null,0.001,null,null,0.001,null,null,0.001,null,null,0.001,null,null],["Sample","Sample","Sample","Sample","Sample","Sample","Sample","Sample","Sample","Sample","Sample","Sample","Sample","Sample","Sample","Sample","Sample","Sample","Sample","Sample","Sample","Sample","Sample","Sample","Sample","Sample","Sample","Sample","Sample","Sample","Sample","Sample","Sample","Sample","Sample","Sample","Sample","Sample","Sample","Sample","Sample","Sample","Sample","Sample","Sample","Sample","Sample","Sample","Sample","Sample","Sample","Sample","Sample","Sample","Sample","Sample","Sample","Sample","Sample","Sample","Sample","Sample","Sample","Sample","Sample","Sample","Sample","Sample","Sample","Sample","Sample","Sample","Sample","Sample","Sample","Sample","Sample","Sample","Sample","Sample","Sample","Sample","Sample","Sample","Sample","Sample","Sample","Sample","Sample","Sample","Sample","Sample","Sample","Sample","Sample","Sample","Sample","Sample","Sample","Sample","Sample","Sample","Sample","Sample","Sample","Sample","Sample","Sample"],["margin","margin","margin","margin","margin","margin","margin","margin","margin","terms","terms","terms","terms","terms","terms","terms","terms","terms","margin","margin","margin","margin","margin","margin","margin","margin","margin","terms","terms","terms","terms","terms","terms","terms","terms","terms","margin","margin","margin","margin","margin","margin","margin","margin","margin","terms","terms","terms","terms","terms","terms","terms","terms","terms","margin","margin","margin","margin","margin","margin","margin","margin","margin","terms","terms","terms","terms","terms","terms","terms","terms","terms","margin","margin","margin","margin","margin","margin","margin","margin","margin","terms","terms","terms","terms","terms","terms","terms","terms","terms","margin","margin","margin","margin","margin","margin","margin","margin","margin","terms","terms","terms","terms","terms","terms","terms","terms","terms"]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>Time<\/th>\n      <th>Distance<\/th>\n      <th>terms<\/th>\n      <th>Df<\/th>\n      <th>SumOfSqs<\/th>\n      <th>R2<\/th>\n      <th>F<\/th>\n      <th>Pr(&gt;F)<\/th>\n      <th>Group<\/th>\n      <th>terms_margin<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":[4,5,6,7,8]},{"orderable":false,"targets":0},{"name":" ","targets":0},{"name":"Time","targets":1},{"name":"Distance","targets":2},{"name":"terms","targets":3},{"name":"Df","targets":4},{"name":"SumOfSqs","targets":5},{"name":"R2","targets":6},{"name":"F","targets":7},{"name":"Pr(>F)","targets":8},{"name":"Group","targets":9},{"name":"terms_margin","targets":10}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
```


At baseline, TP2, TP3: impact of clusterdtp2 in the diff site?



``` r
tmp = ps_up
group = "Time"
Time_out_clust=NULL
# Time_out <- vector("list", length(tmp %>%
#                                get_variable(group) %>%
#                                levels()))
# names(Time_out) <- tmp %>%
#   get_variable(group) %>% levels()

for(tp in tmp %>%
    get_variable(group) %>%
    unique()){
  
  prune_samples(get_variable(tmp, group) == tp,
                tmp) %>%
    subset_taxa(Class != "UNCLASSIFIED") %>% 
    transform_sample_counts(function(x) x/sum(x) * 1) %>% 
    microViz::ps_mutate(cluster_Dtp2 = as.factor(cluster_Dtp2)) %>% 
    compute_plot_beta(ps_up = .,
                      beta = beta,
                      color_group = "Sample",
                      shape_group = "Sample",
                      alpha = NULL,
                      col_pal = sample_pal,
                      fill_pal = sample_pal,
                      path_group = "interaction(Subject)",
                      facet_formula = "Sample ~ .",
                      axis1 = 1,
                      axis2 = 2,
                      seed = 123,
                      permanova_terms = c("Sample", "cluster_Dtp2"),
                      metadata_dist_boxplot = NULL ,
                      strata = "none") -> Time_out_clust[[tp]]
}
```


``` r
Time_out_clust %>% 
  plyr::ldply(., function(x) x$perm) %>%
  rename(Time = '.id') %>% 
  DT::datatable()
```

```{=html}
<div class="datatables html-widget html-fill-item" id="htmlwidget-7531f021030c44f57ef2" style="width:100%;height:auto;"></div>
<script type="application/json" data-for="htmlwidget-7531f021030c44f57ef2">{"x":{"filter":"none","vertical":false,"data":[["1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24","25","26","27","28","29","30","31","32","33","34","35","36","37","38","39","40","41","42","43","44","45","46","47","48","49","50","51","52","53","54","55","56","57","58","59","60","61","62","63","64","65","66","67","68","69","70","71","72","73","74","75","76","77","78","79","80","81","82","83","84","85","86","87","88","89","90","91","92","93","94","95","96","97","98","99","100","101","102","103","104","105","106","107","108","109","110","111","112","113","114","115","116","117","118","119","120","121","122","123","124","125","126","127","128","129","130","131","132","133","134","135","136","137","138","139","140","141","142","143","144"],["TP1","TP1","TP1","TP1","TP1","TP1","TP1","TP1","TP1","TP1","TP1","TP1","TP1","TP1","TP1","TP1","TP1","TP1","TP1","TP1","TP1","TP1","TP1","TP1","TP1","TP1","TP1","TP1","TP1","TP1","TP1","TP1","TP1","TP1","TP1","TP1","TP1","TP1","TP1","TP1","TP1","TP1","TP1","TP1","TP1","TP1","TP1","TP1","TP2","TP2","TP2","TP2","TP2","TP2","TP2","TP2","TP2","TP2","TP2","TP2","TP2","TP2","TP2","TP2","TP2","TP2","TP2","TP2","TP2","TP2","TP2","TP2","TP2","TP2","TP2","TP2","TP2","TP2","TP2","TP2","TP2","TP2","TP2","TP2","TP2","TP2","TP2","TP2","TP2","TP2","TP2","TP2","TP2","TP2","TP2","TP2","TP3","TP3","TP3","TP3","TP3","TP3","TP3","TP3","TP3","TP3","TP3","TP3","TP3","TP3","TP3","TP3","TP3","TP3","TP3","TP3","TP3","TP3","TP3","TP3","TP3","TP3","TP3","TP3","TP3","TP3","TP3","TP3","TP3","TP3","TP3","TP3","TP3","TP3","TP3","TP3","TP3","TP3","TP3","TP3","TP3","TP3","TP3","TP3"],["bjaccard","bjaccard","bjaccard","wjaccard","wjaccard","wjaccard","rAitchison","rAitchison","rAitchison","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","wjaccard","wjaccard","wjaccard","wjaccard","wjaccard","rAitchison","rAitchison","rAitchison","rAitchison","rAitchison","bjaccard","bjaccard","bjaccard","bjaccard","wjaccard","wjaccard","wjaccard","wjaccard","rAitchison","rAitchison","rAitchison","rAitchison","bjaccard","bjaccard","bjaccard","bjaccard","wjaccard","wjaccard","wjaccard","wjaccard","rAitchison","rAitchison","rAitchison","rAitchison","bjaccard","bjaccard","bjaccard","wjaccard","wjaccard","wjaccard","rAitchison","rAitchison","rAitchison","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","wjaccard","wjaccard","wjaccard","wjaccard","wjaccard","rAitchison","rAitchison","rAitchison","rAitchison","rAitchison","bjaccard","bjaccard","bjaccard","bjaccard","wjaccard","wjaccard","wjaccard","wjaccard","rAitchison","rAitchison","rAitchison","rAitchison","bjaccard","bjaccard","bjaccard","bjaccard","wjaccard","wjaccard","wjaccard","wjaccard","rAitchison","rAitchison","rAitchison","rAitchison","bjaccard","bjaccard","bjaccard","wjaccard","wjaccard","wjaccard","rAitchison","rAitchison","rAitchison","bjaccard","bjaccard","bjaccard","bjaccard","bjaccard","wjaccard","wjaccard","wjaccard","wjaccard","wjaccard","rAitchison","rAitchison","rAitchison","rAitchison","rAitchison","bjaccard","bjaccard","bjaccard","bjaccard","wjaccard","wjaccard","wjaccard","wjaccard","rAitchison","rAitchison","rAitchison","rAitchison","bjaccard","bjaccard","bjaccard","bjaccard","wjaccard","wjaccard","wjaccard","wjaccard","rAitchison","rAitchison","rAitchison","rAitchison"],["Sample:cluster_Dtp2","Residual","Total","Sample:cluster_Dtp2","Residual","Total","Sample:cluster_Dtp2","Residual","Total","Sample","cluster_Dtp2","Sample:cluster_Dtp2","Residual","Total","Sample","cluster_Dtp2","Sample:cluster_Dtp2","Residual","Total","Sample","cluster_Dtp2","Sample:cluster_Dtp2","Residual","Total","Sample","cluster_Dtp2","Residual","Total","Sample","cluster_Dtp2","Residual","Total","Sample","cluster_Dtp2","Residual","Total","Sample","cluster_Dtp2","Residual","Total","Sample","cluster_Dtp2","Residual","Total","Sample","cluster_Dtp2","Residual","Total","Sample:cluster_Dtp2","Residual","Total","Sample:cluster_Dtp2","Residual","Total","Sample:cluster_Dtp2","Residual","Total","Sample","cluster_Dtp2","Sample:cluster_Dtp2","Residual","Total","Sample","cluster_Dtp2","Sample:cluster_Dtp2","Residual","Total","Sample","cluster_Dtp2","Sample:cluster_Dtp2","Residual","Total","Sample","cluster_Dtp2","Residual","Total","Sample","cluster_Dtp2","Residual","Total","Sample","cluster_Dtp2","Residual","Total","Sample","cluster_Dtp2","Residual","Total","Sample","cluster_Dtp2","Residual","Total","Sample","cluster_Dtp2","Residual","Total","Sample:cluster_Dtp2","Residual","Total","Sample:cluster_Dtp2","Residual","Total","Sample:cluster_Dtp2","Residual","Total","Sample","cluster_Dtp2","Sample:cluster_Dtp2","Residual","Total","Sample","cluster_Dtp2","Sample:cluster_Dtp2","Residual","Total","Sample","cluster_Dtp2","Sample:cluster_Dtp2","Residual","Total","Sample","cluster_Dtp2","Residual","Total","Sample","cluster_Dtp2","Residual","Total","Sample","cluster_Dtp2","Residual","Total","Sample","cluster_Dtp2","Residual","Total","Sample","cluster_Dtp2","Residual","Total","Sample","cluster_Dtp2","Residual","Total"],[2,74,79,2,74,79,2,74,79,1,2,2,74,79,1,2,2,74,79,1,2,2,74,79,1,2,76,79,1,2,76,79,1,2,76,79,1,2,76,79,1,2,76,79,1,2,76,79,2,74,79,2,74,79,2,74,79,1,2,2,74,79,1,2,2,74,79,1,2,2,74,79,1,2,76,79,1,2,76,79,1,2,76,79,1,2,76,79,1,2,76,79,1,2,76,79,2,73,78,2,73,78,2,73,78,1,2,2,73,78,1,2,2,73,78,1,2,2,73,78,1,2,75,78,1,2,75,78,1,2,75,78,1,2,75,78,1,2,75,78,1,2,75,78],[0.1899081173955963,11.75059825694362,15.60409774637138,0.4758552371934641,21.64301412970805,29.5132891477755,0.1023239468828825,4.212869684694672,5.748058626853797,3.248473994802233,0.4151173772299228,0.1899081173955928,11.75059825694362,15.60409774637138,6.854468472032324,0.5399513088416503,0.4758552371934677,21.64301412970805,29.5132891477755,1.318247494652351,0.1146175006238903,0.102323946882882,4.212869684694672,5.748058626853797,3.248473994802232,0.4151173772299259,11.94050637433922,15.60409774637138,6.854468472032327,0.5399513088416512,22.11886936690152,29.5132891477755,1.318247494652346,0.1146175006238925,4.315193631577555,5.748058626853797,3.248473994802233,0.4151173772299228,11.94050637433922,15.60409774637138,6.854468472032324,0.5399513088416503,22.11886936690152,29.5132891477755,1.318247494652351,0.1146175006238903,4.315193631577555,5.748058626853797,0.1628359892694871,11.16951849610493,14.82864772014526,0.6275550569200803,23.4388368243617,30.20141040560633,0.1032067139931314,4.492109873898399,5.697885820069514,2.980978549295353,0.5153146854754787,0.1628359892694911,11.16951849610493,14.82864772014526,5.276211156400367,0.8588073679241806,0.6275550569200758,23.4388368243617,30.20141040560633,0.9444654786299993,0.1581037535479826,0.1032067139931325,4.492109873898399,5.697885820069514,2.980978549295333,0.5153146854754791,11.33235448537442,14.82864772014526,5.27621115640034,0.858807367924193,24.06639188128177,30.20141040560633,0.9444654786299962,0.1581037535479801,4.595316587891532,5.697885820069514,2.980978549295353,0.5153146854754787,11.33235448537442,14.82864772014526,5.276211156400367,0.8588073679241806,24.06639188128177,30.20141040560633,0.9444654786299993,0.1581037535479826,4.595316587891532,5.697885820069514,0.2214590880396443,11.97524625920603,15.77035934425709,0.5188249453311577,21.6846789874654,29.18787766134289,0.1013507555069912,3.947253763455476,5.274026542934912,3.19338592882711,0.3802680681842956,0.2214590880396514,11.97524625920603,15.77035934425709,6.441306935044498,0.5430667935018239,0.5188249453311791,21.6846789874654,29.18787766134289,1.10982612208635,0.1155959018860966,0.1013507555069926,3.947253763455476,5.274026542934912,3.1961036509494,0.3802680681842929,12.19670534724568,15.77035934425709,6.437153807255001,0.5430667935017972,22.20350393279658,29.18787766134289,1.108519006774032,0.1155959018860964,4.048604518962468,5.274026542934912,3.19338592882711,0.3802680681842956,12.19670534724568,15.77035934425709,6.441306935044498,0.5430667935018239,22.20350393279658,29.18787766134289,1.10982612208635,0.1155959018860966,4.048604518962468,5.274026542934912],[0.0121704003962522,0.7530456709473089,1,0.01612342273376638,0.7333311452118967,1,0.01780147933857967,0.7329204446546483,1,0.2081808283697558,0.02660310028668306,0.01217040039625198,0.7530456709473089,1,0.2322502394670899,0.01829519258724668,0.0161234227337665,0.7333311452118967,1,0.2293378652217914,0.0199402107849805,0.01780147933857959,0.7329204446546483,1,0.2081808283697558,0.02660310028668325,0.7652160713435611,1,0.2322502394670899,0.01829519258724671,0.7494545679456633,1,0.2293378652217905,0.01994021078498089,0.7507219239932281,1,0.2081808283697558,0.02660310028668306,0.7652160713435611,1,0.2322502394670899,0.01829519258724668,0.7494545679456633,1,0.2293378652217914,0.0199402107849805,0.7507219239932281,1,0.01098117591992346,0.7532391831610331,1,0.02077899834782505,0.776084179830579,1,0.01811315938090741,0.7883818693024627,1,0.2010283476655518,0.03475129325349101,0.01098117591992373,0.7532391831610331,1,0.1747008197809509,0.02843600204064506,0.02077899834782491,0.776084179830579,1,0.1657571787948669,0.0277477925217627,0.0181131593809076,0.7883818693024627,1,0.2010283476655505,0.03475129325349104,0.7642203590809568,1,0.17470081978095,0.02843600204064547,0.7968631781784037,1,0.1657571787948664,0.02774779252176227,0.8064950286833704,1,0.2010283476655518,0.03475129325349101,0.7642203590809568,1,0.1747008197809509,0.02843600204064506,0.7968631781784037,1,0.1657571787948669,0.0277477925217627,0.8064950286833704,1,0.01404274203303367,0.7593515149397605,1,0.017775356994123,0.7429344208943666,1,0.01921695969519925,0.7484326692938659,1,0.2024929083172736,0.02411283470993154,0.01404274203303412,0.7593515149397605,1,0.2206843200379559,0.01860590207355413,0.01777535699412373,0.7429344208943666,1,0.2104324111855436,0.02191795982539164,0.0192169596951995,0.7484326692938659,1,0.2026652393379538,0.02411283470993137,0.7733942569727945,1,0.2205420305629319,0.01860590207355322,0.7607097778884903,1,0.210184571076724,0.0219179598253916,0.7676496289890652,1,0.2024929083172736,0.02411283470993154,0.7733942569727945,1,0.2206843200379559,0.01860590207355413,0.7607097778884903,1,0.2104324111855436,0.02191795982539164,0.7676496289890652,1],[0.5979780935396143,null,null,0.8135023925336999,null,null,0.8986715274913711,null,null,20.45743292034655,1.307111571824103,0.5979780935396032,null,null,23.43623045711306,0.9230783802759803,0.813502392533706,null,null,23.15531262661023,1.006641040545572,0.8986715274913673,null,null,20.67617702843294,1.321088054409262,null,null,23.5518188218963,0.9276310373569963,null,null,23.21722224941082,1.009332464674508,null,null,20.67617702843294,1.321088054409252,null,null,23.5518188218963,0.9276310373569947,null,null,23.21722224941092,1.009332464674489,null,null,0.5394083554337687,null,null,0.9906437456789325,null,null,0.8500790330027959,null,null,19.74950063646717,1.707024646518263,0.539408355433782,null,null,16.65780723247387,1.355693239016353,0.9906437456789255,null,null,15.55848974770663,1.302247506292333,0.850079033002805,null,null,19.99181812030653,1.727969070623475,null,null,16.6619096815384,1.356027116241789,null,null,15.62011561184171,1.307405598703237,null,null,19.99181812030667,1.727969070623474,null,null,16.66190968153849,1.35602711624177,null,null,15.62011561184176,1.307405598703257,null,null,0.6749971180954191,null,null,0.8732944820411526,null,null,0.9371838745849378,null,null,19.46658697103361,1.159039587854541,0.6749971180954407,null,null,21.68422260389704,0.9140987502869854,0.8732944820411883,null,null,20.52498059850604,1.068907820902029,0.9371838745849501,null,null,19.65348567474715,1.169172506092496,null,null,21.74370932648184,0.91719779085122,null,null,20.53520543156395,1.070701349174876,null,null,19.63677385353244,1.169172506092505,null,null,21.75773794940347,0.9171977908512651,null,null,20.5594196139976,1.070701349174878,null,null],[0.989,null,null,0.787,null,null,0.581,null,null,0.001,0.118,0.994,null,null,0.001,0.501,0.773,null,null,0.001,0.352,0.5669999999999999,null,null,0.001,0.099,null,null,0.001,0.544,null,null,0.001,0.409,null,null,0.001,0.112,null,null,0.001,0.517,null,null,0.001,0.37,null,null,1,null,null,0.469,null,null,0.805,null,null,0.001,0.027,0.996,null,null,0.001,0.059,0.442,null,null,0.001,0.065,0.791,null,null,0.001,0.015,null,null,0.001,0.063,null,null,0.001,0.082,null,null,0.001,0.019,null,null,0.001,0.064,null,null,0.001,0.065,null,null,0.967,null,null,0.671,null,null,0.517,null,null,0.001,0.228,0.98,null,null,0.001,0.579,0.678,null,null,0.001,0.288,0.524,null,null,0.001,0.204,null,null,0.001,0.546,null,null,0.001,0.294,null,null,0.001,0.201,null,null,0.001,0.555,null,null,0.001,0.269,null,null],["Sample * cluster_Dtp2","Sample * cluster_Dtp2","Sample * cluster_Dtp2","Sample * cluster_Dtp2","Sample * cluster_Dtp2","Sample * cluster_Dtp2","Sample * cluster_Dtp2","Sample * cluster_Dtp2","Sample * cluster_Dtp2","Sample * cluster_Dtp2","Sample * cluster_Dtp2","Sample * cluster_Dtp2","Sample * cluster_Dtp2","Sample * cluster_Dtp2","Sample * cluster_Dtp2","Sample * cluster_Dtp2","Sample * cluster_Dtp2","Sample * cluster_Dtp2","Sample * cluster_Dtp2","Sample * cluster_Dtp2","Sample * cluster_Dtp2","Sample * cluster_Dtp2","Sample * cluster_Dtp2","Sample * cluster_Dtp2","Sample + cluster_Dtp2","Sample + cluster_Dtp2","Sample + cluster_Dtp2","Sample + cluster_Dtp2","Sample + cluster_Dtp2","Sample + cluster_Dtp2","Sample + cluster_Dtp2","Sample + cluster_Dtp2","Sample + cluster_Dtp2","Sample + cluster_Dtp2","Sample + cluster_Dtp2","Sample + cluster_Dtp2","Sample + cluster_Dtp2","Sample + cluster_Dtp2","Sample + cluster_Dtp2","Sample + cluster_Dtp2","Sample + cluster_Dtp2","Sample + cluster_Dtp2","Sample + cluster_Dtp2","Sample + cluster_Dtp2","Sample + cluster_Dtp2","Sample + cluster_Dtp2","Sample + cluster_Dtp2","Sample + cluster_Dtp2","Sample * cluster_Dtp2","Sample * cluster_Dtp2","Sample * cluster_Dtp2","Sample * cluster_Dtp2","Sample * cluster_Dtp2","Sample * cluster_Dtp2","Sample * cluster_Dtp2","Sample * cluster_Dtp2","Sample * cluster_Dtp2","Sample * cluster_Dtp2","Sample * cluster_Dtp2","Sample * cluster_Dtp2","Sample * cluster_Dtp2","Sample * cluster_Dtp2","Sample * cluster_Dtp2","Sample * cluster_Dtp2","Sample * cluster_Dtp2","Sample * cluster_Dtp2","Sample * cluster_Dtp2","Sample * cluster_Dtp2","Sample * cluster_Dtp2","Sample * cluster_Dtp2","Sample * cluster_Dtp2","Sample * cluster_Dtp2","Sample + cluster_Dtp2","Sample + cluster_Dtp2","Sample + cluster_Dtp2","Sample + cluster_Dtp2","Sample + cluster_Dtp2","Sample + cluster_Dtp2","Sample + cluster_Dtp2","Sample + cluster_Dtp2","Sample + cluster_Dtp2","Sample + cluster_Dtp2","Sample + cluster_Dtp2","Sample + cluster_Dtp2","Sample + cluster_Dtp2","Sample + cluster_Dtp2","Sample + cluster_Dtp2","Sample + cluster_Dtp2","Sample + cluster_Dtp2","Sample + cluster_Dtp2","Sample + cluster_Dtp2","Sample + cluster_Dtp2","Sample + cluster_Dtp2","Sample + cluster_Dtp2","Sample + cluster_Dtp2","Sample + cluster_Dtp2","Sample * cluster_Dtp2","Sample * cluster_Dtp2","Sample * cluster_Dtp2","Sample * cluster_Dtp2","Sample * cluster_Dtp2","Sample * cluster_Dtp2","Sample * cluster_Dtp2","Sample * cluster_Dtp2","Sample * cluster_Dtp2","Sample * cluster_Dtp2","Sample * cluster_Dtp2","Sample * cluster_Dtp2","Sample * cluster_Dtp2","Sample * cluster_Dtp2","Sample * cluster_Dtp2","Sample * cluster_Dtp2","Sample * cluster_Dtp2","Sample * cluster_Dtp2","Sample * cluster_Dtp2","Sample * cluster_Dtp2","Sample * cluster_Dtp2","Sample * cluster_Dtp2","Sample * cluster_Dtp2","Sample * cluster_Dtp2","Sample + cluster_Dtp2","Sample + cluster_Dtp2","Sample + cluster_Dtp2","Sample + cluster_Dtp2","Sample + cluster_Dtp2","Sample + cluster_Dtp2","Sample + cluster_Dtp2","Sample + cluster_Dtp2","Sample + cluster_Dtp2","Sample + cluster_Dtp2","Sample + cluster_Dtp2","Sample + cluster_Dtp2","Sample + cluster_Dtp2","Sample + cluster_Dtp2","Sample + cluster_Dtp2","Sample + cluster_Dtp2","Sample + cluster_Dtp2","Sample + cluster_Dtp2","Sample + cluster_Dtp2","Sample + cluster_Dtp2","Sample + cluster_Dtp2","Sample + cluster_Dtp2","Sample + cluster_Dtp2","Sample + cluster_Dtp2"],["margin","margin","margin","margin","margin","margin","margin","margin","margin","terms","terms","terms","terms","terms","terms","terms","terms","terms","terms","terms","terms","terms","terms","terms","margin","margin","margin","margin","margin","margin","margin","margin","margin","margin","margin","margin","terms","terms","terms","terms","terms","terms","terms","terms","terms","terms","terms","terms","margin","margin","margin","margin","margin","margin","margin","margin","margin","terms","terms","terms","terms","terms","terms","terms","terms","terms","terms","terms","terms","terms","terms","terms","margin","margin","margin","margin","margin","margin","margin","margin","margin","margin","margin","margin","terms","terms","terms","terms","terms","terms","terms","terms","terms","terms","terms","terms","margin","margin","margin","margin","margin","margin","margin","margin","margin","terms","terms","terms","terms","terms","terms","terms","terms","terms","terms","terms","terms","terms","terms","terms","margin","margin","margin","margin","margin","margin","margin","margin","margin","margin","margin","margin","terms","terms","terms","terms","terms","terms","terms","terms","terms","terms","terms","terms"]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>Time<\/th>\n      <th>Distance<\/th>\n      <th>terms<\/th>\n      <th>Df<\/th>\n      <th>SumOfSqs<\/th>\n      <th>R2<\/th>\n      <th>F<\/th>\n      <th>Pr(&gt;F)<\/th>\n      <th>Group<\/th>\n      <th>terms_margin<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":[4,5,6,7,8]},{"orderable":false,"targets":0},{"name":" ","targets":0},{"name":"Time","targets":1},{"name":"Distance","targets":2},{"name":"terms","targets":3},{"name":"Df","targets":4},{"name":"SumOfSqs","targets":5},{"name":"R2","targets":6},{"name":"F","targets":7},{"name":"Pr(>F)","targets":8},{"name":"Group","targets":9},{"name":"terms_margin","targets":10}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
```

Saliva


``` r
tmp = ps_up %>% 
  subset_samples(Sample == "Saliva")
group = "Time"
Time_out_saliva_clust=NULL
# Time_out <- vector("list", length(tmp %>%
#                                get_variable(group) %>%
#                                levels()))
# names(Time_out) <- tmp %>%
#   get_variable(group) %>% levels()

for(tp in tmp %>%
    get_variable(group) %>%
    unique()){
  
  prune_samples(get_variable(tmp, group) == tp,
                tmp) %>%
    subset_taxa(Class != "UNCLASSIFIED") %>% 
    transform_sample_counts(function(x) x/sum(x) * 1) %>% 
    microViz::ps_mutate(cluster_Dtp2 = as.factor(cluster_Dtp2)) %>% 
    compute_plot_beta(ps_up = .,
                      beta = beta,
                      color_group = "Sample",
                      shape_group = "Sample",
                      alpha = NULL,
                      col_pal = sample_pal,
                      fill_pal = sample_pal,
                      path_group = "interaction(Subject)",
                      facet_formula = "Sample ~ .",
                      axis1 = 1,
                      axis2 = 2,
                      seed = 123,
                      permanova_terms = c("cluster_Dtp2"),
                      metadata_dist_boxplot = NULL ,
                      strata = "none") -> Time_out_saliva_clust[[tp]]
}
```


``` r
Time_out_saliva_clust %>% 
  plyr::ldply(., function(x) x$perm) %>%
  rename(Time = '.id') %>% 
  DT::datatable()
```

```{=html}
<div class="datatables html-widget html-fill-item" id="htmlwidget-983315589fed18d622a7" style="width:100%;height:auto;"></div>
<script type="application/json" data-for="htmlwidget-983315589fed18d622a7">{"x":{"filter":"none","vertical":false,"data":[["1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24","25","26","27","28","29","30","31","32","33","34","35","36","37","38","39","40","41","42","43","44","45","46","47","48","49","50","51","52","53","54","55","56","57","58","59","60","61","62","63","64","65","66","67","68","69","70","71","72","73","74","75","76","77","78","79","80","81","82","83","84","85","86","87","88","89","90","91","92","93","94","95","96","97","98","99","100","101","102","103","104","105","106","107","108"],["TP1","TP1","TP1","TP1","TP1","TP1","TP1","TP1","TP1","TP1","TP1","TP1","TP1","TP1","TP1","TP1","TP1","TP1","TP1","TP1","TP1","TP1","TP1","TP1","TP1","TP1","TP1","TP1","TP1","TP1","TP1","TP1","TP1","TP1","TP1","TP1","TP2","TP2","TP2","TP2","TP2","TP2","TP2","TP2","TP2","TP2","TP2","TP2","TP2","TP2","TP2","TP2","TP2","TP2","TP2","TP2","TP2","TP2","TP2","TP2","TP2","TP2","TP2","TP2","TP2","TP2","TP2","TP2","TP2","TP2","TP2","TP2","TP3","TP3","TP3","TP3","TP3","TP3","TP3","TP3","TP3","TP3","TP3","TP3","TP3","TP3","TP3","TP3","TP3","TP3","TP3","TP3","TP3","TP3","TP3","TP3","TP3","TP3","TP3","TP3","TP3","TP3","TP3","TP3","TP3","TP3","TP3","TP3"],["bjaccard","bjaccard","bjaccard","wjaccard","wjaccard","wjaccard","rAitchison","rAitchison","rAitchison","bjaccard","bjaccard","bjaccard","wjaccard","wjaccard","wjaccard","rAitchison","rAitchison","rAitchison","bjaccard","bjaccard","bjaccard","wjaccard","wjaccard","wjaccard","rAitchison","rAitchison","rAitchison","bjaccard","bjaccard","bjaccard","wjaccard","wjaccard","wjaccard","rAitchison","rAitchison","rAitchison","bjaccard","bjaccard","bjaccard","wjaccard","wjaccard","wjaccard","rAitchison","rAitchison","rAitchison","bjaccard","bjaccard","bjaccard","wjaccard","wjaccard","wjaccard","rAitchison","rAitchison","rAitchison","bjaccard","bjaccard","bjaccard","wjaccard","wjaccard","wjaccard","rAitchison","rAitchison","rAitchison","bjaccard","bjaccard","bjaccard","wjaccard","wjaccard","wjaccard","rAitchison","rAitchison","rAitchison","bjaccard","bjaccard","bjaccard","wjaccard","wjaccard","wjaccard","rAitchison","rAitchison","rAitchison","bjaccard","bjaccard","bjaccard","wjaccard","wjaccard","wjaccard","rAitchison","rAitchison","rAitchison","bjaccard","bjaccard","bjaccard","wjaccard","wjaccard","wjaccard","rAitchison","rAitchison","rAitchison","bjaccard","bjaccard","bjaccard","wjaccard","wjaccard","wjaccard","rAitchison","rAitchison","rAitchison"],["cluster_Dtp2","Residual","Total","cluster_Dtp2","Residual","Total","cluster_Dtp2","Residual","Total","cluster_Dtp2","Residual","Total","cluster_Dtp2","Residual","Total","cluster_Dtp2","Residual","Total","cluster_Dtp2","Residual","Total","cluster_Dtp2","Residual","Total","cluster_Dtp2","Residual","Total","cluster_Dtp2","Residual","Total","cluster_Dtp2","Residual","Total","cluster_Dtp2","Residual","Total","cluster_Dtp2","Residual","Total","cluster_Dtp2","Residual","Total","cluster_Dtp2","Residual","Total","cluster_Dtp2","Residual","Total","cluster_Dtp2","Residual","Total","cluster_Dtp2","Residual","Total","cluster_Dtp2","Residual","Total","cluster_Dtp2","Residual","Total","cluster_Dtp2","Residual","Total","cluster_Dtp2","Residual","Total","cluster_Dtp2","Residual","Total","cluster_Dtp2","Residual","Total","cluster_Dtp2","Residual","Total","cluster_Dtp2","Residual","Total","cluster_Dtp2","Residual","Total","cluster_Dtp2","Residual","Total","cluster_Dtp2","Residual","Total","cluster_Dtp2","Residual","Total","cluster_Dtp2","Residual","Total","cluster_Dtp2","Residual","Total","cluster_Dtp2","Residual","Total","cluster_Dtp2","Residual","Total","cluster_Dtp2","Residual","Total","cluster_Dtp2","Residual","Total"],[2,37,39,2,37,39,2,37,39,2,37,39,2,37,39,2,37,39,2,37,39,2,37,39,2,37,39,2,37,39,2,37,39,2,37,39,2,37,39,2,37,39,2,37,39,2,37,39,2,37,39,2,37,39,2,37,39,2,37,39,2,37,39,2,37,39,2,37,39,2,37,39,2,37,39,2,37,39,2,37,39,2,37,39,2,37,39,2,37,39,2,37,39,2,37,39,2,37,39,2,37,39,2,37,39,2,37,39],[0.2338107665787152,4.788719897153503,5.02253066373222,0.3945439815163461,8.772918541525685,9.167462523042035,0.09303147350633423,1.919304554073219,2.012336027579555,0.2338107665787189,4.788719897153503,5.02253066373222,0.3945439815163492,8.772918541525685,9.167462523042035,0.09303147350633552,1.919304554073219,2.012336027579555,0.2338107665787152,4.788719897153503,5.02253066373222,0.3945439815163461,8.772918541525685,9.167462523042035,0.09303147350633423,1.919304554073219,2.012336027579555,0.2338107665787189,4.788719897153503,5.02253066373222,0.3945439815163492,8.772918541525685,9.167462523042035,0.09303147350633552,1.919304554073219,2.012336027579555,0.2391203367558168,4.419437377594229,4.658557714350048,0.4828000688934591,9.924799820175593,10.40759988906906,0.1212455515288506,2.18527631368244,2.306521865211293,0.2391203367558209,4.419437377594229,4.658557714350048,0.4828000688934643,9.924799820175593,10.40759988906906,0.1212455515288517,2.18527631368244,2.306521865211293,0.2391203367558168,4.419437377594229,4.658557714350048,0.4828000688934591,9.924799820175593,10.40759988906906,0.1212455515288506,2.18527631368244,2.306521865211293,0.2391203367558209,4.419437377594229,4.658557714350048,0.4828000688934643,9.924799820175593,10.40759988906906,0.1212455515288517,2.18527631368244,2.306521865211293,0.2498259955154341,4.948246985501278,5.198072981016713,0.4388506343102634,8.48609762409437,8.924948258404637,0.09719428882890835,1.83740612572582,1.934600414554729,0.2498259955154368,4.948246985501278,5.198072981016713,0.4388506343102682,8.48609762409437,8.924948258404637,0.09719428882890843,1.83740612572582,1.934600414554729,0.2498259955154341,4.948246985501278,5.198072981016713,0.4388506343102634,8.48609762409437,8.924948258404637,0.09719428882890835,1.83740612572582,1.934600414554729,0.2498259955154368,4.948246985501278,5.198072981016713,0.4388506343102682,8.48609762409437,8.924948258404637,0.09719428882890843,1.83740612572582,1.934600414554729],[0.04655238210232678,0.9534476178976728,1,0.04303742508078721,0.9569625749192124,1,0.04623058586206044,0.9537694141379389,1,0.04655238210232751,0.9534476178976728,1,0.04303742508078755,0.9569625749192124,1,0.04623058586206109,0.9537694141379389,1,0.04655238210232678,0.9534476178976728,1,0.04303742508078721,0.9569625749192124,1,0.04623058586206044,0.9537694141379389,1,0.04655238210232751,0.9534476178976728,1,0.04303742508078755,0.9569625749192124,1,0.04623058586206109,0.9537694141379389,1,0.05132926356568245,0.9486707364343172,1,0.04638918425376214,0.953610815746237,1,0.05256640023993171,0.9474335997600675,1,0.05132926356568332,0.9486707364343172,1,0.04638918425376264,0.953610815746237,1,0.05256640023993217,0.9474335997600675,1,0.05132926356568245,0.9486707364343172,1,0.04638918425376214,0.953610815746237,1,0.05256640023993171,0.9474335997600675,1,0.05132926356568332,0.9486707364343172,1,0.04638918425376264,0.953610815746237,1,0.05256640023993217,0.9474335997600675,1,0.04806127124951784,0.951938728750482,1,0.04917122448267384,0.9508287755173258,1,0.05023998139237387,0.9497600186076259,1,0.04806127124951837,0.951938728750482,1,0.04917122448267438,0.9508287755173258,1,0.05023998139237391,0.9497600186076259,1,0.04806127124951784,0.951938728750482,1,0.04917122448267384,0.9508287755173258,1,0.05023998139237387,0.9497600186076259,1,0.04806127124951837,0.951938728750482,1,0.04917122448267438,0.9508287755173258,1,0.05023998139237391,0.9497600186076259,1],[0.9032683628619379,null,null,0.8319994792500413,null,null,0.8967218132289839,null,null,0.903268362861952,null,null,0.831999479250048,null,null,0.8967218132289962,null,null,0.9032683628619379,null,null,0.8319994792500413,null,null,0.8967218132289839,null,null,0.903268362861952,null,null,0.831999479250048,null,null,0.8967218132289962,null,null,1.00097045212364,null,null,0.8999477507215825,null,null,1.026434364038822,null,null,1.000970452123657,null,null,0.8999477507215924,null,null,1.026434364038831,null,null,1.00097045212364,null,null,0.8999477507215825,null,null,1.026434364038822,null,null,1.000970452123657,null,null,0.8999477507215924,null,null,1.026434364038831,null,null,0.9340238938310239,null,null,0.956710268296766,null,null,0.978604739670449,null,null,0.9340238938310342,null,null,0.9567102682967764,null,null,0.9786047396704498,null,null,0.9340238938310239,null,null,0.956710268296766,null,null,0.978604739670449,null,null,0.9340238938310342,null,null,0.9567102682967764,null,null,0.9786047396704498,null,null],[0.669,null,null,0.787,null,null,0.877,null,null,0.68,null,null,0.803,null,null,0.843,null,null,0.696,null,null,0.803,null,null,0.867,null,null,0.678,null,null,0.825,null,null,0.869,null,null,0.442,null,null,0.709,null,null,0.34,null,null,0.429,null,null,0.706,null,null,0.343,null,null,0.456,null,null,0.719,null,null,0.357,null,null,0.445,null,null,0.731,null,null,0.397,null,null,0.625,null,null,0.521,null,null,0.552,null,null,0.619,null,null,0.534,null,null,0.5629999999999999,null,null,0.631,null,null,0.513,null,null,0.532,null,null,0.643,null,null,0.525,null,null,0.5629999999999999,null,null],["cluster_Dtp2","cluster_Dtp2","cluster_Dtp2","cluster_Dtp2","cluster_Dtp2","cluster_Dtp2","cluster_Dtp2","cluster_Dtp2","cluster_Dtp2","cluster_Dtp2","cluster_Dtp2","cluster_Dtp2","cluster_Dtp2","cluster_Dtp2","cluster_Dtp2","cluster_Dtp2","cluster_Dtp2","cluster_Dtp2","cluster_Dtp2","cluster_Dtp2","cluster_Dtp2","cluster_Dtp2","cluster_Dtp2","cluster_Dtp2","cluster_Dtp2","cluster_Dtp2","cluster_Dtp2","cluster_Dtp2","cluster_Dtp2","cluster_Dtp2","cluster_Dtp2","cluster_Dtp2","cluster_Dtp2","cluster_Dtp2","cluster_Dtp2","cluster_Dtp2","cluster_Dtp2","cluster_Dtp2","cluster_Dtp2","cluster_Dtp2","cluster_Dtp2","cluster_Dtp2","cluster_Dtp2","cluster_Dtp2","cluster_Dtp2","cluster_Dtp2","cluster_Dtp2","cluster_Dtp2","cluster_Dtp2","cluster_Dtp2","cluster_Dtp2","cluster_Dtp2","cluster_Dtp2","cluster_Dtp2","cluster_Dtp2","cluster_Dtp2","cluster_Dtp2","cluster_Dtp2","cluster_Dtp2","cluster_Dtp2","cluster_Dtp2","cluster_Dtp2","cluster_Dtp2","cluster_Dtp2","cluster_Dtp2","cluster_Dtp2","cluster_Dtp2","cluster_Dtp2","cluster_Dtp2","cluster_Dtp2","cluster_Dtp2","cluster_Dtp2","cluster_Dtp2","cluster_Dtp2","cluster_Dtp2","cluster_Dtp2","cluster_Dtp2","cluster_Dtp2","cluster_Dtp2","cluster_Dtp2","cluster_Dtp2","cluster_Dtp2","cluster_Dtp2","cluster_Dtp2","cluster_Dtp2","cluster_Dtp2","cluster_Dtp2","cluster_Dtp2","cluster_Dtp2","cluster_Dtp2","cluster_Dtp2","cluster_Dtp2","cluster_Dtp2","cluster_Dtp2","cluster_Dtp2","cluster_Dtp2","cluster_Dtp2","cluster_Dtp2","cluster_Dtp2","cluster_Dtp2","cluster_Dtp2","cluster_Dtp2","cluster_Dtp2","cluster_Dtp2","cluster_Dtp2","cluster_Dtp2","cluster_Dtp2","cluster_Dtp2"],["margin","margin","margin","margin","margin","margin","margin","margin","margin","terms","terms","terms","terms","terms","terms","terms","terms","terms","margin","margin","margin","margin","margin","margin","margin","margin","margin","terms","terms","terms","terms","terms","terms","terms","terms","terms","margin","margin","margin","margin","margin","margin","margin","margin","margin","terms","terms","terms","terms","terms","terms","terms","terms","terms","margin","margin","margin","margin","margin","margin","margin","margin","margin","terms","terms","terms","terms","terms","terms","terms","terms","terms","margin","margin","margin","margin","margin","margin","margin","margin","margin","terms","terms","terms","terms","terms","terms","terms","terms","terms","margin","margin","margin","margin","margin","margin","margin","margin","margin","terms","terms","terms","terms","terms","terms","terms","terms","terms"]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>Time<\/th>\n      <th>Distance<\/th>\n      <th>terms<\/th>\n      <th>Df<\/th>\n      <th>SumOfSqs<\/th>\n      <th>R2<\/th>\n      <th>F<\/th>\n      <th>Pr(&gt;F)<\/th>\n      <th>Group<\/th>\n      <th>terms_margin<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":[4,5,6,7,8]},{"orderable":false,"targets":0},{"name":" ","targets":0},{"name":"Time","targets":1},{"name":"Distance","targets":2},{"name":"terms","targets":3},{"name":"Df","targets":4},{"name":"SumOfSqs","targets":5},{"name":"R2","targets":6},{"name":"F","targets":7},{"name":"Pr(>F)","targets":8},{"name":"Group","targets":9},{"name":"terms_margin","targets":10}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
```

Plaque


``` r
tmp = ps_up %>% 
  subset_samples(Sample == "Plaque")
group = "Time"
Time_out_plaque_clust=NULL
# Time_out <- vector("list", length(tmp %>%
#                                get_variable(group) %>%
#                                levels()))
# names(Time_out) <- tmp %>%
#   get_variable(group) %>% levels()

for(tp in tmp %>%
    get_variable(group) %>%
    unique()){
  
  prune_samples(get_variable(tmp, group) == tp,
                tmp) %>%
    subset_taxa(Class != "UNCLASSIFIED") %>% 
    transform_sample_counts(function(x) x/sum(x) * 1) %>% 
    microViz::ps_mutate(cluster_Dtp2 = as.factor(cluster_Dtp2)) %>% 
    compute_plot_beta(ps_up = .,
                      beta = beta,
                      color_group = "Sample",
                      shape_group = "Sample",
                      alpha = NULL,
                      col_pal = sample_pal,
                      fill_pal = sample_pal,
                      path_group = "interaction(Subject)",
                      facet_formula = "Sample ~ .",
                      axis1 = 1,
                      axis2 = 2,
                      seed = 123,
                      permanova_terms = c("cluster_Dtp2"),
                      metadata_dist_boxplot = NULL ,
                      strata = "none") -> Time_out_plaque_clust[[tp]]
}
```


``` r
Time_out_plaque_clust %>% 
  plyr::ldply(., function(x) x$perm) %>%
  rename(Time = '.id') %>% 
  DT::datatable()
```

```{=html}
<div class="datatables html-widget html-fill-item" id="htmlwidget-6f1d8b9e786d9d3fe28d" style="width:100%;height:auto;"></div>
<script type="application/json" data-for="htmlwidget-6f1d8b9e786d9d3fe28d">{"x":{"filter":"none","vertical":false,"data":[["1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24","25","26","27","28","29","30","31","32","33","34","35","36","37","38","39","40","41","42","43","44","45","46","47","48","49","50","51","52","53","54","55","56","57","58","59","60","61","62","63","64","65","66","67","68","69","70","71","72","73","74","75","76","77","78","79","80","81","82","83","84","85","86","87","88","89","90","91","92","93","94","95","96","97","98","99","100","101","102","103","104","105","106","107","108"],["TP1","TP1","TP1","TP1","TP1","TP1","TP1","TP1","TP1","TP1","TP1","TP1","TP1","TP1","TP1","TP1","TP1","TP1","TP1","TP1","TP1","TP1","TP1","TP1","TP1","TP1","TP1","TP1","TP1","TP1","TP1","TP1","TP1","TP1","TP1","TP1","TP2","TP2","TP2","TP2","TP2","TP2","TP2","TP2","TP2","TP2","TP2","TP2","TP2","TP2","TP2","TP2","TP2","TP2","TP2","TP2","TP2","TP2","TP2","TP2","TP2","TP2","TP2","TP2","TP2","TP2","TP2","TP2","TP2","TP2","TP2","TP2","TP3","TP3","TP3","TP3","TP3","TP3","TP3","TP3","TP3","TP3","TP3","TP3","TP3","TP3","TP3","TP3","TP3","TP3","TP3","TP3","TP3","TP3","TP3","TP3","TP3","TP3","TP3","TP3","TP3","TP3","TP3","TP3","TP3","TP3","TP3","TP3"],["bjaccard","bjaccard","bjaccard","wjaccard","wjaccard","wjaccard","rAitchison","rAitchison","rAitchison","bjaccard","bjaccard","bjaccard","wjaccard","wjaccard","wjaccard","rAitchison","rAitchison","rAitchison","bjaccard","bjaccard","bjaccard","wjaccard","wjaccard","wjaccard","rAitchison","rAitchison","rAitchison","bjaccard","bjaccard","bjaccard","wjaccard","wjaccard","wjaccard","rAitchison","rAitchison","rAitchison","bjaccard","bjaccard","bjaccard","wjaccard","wjaccard","wjaccard","rAitchison","rAitchison","rAitchison","bjaccard","bjaccard","bjaccard","wjaccard","wjaccard","wjaccard","rAitchison","rAitchison","rAitchison","bjaccard","bjaccard","bjaccard","wjaccard","wjaccard","wjaccard","rAitchison","rAitchison","rAitchison","bjaccard","bjaccard","bjaccard","wjaccard","wjaccard","wjaccard","rAitchison","rAitchison","rAitchison","bjaccard","bjaccard","bjaccard","wjaccard","wjaccard","wjaccard","rAitchison","rAitchison","rAitchison","bjaccard","bjaccard","bjaccard","wjaccard","wjaccard","wjaccard","rAitchison","rAitchison","rAitchison","bjaccard","bjaccard","bjaccard","wjaccard","wjaccard","wjaccard","rAitchison","rAitchison","rAitchison","bjaccard","bjaccard","bjaccard","wjaccard","wjaccard","wjaccard","rAitchison","rAitchison","rAitchison"],["cluster_Dtp2","Residual","Total","cluster_Dtp2","Residual","Total","cluster_Dtp2","Residual","Total","cluster_Dtp2","Residual","Total","cluster_Dtp2","Residual","Total","cluster_Dtp2","Residual","Total","cluster_Dtp2","Residual","Total","cluster_Dtp2","Residual","Total","cluster_Dtp2","Residual","Total","cluster_Dtp2","Residual","Total","cluster_Dtp2","Residual","Total","cluster_Dtp2","Residual","Total","cluster_Dtp2","Residual","Total","cluster_Dtp2","Residual","Total","cluster_Dtp2","Residual","Total","cluster_Dtp2","Residual","Total","cluster_Dtp2","Residual","Total","cluster_Dtp2","Residual","Total","cluster_Dtp2","Residual","Total","cluster_Dtp2","Residual","Total","cluster_Dtp2","Residual","Total","cluster_Dtp2","Residual","Total","cluster_Dtp2","Residual","Total","cluster_Dtp2","Residual","Total","cluster_Dtp2","Residual","Total","cluster_Dtp2","Residual","Total","cluster_Dtp2","Residual","Total","cluster_Dtp2","Residual","Total","cluster_Dtp2","Residual","Total","cluster_Dtp2","Residual","Total","cluster_Dtp2","Residual","Total","cluster_Dtp2","Residual","Total","cluster_Dtp2","Residual","Total","cluster_Dtp2","Residual","Total","cluster_Dtp2","Residual","Total","cluster_Dtp2","Residual","Total"],[2,37,39,2,37,39,2,37,39,2,37,39,2,37,39,2,37,39,2,37,39,2,37,39,2,37,39,2,37,39,2,37,39,2,37,39,2,37,39,2,37,39,2,37,39,2,37,39,2,37,39,2,37,39,2,37,39,2,37,39,2,37,39,2,37,39,2,37,39,2,37,39,2,36,38,2,36,38,2,36,38,2,36,38,2,36,38,2,36,38,2,36,38,2,36,38,2,36,38,2,36,38,2,36,38,2,36,38],[0.3712147280467937,6.961878359790128,7.333093087836922,0.621262564518771,12.87009558818236,13.49135815270113,0.1239099740004366,2.293565130621456,2.417475104621893,0.3712147280467973,6.961878359790128,7.333093087836922,0.6212625645187718,12.87009558818236,13.49135815270113,0.123909974000437,2.293565130621456,2.417475104621893,0.3712147280467937,6.961878359790128,7.333093087836922,0.621262564518771,12.87009558818236,13.49135815270113,0.1239099740004366,2.293565130621456,2.417475104621893,0.3712147280467973,6.961878359790128,7.333093087836922,0.6212625645187718,12.87009558818236,13.49135815270113,0.123909974000437,2.293565130621456,2.417475104621893,0.439030337989144,6.750081118510703,7.189111456499852,1.003562355950793,13.51403700418611,14.5175993601369,0.1400649160122622,2.306833560215958,2.446898476228222,0.4390303379891511,6.750081118510703,7.189111456499852,1.003562355950794,13.51403700418611,14.5175993601369,0.1400649160122639,2.306833560215958,2.446898476228222,0.439030337989144,6.750081118510703,7.189111456499852,1.003562355950793,13.51403700418611,14.5175993601369,0.1400649160122622,2.306833560215958,2.446898476228222,0.4390303379891511,6.750081118510703,7.189111456499852,1.003562355950794,13.51403700418611,14.5175993601369,0.1400649160122639,2.306833560215958,2.446898476228222,0.3519011607085032,7.026999273704757,7.378900434413267,0.6230411045227378,13.19858136337102,13.82162246789376,0.1197523685641819,2.109847637729654,2.229600006293836,0.3519011607085085,7.026999273704757,7.378900434413267,0.6230411045227362,13.19858136337102,13.82162246789376,0.1197523685641818,2.109847637729654,2.229600006293836,0.3519011607085032,7.026999273704757,7.378900434413267,0.6230411045227378,13.19858136337102,13.82162246789376,0.1197523685641819,2.109847637729654,2.229600006293836,0.3519011607085085,7.026999273704757,7.378900434413267,0.6230411045227362,13.19858136337102,13.82162246789376,0.1197523685641818,2.109847637729654,2.229600006293836],[0.05062184859790082,0.949378151402099,1,0.04604892683798382,0.9539510731620162,1,0.05125594624057848,0.9487440537594213,1,0.05062184859790131,0.949378151402099,1,0.04604892683798388,0.9539510731620162,1,0.05125594624057866,0.9487440537594213,1,0.05062184859790082,0.949378151402099,1,0.04604892683798382,0.9539510731620162,1,0.05125594624057848,0.9487440537594213,1,0.05062184859790131,0.949378151402099,1,0.04604892683798388,0.9539510731620162,1,0.05125594624057866,0.9487440537594213,1,0.06106879002302933,0.93893120997697,1,0.06912729378015631,0.9308727062198434,1,0.05724181749794772,0.9427581825020516,1,0.06106879002303032,0.93893120997697,1,0.0691272937801564,0.9308727062198434,1,0.0572418174979484,0.9427581825020516,1,0.06106879002302933,0.93893120997697,1,0.06912729378015631,0.9308727062198434,1,0.05724181749794772,0.9427581825020516,1,0.06106879002303032,0.93893120997697,1,0.0691272937801564,0.9308727062198434,1,0.0572418174979484,0.9427581825020516,1,0.0476901895934695,0.9523098104065295,1,0.04507727699624265,0.9549227230037574,1,0.05371024767946647,0.9462897523205335,1,0.04769018959347023,0.9523098104065295,1,0.04507727699624254,0.9549227230037574,1,0.05371024767946642,0.9462897523205335,1,0.0476901895934695,0.9523098104065295,1,0.04507727699624265,0.9549227230037574,1,0.05371024767946647,0.9462897523205335,1,0.04769018959347023,0.9523098104065295,1,0.04507727699624254,0.9549227230037574,1,0.05371024767946642,0.9462897523205335,1],[0.9864395948843767,null,null,0.8930281336955063,null,null,0.9994634503302527,null,null,0.9864395948843864,null,null,0.8930281336955075,null,null,0.9994634503302561,null,null,0.9864395948843767,null,null,0.8930281336955063,null,null,0.9994634503302527,null,null,0.9864395948843864,null,null,0.8930281336955075,null,null,0.9994634503302561,null,null,1.203253873575843,null,null,1.373823645690677,null,null,1.12327173963269,null,null,1.203253873575863,null,null,1.373823645690679,null,null,1.123271739632704,null,null,1.203253873575843,null,null,1.373823645690677,null,null,1.12327173963269,null,null,1.203253873575863,null,null,1.373823645690679,null,null,1.123271739632704,null,null,0.9014119179513657,null,null,0.8496928247555956,null,null,1.021657960323046,null,null,0.9014119179513792,null,null,0.8496928247555934,null,null,1.021657960323045,null,null,0.9014119179513657,null,null,0.8496928247555956,null,null,1.021657960323046,null,null,0.9014119179513792,null,null,0.8496928247555934,null,null,1.021657960323045,null,null],[0.457,null,null,0.725,null,null,0.434,null,null,0.451,null,null,0.73,null,null,0.462,null,null,0.459,null,null,0.733,null,null,0.459,null,null,0.448,null,null,0.718,null,null,0.47,null,null,0.134,null,null,0.062,null,null,0.176,null,null,0.14,null,null,0.057,null,null,0.184,null,null,0.14,null,null,0.058,null,null,0.158,null,null,0.151,null,null,0.061,null,null,0.167,null,null,0.665,null,null,0.852,null,null,0.38,null,null,0.662,null,null,0.827,null,null,0.395,null,null,0.637,null,null,0.845,null,null,0.368,null,null,0.68,null,null,0.856,null,null,0.38,null,null],["cluster_Dtp2","cluster_Dtp2","cluster_Dtp2","cluster_Dtp2","cluster_Dtp2","cluster_Dtp2","cluster_Dtp2","cluster_Dtp2","cluster_Dtp2","cluster_Dtp2","cluster_Dtp2","cluster_Dtp2","cluster_Dtp2","cluster_Dtp2","cluster_Dtp2","cluster_Dtp2","cluster_Dtp2","cluster_Dtp2","cluster_Dtp2","cluster_Dtp2","cluster_Dtp2","cluster_Dtp2","cluster_Dtp2","cluster_Dtp2","cluster_Dtp2","cluster_Dtp2","cluster_Dtp2","cluster_Dtp2","cluster_Dtp2","cluster_Dtp2","cluster_Dtp2","cluster_Dtp2","cluster_Dtp2","cluster_Dtp2","cluster_Dtp2","cluster_Dtp2","cluster_Dtp2","cluster_Dtp2","cluster_Dtp2","cluster_Dtp2","cluster_Dtp2","cluster_Dtp2","cluster_Dtp2","cluster_Dtp2","cluster_Dtp2","cluster_Dtp2","cluster_Dtp2","cluster_Dtp2","cluster_Dtp2","cluster_Dtp2","cluster_Dtp2","cluster_Dtp2","cluster_Dtp2","cluster_Dtp2","cluster_Dtp2","cluster_Dtp2","cluster_Dtp2","cluster_Dtp2","cluster_Dtp2","cluster_Dtp2","cluster_Dtp2","cluster_Dtp2","cluster_Dtp2","cluster_Dtp2","cluster_Dtp2","cluster_Dtp2","cluster_Dtp2","cluster_Dtp2","cluster_Dtp2","cluster_Dtp2","cluster_Dtp2","cluster_Dtp2","cluster_Dtp2","cluster_Dtp2","cluster_Dtp2","cluster_Dtp2","cluster_Dtp2","cluster_Dtp2","cluster_Dtp2","cluster_Dtp2","cluster_Dtp2","cluster_Dtp2","cluster_Dtp2","cluster_Dtp2","cluster_Dtp2","cluster_Dtp2","cluster_Dtp2","cluster_Dtp2","cluster_Dtp2","cluster_Dtp2","cluster_Dtp2","cluster_Dtp2","cluster_Dtp2","cluster_Dtp2","cluster_Dtp2","cluster_Dtp2","cluster_Dtp2","cluster_Dtp2","cluster_Dtp2","cluster_Dtp2","cluster_Dtp2","cluster_Dtp2","cluster_Dtp2","cluster_Dtp2","cluster_Dtp2","cluster_Dtp2","cluster_Dtp2","cluster_Dtp2"],["margin","margin","margin","margin","margin","margin","margin","margin","margin","terms","terms","terms","terms","terms","terms","terms","terms","terms","margin","margin","margin","margin","margin","margin","margin","margin","margin","terms","terms","terms","terms","terms","terms","terms","terms","terms","margin","margin","margin","margin","margin","margin","margin","margin","margin","terms","terms","terms","terms","terms","terms","terms","terms","terms","margin","margin","margin","margin","margin","margin","margin","margin","margin","terms","terms","terms","terms","terms","terms","terms","terms","terms","margin","margin","margin","margin","margin","margin","margin","margin","margin","terms","terms","terms","terms","terms","terms","terms","terms","terms","margin","margin","margin","margin","margin","margin","margin","margin","margin","terms","terms","terms","terms","terms","terms","terms","terms","terms"]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>Time<\/th>\n      <th>Distance<\/th>\n      <th>terms<\/th>\n      <th>Df<\/th>\n      <th>SumOfSqs<\/th>\n      <th>R2<\/th>\n      <th>F<\/th>\n      <th>Pr(&gt;F)<\/th>\n      <th>Group<\/th>\n      <th>terms_margin<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":[4,5,6,7,8]},{"orderable":false,"targets":0},{"name":" ","targets":0},{"name":"Time","targets":1},{"name":"Distance","targets":2},{"name":"terms","targets":3},{"name":"Df","targets":4},{"name":"SumOfSqs","targets":5},{"name":"R2","targets":6},{"name":"F","targets":7},{"name":"Pr(>F)","targets":8},{"name":"Group","targets":9},{"name":"terms_margin","targets":10}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script>
```


``` r
save(Time_out_plaque_clust, Time_out_saliva_clust, Time_out_clust, Time_out, plaque, saliva, d_box, d_box_leg, 
     ps_up,
     sample_pal,
     time_pal, sub_pal, sex_pal, sample_pal, period_pal,
     file = here::here("../../data/processed_data/metaphlan/03_data.Rdata"))
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
##  [4] vegan_2.6-8          lattice_0.22-6       permute_0.9-7       
##  [7] GUniFrac_1.8         ape_5.8              reshape2_1.4.4      
## [10] scales_1.3.0         microViz_0.12.4      speedyseq_0.5.3.9021
## [13] phyloseq_1.50.0      readxl_1.4.3         lubridate_1.9.3     
## [16] forcats_1.0.0        stringr_1.5.1        dplyr_1.1.4         
## [19] purrr_1.0.2          readr_2.1.5          tidyr_1.3.1         
## [22] tibble_3.2.1         ggplot2_3.5.1        tidyverse_2.0.0     
## 
## loaded via a namespace (and not attached):
##   [1] RColorBrewer_1.1-3      rstudioapi_0.17.1       jsonlite_1.8.9         
##   [4] magrittr_2.0.3          farver_2.1.2            rmarkdown_2.29         
##   [7] zlibbioc_1.51.2         vctrs_0.6.5             multtest_2.61.0        
##  [10] htmltools_0.5.8.1       broom_1.0.7             cellranger_1.1.0       
##  [13] Rhdf5lib_1.27.0         Formula_1.2-5           rhdf5_2.49.0           
##  [16] sass_0.4.9              bslib_0.8.0             htmlwidgets_1.6.4      
##  [19] plyr_1.8.9              plotly_4.10.4           cachem_1.1.0           
##  [22] igraph_2.1.1            lifecycle_1.0.4         iterators_1.0.14       
##  [25] pkgconfig_2.0.3         Matrix_1.7-1            R6_2.5.1               
##  [28] fastmap_1.2.0           GenomeInfoDbData_1.2.13 clue_0.3-65            
##  [31] digest_0.6.37           colorspace_2.1-1        GGally_2.2.1           
##  [34] spatial_7.3-17          S4Vectors_0.43.2        rprojroot_2.0.4        
##  [37] crosstalk_1.2.1         labeling_0.4.3          fansi_1.0.6            
##  [40] timechange_0.3.0        httr_1.4.7              abind_1.4-8            
##  [43] mgcv_1.9-1              compiler_4.4.0          here_1.0.1             
##  [46] withr_3.0.2             backports_1.5.0         inline_0.3.19          
##  [49] carData_3.0-5           ggstats_0.7.0           highr_0.11             
##  [52] ggsignif_0.6.4          MASS_7.3-61             biomformat_1.34.0      
##  [55] fBasics_4041.97         tools_4.4.0             glue_1.8.0             
##  [58] stabledist_0.7-2        nlme_3.1-166            rhdf5filters_1.17.0    
##  [61] grid_4.4.0              Rtsne_0.17              cluster_2.1.6          
##  [64] ade4_1.7-22             generics_0.1.3          gtable_0.3.6           
##  [67] microbiome_1.28.0       tzdb_0.4.0              data.table_1.16.2      
##  [70] hms_1.1.3               car_3.1-3               utf8_1.2.4             
##  [73] XVector_0.45.0          rmutil_1.1.10           BiocGenerics_0.52.0    
##  [76] ggrepel_0.9.6           foreach_1.5.2           pillar_1.9.0           
##  [79] splines_4.4.0           survival_3.7-0          tidyselect_1.2.1       
##  [82] Biostrings_2.73.2       knitr_1.48              IRanges_2.39.2         
##  [85] stats4_4.4.0            xfun_0.49               Biobase_2.66.0         
##  [88] statmod_1.5.0           timeDate_4041.110       matrixStats_1.4.1      
##  [91] DT_0.33                 stringi_1.8.4           UCSC.utils_1.2.0       
##  [94] lazyeval_0.2.2          yaml_2.3.10             evaluate_1.0.1         
##  [97] codetools_0.2-20        timeSeries_4041.111     cli_3.6.3              
## [100] rpart_4.1.23            munsell_0.5.1           jquerylib_0.1.4        
## [103] Rcpp_1.0.13-1           GenomeInfoDb_1.42.0     stable_1.1.6           
## [106] parallel_4.4.0          modeest_2.4.0           viridisLite_0.4.2      
## [109] statip_0.2.3            crayon_1.5.3            rlang_1.1.4            
## [112] cowplot_1.1.3
```
