#' Import and Process HUMAnN Pathway Data for Microbial Ecology Analysis
#'
#' This function imports HUMAnN pathway data, processes it into a `phyloseq` object, 
#' and cleans and annotates the data to facilitate downstream analysis in microbial ecology studies.
#' 
#' @param pathway A file path or object containing HUMAnN pathway or gene family level data.
#' @param rank The taxonomic rank used for filtering or organizing the taxonomy table. Default is `"pathway"`.
#' @param tax_clean_gsub A regular expression pattern to clean unwanted characters or prefixes in taxonomy labels. Default is `"[a-t]__"`.
#' @param no_sel A vector of values to exclude during filtering, such as `NA` or `"unclassified"`. Default is `c(NA, "<NA>", "NA", "unclassified")`.
#' @param sub_pat_a A string pattern for initial cleaning of sample names. Default is `"_trimmed_HF_merged_Abundance-RELAB"`.
#' @param sub_pat_b A string pattern for additional cleaning of sample names. Default is `"_subset50"`.
#' @param tax_na_if A string indicating taxonomy entries to replace with `NA`. Default is `"unclassified"`.
#' @param sample_name_clean_a A regular expression pattern for extracting meaningful parts of sample names. Default is `"[^_]+"`.
#' @param sample_name_clean_b A regular expression for removing specific patterns in sample names. Default is `"[A-Z]"`.
#' @param meta File path to metadata (in Excel format). Default points to a metadata file for the Deerland Study.
#' @param sample_column The column name in the metadata file corresponding to sample IDs. Default is `"ID"`.
#' @param sample_nam_prefix A prefix to add to sample names. Default is `"S_"`.
#' @param rm_un Logical; whether to remove taxonomy entries starting with `"UN"`. Default is `TRUE`.
#' @param return_unstrat Logical; whether to return only unstratified pathways. Default is `TRUE`.
#' @param str_trun An integer specifying the maximum length of pathway names for truncation. Default is `40`.
#' @param add_MetaCyc_pathway_map Logical; whether to include MetaCyc pathway mapping. Default is `TRUE`.
#' @param add_CHOCOPhlAn_taxonomy Logical; whether to annotate with CHOCOPhlAn taxonomy. Default is `TRUE`.
#' @param final_tax_table_ranks A vector of column names defining the final structure of the taxonomy table. Default includes taxonomic and pathway features.
#' 
#' @return A `phyloseq` object with processed HUMAnN pathway data, cleaned sample names, taxonomy table, and metadata.
#' 
#' @details 
#' This function performs the following operations:
#' - Imports HUMAnN pathway data using `mia::importHUMAnN`.
#' - Converts the data into a `phyloseq` object for downstream microbial ecology analysis.
#' - Cleans sample names based on specified substitution patterns and regular expressions.
#' - Adds metadata to the `phyloseq` object for sample-level annotations.
#' - Processes and cleans the taxonomy table, replacing "unclassified" entries with `NA`, and optionally removing entries starting with "UN".
#' - Annotates pathways using MetaCyc and/or CHOCOPhlAn taxonomy data, with abbreviations for long pathway names and features.
#' 
#' @examples
#' \dontrun{
#' pathway_data <- "path/to/humann_pathway_data.tsv"
#' processed_data <- import_humann_pathway_mia(
#'   pathway = pathway_data,
#'   meta = "path/to/metadata.xlsx",
#'   sample_column = "SampleID"
#' )
#' }
#' 
#' @import microeco file2meco mia phyloseq microViz stringr dplyr readxl here
#' @export
#' 
#' 
import_humann_pathway_mia <- function(pathway = NULL,
                                      rank = "pathway",
                                      tax_clean_gsub = "[a-t]__",
                                      no_sel = c(NA, "<NA>","NA", "unclassified"),
                                      sub_pat_a = "_trimmed_HF_merged_Abundance-RELAB",
                                      sub_pat_b = "_subset50",
                                      tax_na_if = "unclassified",
                                      sample_name_clean_a = "[^_]+",
                                      sample_name_clean_b = "[A-Z]",
                                      meta = here::here("../../data/processed_data/metaphlan/Metadata_Deerland_Study.xlsx") %>% 
                                        read_excel(sheet = "Trial2"),
                                      sample_column = "ID",
                                      sample_nam_prefix = "S_",
                                      rm_un = TRUE,
                                      return_unstrat = TRUE,
                                      str_trun = 40,
                                      add_MetaCyc_pathway_map = TRUE,
                                      add_CHOCOPhlAn_taxonomy = TRUE,
                                      final_tax_table_ranks = c("Superclass1", "Superclass2", "feature", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")){
  
  
  suppressMessages(suppressWarnings({
    # Load required libraries
    require(microeco); require(file2meco)# For ecological and microbial ecology data analysis
    # require(microViz)   # An extension of phyloseq for speedy data manipulation
    
    pathway %>% 
      mia::importHUMAnN()  %>% 
      mia::convertToPhyloseq() %>% 
      # Clean sample names by applying various substitution patterns iteratively
      clean_phyloseq_sample_names(sub_pat = sub_pat_a) %>% 
      clean_phyloseq_sample_names(sub_pat = sub_pat_b) -> pw_mia_ps
    
    colnames(tax_table(pw_mia_ps)) = str_to_title(rank_names(pw_mia_ps))
    # taxa_sums(pw_mia_ps) %>%  sort(decreasing = T) %>%  tail()
    
    # Extract and clean sample names
    pw_mia_ps %>% 
      sample_names() %>% 
      str_extract(., sample_name_clean_a) %>%  # Extract meaningful parts of sample names
      str_remove_all(., sample_name_clean_b) -> sample_names(pw_mia_ps)
    
    # Add metadata to the phyloseq object
    pw_mia_ps %>% 
      physeq_add_metadata(metadata = meta, 
                          sample_column = sample_column) -> pw_mia_ps
    
    # Add a prefix to sample names
    pw_mia_ps %>% 
      sample_names() %>% 
      paste0(sample_nam_prefix, .) -> sample_names(pw_mia_ps)
    
    
    # Clean up the taxonomy table by removing unwanted patterns
    tax_table(pw_mia_ps) <- tax_table(pw_mia_ps) %>% 
      gsub(pattern = tax_clean_gsub, replacement = "") %>%  # Remove specific patterns
      # as.matrix() %>%  
      tax_table()
    
    # Replace "unclassified" entries with NA in the taxonomy table
    tax_table(pw_mia_ps) <- tax_table(pw_mia_ps) %>% 
      data.frame() %>% 
      mutate_all(na_if, tax_na_if) %>%  # Convert "unclassified" to NA
      as.matrix() %>%  
      tax_table()
    
    # Filter taxonomy table based on conditions (commented lines suggest alternative filters)
    pw_mia_ps %>%
      microViz::tax_mutate(feature = taxa_names(.)) -> pw_mia_ps
    
    # rank="feature"
    # pw_mia_ps %>% 
    #   speedyseq::filter_tax_table(get(rank) %!in% eval(as.character(no_sel))) -> pw_mia_ps
    # tax_mutate(Strain = NULL) %>%
    # microViz::tax_select(ranks_searched = rank, tax_list = eval(as.character(no_sel)),
    # deselect = TRUE) -> pw_me_ps_2              # Optional taxonomic selection step
    ## speedyseq::filter_tax_table(!eval(as.name((rank))) %in% no_sel) -> pw_me_ps# Exclude unwanted taxonomic ranks
    ## speedyseq::filter_tax_table(is.na("pathway")) %>%     # Alternative: filter NA pathways
    
    # pw_me_ps %>% 
    #   subset_taxa(is.na(Genus)) -> pw_me_ps
    
    if(return_unstrat)
    {  
      rank="Genus"
      no_sel=c(NA,"NA", "<NA>", NA_character_)
      
      pw_mia_ps %>% 
        speedyseq::filter_tax_table(get(rank) %in% eval(as.character(no_sel))) %>% 
        subset_taxa(!(grepl("unclassified", feature))) -> pw_mia_ps
    }
    
    if(rm_un)
    {
      pw_mia_ps %>%
        subset_taxa(!(grepl("^UN", feature)))  -> pw_mia_ps
      
      pw_mia_ps %>%
        microViz::tax_mutate(feature_code = taxa_names(.)) -> pw_mia_ps
      # Return the cleaned and updated phyloseq object
    }
    
    if(isTRUE(rm_un) & isTRUE(return_unstrat)){
      taxa_names(pw_mia_ps) <- taxa_names(pw_mia_ps) %>% 
        stringi::stri_extract(., regex = '[^:]*') 
    }
    
    if(add_MetaCyc_pathway_map)
    {
      data(MetaCyc_pathway_map)  
      
      pw_mia_ps %>%
        microViz::tax_mutate(feature_code = taxa_names(.)) -> pw_mia_ps
      
      pw_mia_ps %>% 
        tax_table() %>% 
        as.data.frame() %>% 
        rownames_to_column("tmp") %>% 
        # select(feature) %>% 
        dplyr::left_join(.,
                         MetaCyc_pathway_map %>% 
                           rownames_to_column("tmp") %>% 
                           mutate_all(funs(str_replace_all(., "&&", "-"))) %>% 
                           mutate_all(funs(str_replace_all(., "[&$|;]", "-"))),
                         by = c("tmp" = "tmp") 
        ) %>%  select(one_of(c("tmp", colnames(MetaCyc_pathway_map), rank_names(pw_mia_ps)))) %>% 
        column_to_rownames('tmp') %>% 
        mutate(feature = 
                 stringr::str_replace_all(pathway, 
                                          c(biosynthesis="Bios.", 
                                            vitamin="Vit.",
                                            metabolism = "Met.",
                                            nucleotides = "Nucleot.",
                                            nucleosides = "Nucleos.",
                                            carbohydrate = "Carb.",
                                            glycolysis = "Glycol.",
                                            assimilation = "Assim.",
                                            degradation = "Degrad.",
                                            `amino acid` = "AA",
                                            superpathway = "SPWY.",
                                            fermentation = "Ferm.",
                                            nutrient = "Nutr.",
                                            photosynthesis = "Photosynth.",
                                            pathways = "PWY.",
                                            degradation = "Deg.",
                                            utilization = "Util.",
                                            structure = "Stru.",
                                            vitamin = "Vit.",
                                            `fatty Acid` = "F.Acid",
                                            inorganic = "Inorg.",
                                            aromatic = "Aroma.",
                                            compound = "Comp.",
                                            decarboxylation = "decarboxyl.",
                                            reduction = "Red.",
                                            transformations = "transf.",
                                            photosynthetic = "Photosynth."))) %>% 
        mutate(across(c("Superclass1", "Superclass2"), 
                      ~ stringr::str_replace_all(., c(Precursor = "Precurs.",
                                                      Biosynthesis="Bios.", 
                                                      Transfer = "Trans.",
                                                      Vitamin="Vit.",
                                                      Generation = "Gen.",
                                                      Metabolism = "Met.",
                                                      Nucleotides = "Nucleot.",
                                                      Nucleosides = "Nucleos.",
                                                      Carbohydrate = "Carb.",
                                                      Glycolysis = "Glycol.",
                                                      Assimilation = "Assim.",
                                                      Degradation = "Degrad.",
                                                      `Amino acid` = "AA",
                                                      Superpathway = "SPWY.",
                                                      Fermentation = "Ferm.",
                                                      Nutrient = "Nutr.",
                                                      Photosynthesis = "Photosynth.",
                                                      Pathways = "PWY.",
                                                      Degradation = "Deg.",
                                                      Utilization = "Util.",
                                                      Structure = "Stru.",
                                                      Vitamin = "Vit.",
                                                      `Fatty Acid` = "F.Acid",
                                                      Inorganic = "Inorg.",
                                                      Aromatic = "Aroma.",
                                                      Compound = "Comp.",
                                                      Decarboxylation = "decarboxyl.",
                                                      Reduction = "Red.",
                                                      Transformations = "transf.",
                                                      Photosynthetic = "Photosynth.")
                      ))) %>% 
        mutate(across(c("Superclass1", "Superclass2","pathway"), 
                      ~ stringr::str_trunc(., width = str_trun, side ="center"))) %>% 
        select(one_of(final_tax_table_ranks)) %>% 
        as.matrix() -> tax_table(pw_mia_ps)
    }
    
    if(add_CHOCOPhlAn_taxonomy){
      data(CHOCOPhlAn_taxonomy)
      
      CHOCOPhlAn_taxonomy %>%
        as.matrix() %>%
        gsub(pattern = tax_clean_gsub, replacement = "") %>%  # Remove specific patterns
        as.data.frame() -> CHOCOPhlAn_taxonomy
      
      pw_mia_ps %>% 
        tax_table() %>% 
        as.data.frame() %>% 
        rownames_to_column("tmp") %>% 
        dplyr::left_join(.,
                         CHOCOPhlAn_taxonomy,
                         # rownames_to_column("tmp"),
                         by = c("Genus" = "Genus")) %>% 
        select(one_of(c("tmp",final_tax_table_ranks))) %>% 
        column_to_rownames('tmp') %>% 
        as.matrix() -> tax_table(pw_mia_ps)
      
    }
    
    if(!isTRUE(add_CHOCOPhlAn_taxonomy) & !isTRUE(add_MetaCyc_pathway_map)){
      pw_mia_ps %>% 
        physeq_sel_tax_table(final_tax_table_ranks) -> pw_mia_ps
    }
    
  }))
  return(pw_mia_ps)
}



#' @title Differiential feature analyses wrapper for Microbiome analyses
#' @author Florentin Constancias
#' @description 
#' This script performs a comparative analysis of various approaches for microbiome studies, 
#' including metagenomics, metabarcoding, metabolomics, and transcriptomics. The script 
#' calculates the number of publications, environments analyzed, and the methods applied.
#' 
#' @param ps_tmp A phyloseq object. Defaults to saliva samples (`ps_up %>% subset_samples(Sample == "Saliva")`).
#' @param approach Character vector specifying the analysis approach. Options include:
#'   - "run_lefse": Linear discriminant analysis (LEfSe)
#'   - "run_ancom": Analysis of composition of microbiomes (ANCOM)
#'   - "ancombc2": ANCOMBC2 method for differential abundance
#'   - "maaslin3": Multivariate analysis using Maaslin3
#'   - "trans_diff_rf": Transformation-based Random Forest
#'   - "classifier_rf": Random Forest classifier
#' @param glom Taxonomic rank for agglomeration (e.g., "Species"). Set to `NULL` to skip agglomeration.
#' @param unclassified_name Name for unclassified taxa (default: "UNCLASSIFIED").
#' @param taxa_rank Taxonomic rank for differential abundance (default: "all").
#' @param density Density metric to normalize counts. Default is "Quant".
#' @param comp_group Grouping variable for comparison (e.g., "Time").
#' @param palette Color palette for plots.
#' @param pvalue_cutoff p-value cutoff for significance.
#' @param p_adjust Method for p-value adjustment (e.g., "BH" for Benjamini-Hochberg).
#' @param lefse_* Parameters specific to the LEfSe method.
#' @param maaslin3_* Parameters specific to Maaslin3 analysis.
#' @param ancombc2_* Parameters specific to ANCOMBC2 analysis.
#' @param linda_* Parameters specific to LINDA (if used in the future).
#' @param rf_* Parameters specific to Random Forest methods.
#' @param trans_diff_rf_MeanDecreaseGini_cutoff Gini index cutoff for feature importance in Random Forest.
#' @import tidyverse dplyr ggplot2

#' @return A list containing results of the selected approaches, including plots, tables, and significance testing results.
#' @note To improve this function, consider outputting all results in a structured list for easier downstream processing.
#' @examples
#' # Run LEfSe and ANCOM on saliva samples
#' results <- phyloseq_diff_abundance(ps_tmp, approach = c("run_lefse", "run_ancom"))
#' 
#' # Access LEfSe results
#' results$mmlefse

phyloseq_diff_abundance <- function(ps_tmp = ps_up %>%  subset_samples(Sample == "Saliva"),
                                    approach = c("run_lefse", 
                                                 "run_ancom",
                                                 "ancombc2",
                                                 "maaslin3",
                                                 "trans_diff_rf",
                                                 "classifier_rf"),
                                    glom = "Species",
                                    unclassified_name = "UNCLASSIFIED",
                                    taxa_rank = "all", #"OTU"
                                    density = "Quant",
                                    comp_group = "Time",
                                    
                                    plot_top_sig = NULL,
                                    
                                    palette = time_pal,
                                    pvalue_cutoff = 0.05,
                                    p_adjust = "BH",
                                    
                                    
                                    lefse_taxa_rank = taxa_rank,
                                    
                                    lefse_prv_cut = 0.2,
                                    lefse_sample_min = 10,
                                    lefse_lda_cutoff = 2,
                                    lefse_multigrp_strat = FALSE,
                                    lefse_strict = "0",
                                    lefse_bootstrap_n = 30,
                                    
                                    maaslin3_formula =  formula,
                                    maaslin3_fixed_effects = NULL,
                                    maaslin3_strata_effects = NULL,
                                    masslin3_output_dir = "~/test_masslin3/",
                                    maaslin3_prv_cut = 0,
                                    maaslin3_lefse_plot = TRUE,
                                    
                                    formula =  "~ Time +  (1|Subject)",
                                    
                                    ancom_prv_cut = 0.1,
                                    ancom_confounders = character(0), # Subject
                                    
                                    ancombc2_forumla = formula,
                                    ancombc2_fix_formula = "Time",
                                    ancombc2_rand_formula = NULL,
                                    ancombc2_group = "Time",
                                    ancombc2_pairwise = TRUE,
                                    ancombc2_dunnet = TRUE,
                                    ancombc2_trend = TRUE,
                                    ancombc2_global = TRUE,
                                    
                                    linda_prv_cut = 0.1,
                                    linda_formula = formula,
                                    linda_comp_group = comp_group,
                                    
                                    rf_prv_cut = 0.1,
                                    rf_prop_train = 3/4,
                                    
                                    ref_train_max_mtry=3,
                                    ref_train_ntree = c(100, 500, 1000),
                                    feature_imp_nrep = 99,
                                    
                                    trans_diff_rf_MeanDecreaseGini_cutoff = 0){
  # To improve: store in list output of test so we can plot those in a loop or laply manner
  #
  # - Lefse
  # - ANCOMBC (en count et pas %)
  # - linda (count? percent?)
  # - maaslin3 (avec lefse plot) - check maaslin3 et low replication
  # 
  # https://microbiome.github.io/course_2022_radboud/differential-abundance-analysis-demo.html
  # Boxplot of the significant taxa + heatmap + vulcano?
  # Ideallement un modele global Sample Type + Time + 1:Subject
  # puis lefse sur les significant pour les pairwise comp.
  
  
  ########## ----- Source function and load packages
  require(microbiomeMarker);require(tidyverse);require(magrittr);require(ANCOMBC);require(microeco)
  
  source("https://raw.githubusercontent.com/fconstancias/DivComAnalyses/master/R/phyloseq_taxa_tests.R")
  source("https://raw.githubusercontent.com/fconstancias/DivComAnalyses/master/R/phyloseq_microeco.R")  
  
  suppressMessages(suppressWarnings({
    
    ########## ----- 
    
    out <- list()
    sumfilter = 0
    
    ps_tmp %>% 
      filter_taxa(function(x) sum(x > 0) > 0, TRUE) -> ps_tmp
    
    # ps_tmp %>%
    # microbiome::core(detection = 0, prevalence = prv_cut) #-> ps_filtered
    
    # ps_tmp %>%
    # filter_taxa(function(x){sum(x > sumfilter) >  prv_cut*nsamples(ps_tmp)}, prune = TRUE) -> ps_tmp
    
    ########## ----- Agglomerate Taxa
    
    if(!is.null(glom))
    {
      ps_tmp %>% 
        physeq_glom_rename(speedyseq = TRUE, taxrank = glom, rename_ASV = glom) -> ps_tmp
    }
    
    ########## ----- Prepare count object
    
    if(!is.null(density))
    {
      ps_tmp %>%  
        phyloseq_density_normalize(value_idx = density) -> ps_count
    }
    
    ########## ----- 
    
    if ("run_lefse"  %in% approach)
    {
      ########## ----- microbiomeMarker::run_lefse see:  https://github.com/yiluheihei/microbiomeMarker/issues/95
      # increase permutation to 999
      
      ps_tmp %>%
        physeq_sel_tax_table(c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")) %>% 
        filter_tax_table(Kingdom != unclassified_name) %>% 
        # subset_taxa(Kingdom != unclassified_name) %>% 
        transform_sample_counts(function(x) x/sum(x) * 100) %>% 
        filter_taxa(function(x){sum(x > 0) >  lefse_prv_cut*nsamples(ps_tmp)}, prune = TRUE) %>% 
        microbiomeMarker::run_lefse(group = comp_group, 
                                    norm = "CPM",
                                    wilcoxon_cutoff = pvalue_cutoff,
                                    kw_cutoff = pvalue_cutoff,
                                    taxa_rank = taxa_rank, 
                                    sample_min = lefse_sample_min,
                                    multigrp_strat = lefse_multigrp_strat,
                                    lda_cutoff = lefse_lda_cutoff,
                                    bootstrap_n =  lefse_bootstrap_n,
                                    strict = lefse_strict) -> out$mmlefse
      
      microbiomeMarker::plot_abundance(out$mmlefse, group = comp_group) -> out$mmlefse_abp
      
      # out$mmlefse_abp$data %>% 
      #   mutate(prop = abd / 1e+6) %>% 
      #   ggplot()
      
      # reoder plot so it is in the same order of lefse see: https://stackoverflow.com/questions/12774210/how-do-you-specifically-order-ggplot2-x-axis-instead-of-alphabetical-order
      
      microbiomeMarker::plot_ef_bar(out$mmlefse) + 
        scale_fill_manual(values = palette, na.value = "grey10") -> out$mmlefse_p
      
      microbiomeMarker::marker_table(out$mmlefse) %>%
        data.frame() -> out$mmlefse_df
      
      # save plot, mmlesfse and mmlefse_df
    }
    
    if ("lefse"  %in% approach)
    {
      ########## ----- microeco lefse https://chiliubio.github.io/microeco_tutorial/model-based-class.html
      
      ps_tmp %>%
        physeq_sel_tax_table(c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")) %>% 
        # subset_taxa(Kingdom != unclassified_name) %>% 
        filter_tax_table(Kingdom != unclassified_name) %>% 
        transform_sample_counts(function(x) x/sum(x) * 100) %>% 
        filter_taxa(function(x){sum(x > 0) >  lefse_prv_cut*nsamples(ps_tmp)}, prune = TRUE) %>% 
        file2meco::phyloseq2meco(.) -> data
      
      t1 <- trans_diff$new(dataset = data, 
                           lefse_min_subsam = lefse_sample_min, 
                           lefse_sub_strict = lefse_multigrp_strat, 
                           lefse_sub_alpha = NULL, 
                           lefse_bootstrap_n= lefse_bootstrap_n, 
                           boots = lefse_bootstrap_n,
                           group_choose_paired = NULL,
                           method = "lefse",
                           filter_thres = 0, #default 0; the abundance threshold,
                           alpha = pvalue_cutoff,
                           group = comp_group, 
                           lefse_taxa_rank = taxa_rank, 
                           p_adjust_method = p_adjust,
                           plot_pal = palette,
                           n.cores = 4,
                           lefse_subgroup = NULL)
      
      out$lefse$res_diff  <-  t1$res_diff
      
      if(is.null(plot_top_sig)){
        t1$res_diff %>%  nrow() -> plot_top_sig
        # print("toto")
      }
      
      
      # out$lefse$res_diff <- t1$res_diff
      
      out$lefse$g1 <- t1$plot_diff_bar(use_number = 1:plot_top_sig, color_values = palette, threshold = lefse_lda_cutoff)
      # plot the abundance using same taxa in g
      out$lefse$g2 <- t1$plot_diff_abund(select_taxa = t1$plot_diff_bar_taxa, plot_type = "barerrorbar",
                                         add_sig = TRUE, errorbar_addpoint = FALSE, 
                                         errorbar_color_black = TRUE, 
                                         color_values = palette) # barerrorbar  ggboxplot
      # now the y axis in g1 and g2 is same, so we can merge them
      # remove g1 legend; remove g2 y axis text and ticks
      out$lefse$g1 <- out$lefse$g1 + theme(legend.position = "none")
      out$lefse$g2 <- out$lefse$g2 + theme(axis.text.y = element_blank(), 
                                           axis.ticks.y = element_blank(), 
                                           panel.border = element_blank()) + scale_y_continuous(trans='sqrt')
      out$lefse$p <- out$lefse$g1 %>% aplot::insert_right(out$lefse$g2 )
      
      out$lefse$all <- clone(t1)
      
      
    }
    
    
    
    if ("run_ancom"  %in% approach)
    {
      
      ########## ----- microbiomeMarker::run_ancom 
      
      ps_count %>%
        physeq_sel_tax_table(c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")) %>% 
        # subset_taxa(Kingdom != unclassified_name) %>% 
        filter_tax_table(Kingdom != unclassified_name) %>% 
        filter_taxa(function(x){sum(x > 0) >  ancom_prv_cut*nsamples(ps_tmp)}, prune = TRUE) %>%   
        microbiomeMarker::run_ancom(group = comp_group, 
                                    pvalue_cutoff = pvalue_cutoff,
                                    taxa_rank = taxa_rank, 
                                    confounders = ancom_confounders,
                                    p_adjust = p_adjust) -> out$mmancom
      
      microbiomeMarker::plot_ef_bar(out$mmancom) + 
        scale_fill_manual(values = palette, na.value = "grey10") -> out$mmancom_p
      
      microbiomeMarker::marker_table(out$mmancom) %>%
        data.frame() -> out$mmancom_df
      
      microbiomeMarker::plot_abundance(out$mmancom, group = comp_group) -> out$mmancom_abp
      
      # out$mmlefse_abp$data %>% 
      #   mutate(prop = abd / 1e+6) %>% 
      #   ggplot()
      
      # reoder plot so it is in the same order of lefse see: https://stackoverflow.com/questions/12774210/how-do-you-specifically-order-ggplot2-x-axis-instead-of-alphabetical-order
      
    }
    
    ########## ----- microbiomeMarker::run_ancombc see: https://github.com/yiluheihei/microbiomeMarker/issues/97
    
    # ps_count %>%
    #   physeq_sel_tax_table(c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")) %>% 
    #   subset_taxa(Kingdom != unclassified_name) %>% 
    #   # subset_taxa(Order != "unassigned") %>% 
    #   microbiomeMarker::run_ancombc(group = comp_group, pvalue_cutoff = pvalue_cutoff,
    #                               taxa_rank = taxa_rank, 
    #                               confounders = ancom_confounders,
    #                               p_adjust = p_adjust) -> mmancombc
    # 
    # microbiomeMarker::plot_ef_bar(mmancombc) + 
    #   scale_fill_manual(values = palette, na.value = "grey10") -> mmancombc_p
    # 
    # microbiomeMarker::marker_table(mmancombc) %>%
    #   data.frame() -> mmancombc_df
    
    
    ########## ----- microbiomeMarker::run_ancombc see: https://github.com/yiluheihei/microbiomeMarker/issues/97
    
    # ps_count %>%
    #   physeq_sel_tax_table(c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")) %>% 
    #   subset_taxa(Kingdom != unclassified_name) %>% 
    #   # subset_taxa(Order != "unassigned") %>% 
    #   microbiomeMarker::run_ancombc(group = comp_group, pvalue_cutoff = pvalue_cutoff,
    #                               taxa_rank = taxa_rank, 
    #                               confounders = ancom_confounders,
    #                               p_adjust = p_adjust) -> mmancombc
    # 
    # microbiomeMarker::plot_ef_bar(mmancombc) + 
    #   scale_fill_manual(values = palette, na.value = "grey10") -> mmancombc_p
    # 
    # microbiomeMarker::marker_table(mmancombc) %>%
    #   data.frame() -> mmancombc_df
    
    if ("ancombc2"  %in% approach)
    {
      ########## ----- ancombc2 https://www.bioconductor.org/packages/release/bioc/vignettes/ANCOMBC/inst/doc/ANCOMBC2.html
      
      ps_count %>% # accepts count & 
        physeq_sel_tax_table(c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")) %>% 
        # subset_taxa(Kingdom != unclassified_name) %>% 
        filter_tax_table(Kingdom != unclassified_name) %>% 
        physeq_glom_rename(speedyseq = TRUE, taxrank = ifelse(taxa_rank == "all", "Species", taxa_rank), rename_ASV = ifelse(taxa_rank == "all", "Species", taxa_rank)) %>%
        ANCOMBC::ancombc2(data = ., tax_level = ifelse(taxa_rank == "all", "Species", taxa_rank), 
                          group = ancombc2_group , pairwise = ancombc2_pairwise, #group = "Time",
                          fix_formula = ancombc2_fix_formula, 
                          rand_formula = ancombc2_rand_formula, 
                          p_adj_method = p_adjust, prv_cut = 0, #lib_cut = 1000, 
                          struc_zero = TRUE, neg_lb = TRUE, alpha = pvalue_cutoff, 
                          n_cl = 2, verbose = TRUE,
                          dunnet = ancombc2_dunnet, trend = ancombc2_trend,global = ancombc2_global,
                          iter_control = list(tol = 1e-2, max_iter = 20, 
                                              verbose = TRUE),
                          em_control = list(tol = 1e-5, max_iter = 100),
                          lme_control = lme4::lmerControl(),
                          mdfdr_control = list(fwer_ctrl_method = "holm", B = 100),
                          trend_control = list(contrast = list(matrix(c(1, 0, -1, 1),
                                                                      nrow = 2, 
                                                                      byrow = TRUE)),
                                               node = list(2),
                                               solver = "ECOS",
                                               B = 10)) -> out$ancombc2
      # ancombc2$res %>%
      #   mutate_if(is.numeric, function(x) round(x, 2))
      # 
      # out$ancombc2$res_global %>%
      # mutate_if(is.numeric, function(x) round(x, 2)) %>%
      # dplyr::filter(passed_ss == "TRUE" & diff_abn == "TRUE") -> out$ancombc2global_signif
      
      # out$ancombc2$res_pair %>%
      # mutate_if(is.numeric, function(x) round(x, 2))
      
      # out$ancombc2$res_trend %>%
      # mutate_if(is.numeric, function(x) round(x, 2)) %>%
      # filter(diff_abn == TRUE & passed_ss == TRUE) -> out$ancombc2trend_signif
      
    }
    
    if ("maaslin3"  %in% approach)
    {
      
      ########## ----- maaslin3 see https://github.com/biobakery/biobakery/wiki/MaAsLin3#44-level-contrasts
      
      # By default, the (up to) two metadata variables with the most significant associations will be plotted in the coefficient plot, and the rest will be plotted in the heatmap. Because predicting the output variable names can be tricky, it is recommended to first run maaslin3 without setting coef_plot_vars or heatmap_vars, look at the names of the variables in the summary plot, and then rerun with maaslin_plot_results_from_output after updating coef_plot_vars and heatmap_vars with the desired variables.
      
      ps_tmp %>% 
        transform_sample_counts(function(x) x/sum(x) * 100) %>% 
        filter_taxa(function(x){sum(x > 0) > maaslin3_prv_cut*nsamples(ps_tmp)}, prune = TRUE) %>% 
        physeq_glom_rename(speedyseq = TRUE, taxrank = ifelse(taxa_rank == "all", "Species", taxa_rank), rename_ASV = ifelse(taxa_rank == "all", "Species", taxa_rank)) %>%
        phyloseq_maaslin3(formula = maaslin3_formula,
                          fixed_effects =  maaslin3_fixed_effects, 
                          strata_effects = maaslin3_strata_effects,
                          correction = p_adjust,
                          max_significance = pvalue_cutoff,
                          output_dir = masslin3_output_dir,
                          min_prevalence = 0,
                          # coef_plot_vars = TRUE,
                          augment = TRUE,
                          plot_associations = TRUE, 
                          save_models = TRUE, 
                          plot_summary_plot = TRUE,
                          cores = 2, 
                          verbosity = "error") -> out$maaslin3 # random_effects = , group_effects = )
      
      # Save results in LEfSe format
      # maaslin_write_results_lefse_format(masslin3_output_dir, maaslin3$fit_data_abundance, maaslin3$fit_data_prevalence)
      # maaslin3::maaslin_plot_results_from_output()
      
      maaslin3:::preprocess_merged_results(rbind(out$maaslin3$fit_data_prevalence$results,
                                                 out$maaslin3$fit_data_abundance$results)) -> out$maaslin3$merged_res
      if(maaslin3_lefse_plot){
        run_make_coef_plot(merged_results_sig = out$maaslin3$merged_res %>% filter(qval_joint <= pvalue_cutoff), 
                           max_significance = pvalue_cutoff, 
                           class = comp_group) -> out$maaslin3$maaslin3_lefse_plot 
      }
      # maaslin3::maaslin_contrast_test()
      
    }
    
    if ("linda"  %in% approach)
    {
      ########## ----- microeco::linda
      
      ps_tmp %>% 
        # subset_taxa(Kingdom != unclassified_name) %>% 
        filter_tax_table(Kingdom != unclassified_name) %>% 
        transform_sample_counts(function(x) x/sum(x) * 100) %>% 
        filter_taxa(function(x){sum(x > 0) > linda_prv_cut*nsamples(ps_tmp)}, prune = TRUE) %>% 
        file2meco::phyloseq2meco(.) -> data
      
      
      t1 <- trans_diff$new(dataset = data, 
                           method = "linda", 
                           # formula = linda_formula,
                           remove_unknown = TRUE,
                           alpha = pvalue_cutoff, 
                           taxa_level = ifelse(taxa_rank == "all", "Species",taxa_rank),
                           group = linda_comp_group,
                           filter_thres = 0,
                           p_adjust_method = p_adjust)
      
      t1$res_diff -> out$linda$res_diff
      
      t1$res_diff  %<>%   subset(P.adj <= 0.05) # subset(Significance %in% c("*","**","***"))
      
      out$linda$diff_bar  <- t1$plot_diff_bar(keep_full_name = FALSE, 
                                              heatmap_cell =  "P.adj",
                                              heatmap_sig = "Significance",
                                              heatmap_x = "Factors",
                                              heatmap_y = "Taxa",
                                              heatmap_lab_fill = "P.adj")
      
      
      out$linda$all <-  clone(t1)
      
      # uses CLR?
      # using directly the linda package?
      
    }
    
    if ("classifier_rf"  %in% approach)
    {
      ########## ----- microeco::rf
      # or https://readingradio.github.io/J.nigra.Rmds/RF.GMW.Jnigra.html
      # error Error in serialize(data, node$con) : connection is not open
      
      ps_tmp %>% 
        # subset_taxa(Kingdom != unclassified_name) %>% 
        filter_tax_table(Kingdom != unclassified_name) %>% 
        transform_sample_counts(function(x) x/sum(x) * 100) %>% 
        filter_taxa(function(x){sum(x > 0) > rf_prv_cut*nsamples(ps_tmp)}, prune = TRUE) %>% 
        file2meco::phyloseq2meco(.) -> data
      
      # initialize: use "genotype" as response variable
      # x.predictors parameter is used to select the taxa; here we use all the taxa data in d1$taxa_abund
      t1 <- trans_classifier$new(dataset = data, y.response = comp_group, x.predictors = taxa_rank,  n.cores = 4)
      
      # generate train and test set
      t1$cal_split(prop.train = rf_prop_train)
      
      # Before training the model, we run the set_trainControl to invoke the trainControl function of caret package to generate the parameters used for training. 
      #Here we use the default parameters in trainControl function.
      t1$set_trainControl(method = "repeatedcv",
                          classProbs = TRUE,
                          savePredictions = TRUE)
      
      t1$cal_feature_sel(
        boruta.maxRuns = 300,
        boruta.pValue = 0.05,
        boruta.repetitions = 4)
      
      # use default parameter method = "rf"
      # require(doParallel)
      # library(caret)
      # library(randomForest)
      n <- parallel::detectCores()/2 # experiment!
      cl <- parallel::makeCluster(n)
      doParallel::registerDoParallel(cl)
      
      t1$cal_train(method = "rf", max.mtry = ref_train_max_mtry, ntree = ref_train_ntree)
      
      # t1$
      
      # t1$cal_caretList()
      
      # t1$cal_caretList_resamples()
      
      t1$cal_predict()
      
      out$classifier_rf$res_train <- t1$res_train
      
      # plot the confusionMatrix to check out the performance
      
      out$classifier_rf$res_confusion_stats <- t1$res_confusion_stats
      out$classifier_rf$res_confusion_fit <- t1$res_confusion_fit
      
      t1$plot_confusionMatrix()
      # t1$plot_confusion()
      
      t1$cal_ROC()
      #Using cal_ROC and plot_ROC can get the ROC (Receiver Operator Characteristic) curve.
      # out$Specificitysensitivity() <- t1$res_ROC$res_roc
      # out$RecallPrecision() <- t1$res_ROC$res_pr
      
      out$classifier_rf$plotROC  <- t1$plot_ROC(plot_method = FALSE)
      
      out$classifier_rf$plotROC2  <- t1$plot_ROC(plot_method = FALSE, plot_type = "PR")
      
      # default all groups
      #t1$plot_ROC(size = 0.5, alpha = 0.7)
      
      
      # default method in caret package without significance
      # t1$cal_feature_imp()
      
      # out$res_feature_imp <-  t1$res_feature_imp()
      
      # generate significance with rfPermute package
      t1$cal_feature_imp(rf_feature_sig = TRUE,group = comp_group,  num.rep = feature_imp_nrep, n.cores = 4) #, 
      
      out$classifier_rf$res_feature_imp <- t1$res_feature_imp
      
      out$classifier_rf$plot_feature_imp1 <- t1$plot_feature_imp(colour = "red", fill = "red", width = 0.6)
      
      # add_sig = TRUE: add significance label
      
      out$classifier_rf$plot_feature_imp2 <- t1$plot_feature_imp(coord_flip = TRUE, colour = "red", fill = "red", width = 0.6, add_sig = TRUE)
      
      out$classifier_rf$plot_feature_imp3 <- t1$plot_feature_imp(show_sig_group = TRUE, rf_sig_show = "MeanDecreaseGini", coord_flip = TRUE, width = 0.6, add_sig = TRUE, group_aggre = FALSE)
      
      out$classifier_rf$all <- clone(t1)
      
    } 
    
    if ("trans_diff_rf"  %in% approach)
    {
      ########## ----- microeco::rf
      # https://chiliubio.github.io/microeco_tutorial/model-based-class.html#trans_diff-class
      # The ‘rf’ method depends on the random forest(Beck and Foster 2014; Yatsunenko et al. 2012) and the non-parametric test. The current method implements random forest by bootstrapping like the operation in LEfSe and employs the significant features as input. MeanDecreaseGini is selected as the indicator value in the analysis.
      
      ps_tmp %>% 
        # subset_taxa(Kingdom != unclassified_name) %>% 
        filter_tax_table(Kingdom != unclassified_name) %>% 
        transform_sample_counts(function(x) x/sum(x) * 100) %>% 
        filter_taxa(function(x){sum(x > 0) > rf_prv_cut*nsamples(ps_tmp)}, prune = TRUE) %>% 
        file2meco::phyloseq2meco(.) -> data
      
      t1 <- trans_diff$new(dataset = data, 
                           method = "rf",
                           filter_thres = 0, #default 0; the abundance threshold,
                           alpha = pvalue_cutoff,
                           group = comp_group, 
                           taxa_level = taxa_rank, 
                           p_adjust_method = p_adjust,
                           plot_pal = palette,
                           nresam = rf_prop_train, #boots = 999, 
                           n.cores = 4)
      
      t1$res_diff %>% 
        dplyr::filter(MeanDecreaseGini >= trans_diff_rf_MeanDecreaseGini_cutoff) -> t1$res_diff
      
      
      if(is.null(plot_top_sig)){
        t1$res_diff %>%  nrow() -> plot_top_sig
      }
      
      
      out$trans_diff_rf$res_diff <- t1$res_diff
      # plot the MeanDecreaseGini bar
      # group_order is designed to sort the groups
      out$trans_diff_rf$g1 <- t1$plot_diff_bar(use_number = 1:plot_top_sig, color_values = palette)
      # plot the abundance using same taxa in g
      out$trans_diff_rf$g2 <- t1$plot_diff_abund(select_taxa = t1$plot_diff_bar_taxa, plot_type = "barerrorbar", add_sig = TRUE, errorbar_addpoint = FALSE, errorbar_color_black = TRUE, color_values = palette) # barerrorbar  ggboxplot
      # now the y axis in g1 and g2 is same, so we can merge them
      # remove g1 legend; remove g2 y axis text and ticks
      out$trans_diff_rf$g1 <- out$trans_diff_rf$g1 + theme(legend.position = "none")
      out$trans_diff_rf$g2 <- out$trans_diff_rf$g2 + theme(axis.text.y = element_blank(), 
                                                           axis.ticks.y = element_blank(), 
                                                           panel.border = element_blank()) + scale_y_continuous(trans='sqrt')
      out$trans_diff_rf$p <- out$trans_diff_rf$g1 %>% aplot::insert_right(out$trans_diff_rf$g2 )
      
      out$trans_diff_rf$all <- clone(t1)
      
    } 
    
    
  }))
  
  return(out)
}



#' @title Create Top Heatmap and Barplot Visualization for Microbiome Data
#' @author Florentin Constancias
#' @description
#' This function generates heatmaps and bar plots for a given phyloseq object
#' by displaying the most abundant taxa at various taxonomic levels (Phylum, Family, Genus, etc.).
#' 
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
#' 
#' @return A list of ggplot objects including the combined heatmap plot, legend, taxa information, bar plot, nested legend, and p plot.
#' 
#' @import ggplot2, phyloseq, ggpubr, microViz
#' 

phyloseq_top_heatmap_barplot <- function(
    ps_up,
    group_var = "Sample_Time", 
    tax_levels = c("Phylum", "Family", "Genus", "Species"), 
    ntax = 5, 
    ntax_species = 14,
    plot_heights = c(1.4, 1.5, 4),
    plot_x = "Subject",
    facet_by = c("Sample_Type", "Time"),
    group_by = c("Sample_Type", "Time"),
    facet_heat = "~ Sample_Type + Time",
    facet_formula = "Sample_Type ~ Time",
    rm_unclassified = TRUE,
    barplot_level = "Species",
    boxplot_main_group = "Class",
    per_sub_plot = TRUE) {
  
  require(ggpubr); require(tidyverse);require(speedyseq);require(ampvis2);require(microViz);require(rstatix);require(ggnested)
  source("https://raw.githubusercontent.com/fconstancias/DivComAnalyses/master/R/phyloseq_heatmap.R")
  
  
  suppressMessages({
    suppressWarnings({
      
      # Initialize output container
      out <- list()
      
      # Loop through the specified taxonomic levels
      for (tax in rank_names(ps_up)[rank_names(ps_up) %in% tax_levels]) {
        # print(tax)
        
        # Remove unclassified taxa (if specified) from all levels except "Phylum"
        if (rm_unclassified && tax != "Phylum") {
          ps_up <- ps_up %>% subset_taxa(Kingdom != "UNCLASSIFIED")
        }
        
        # Aggregate to taxonomic level and retain most abundant taxa
        out$most_ab_treat[[tax]] <- ps_up %>%
          tax_glom(taxrank = tax) %>%
          physeq_most_abundant(
            physeq = .,
            group_var = group_var,
            tax_level = tax,
            ntax = ifelse(tax == "Species", ntax_species, ntax)
          )
        
 
        # Transform to relative abundance (percentage), filter to most abundant taxa
        out$heat[[tax]] <- ps_up %>%
          transform_sample_counts(function(x) x / sum(x) * 100) %>% # conditional
          tax_glom(taxrank = tax) %>%
          filter_tax_table(get(tax) %in% out$most_ab_treat[[tax]]) %>%
          tax_mutate(Strain = NULL) %>%
          phyloseq_ampvis_heatmap(
            tax_aggregate = tax,
            physeq = .,
            tax_add = NULL,
            transform = FALSE,
            facet_by = facet_by,
            group_by = group_by,
            ntax = Inf
          )
        
        # Replace zero values in abundance data with NA
        out$heat[[tax]]$data <- out$heat[[tax]]$data %>%
          mutate(Abundance = na_if(Abundance, 0))
        
        # Apply color scaling and other theme adjustments to the heatmap plot
        out$heat[[tax]] <- out$heat[[tax]] +
          facet_grid(as.formula(facet_heat), scales = "free", space = "free") +
          scale_fill_viridis_c(
            breaks = c(0, 0.01, 1, 10, 50, 75, 100),
            labels = c(0, 0.01, 1, 10, 50, 75, 100),
            trans = scales::pseudo_log_trans(sigma = 0.001),
            na.value = 'transparent'
          ) +
          ylab(NULL) +
          theme(
            axis.title.x = element_blank(),
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank()
          )
      }
      
      # Arrange heatmaps for Phylum, Genus, and Species levels
      
      ggpubr::ggarrange(out$heat$Phylum + ylab(NULL) +
                          theme(axis.title.x=element_blank(),
                                axis.text.x=element_blank(),
                                axis.ticks.x=element_blank(),
                                legend.position = "none"),
                        out$heat$Genus +
                          theme(
                            strip.background = element_blank(),
                            strip.text.x = element_blank(),
                            legend.position = "none"
                          ) +
                          theme(axis.title.x=element_blank(),
                                axis.text.x=element_blank(),
                                axis.ticks.x=element_blank(),
                                legend.position = "none"),
                        out$heat$Species +
                          theme(
                            strip.background = element_blank(),
                            strip.text.x = element_blank(),
                            legend.position = "none"
                          ),
                        align = "v",
                        ncol = 1,
                        heights = plot_heights,
                        common.legend = FALSE) -> out$heat_all
      
      out$heat[[tax]] %>% 
        ggpubr::get_legend() %>% 
        ggpubr::as_ggplot() -> out$heat_legend
      
      
      # Extract legend for the heatmap
      legend <- ggpubr::get_legend(out$heat[[tax]]) %>%
        ggpubr::as_ggplot()
      
      # Create a bar plot with abundance proportions at specified taxonomic level
      bar_plot <- ps_up %>%
        subset_taxa(Class != "UNCLASSIFIED") %>%
        tax_glom(barplot_level) %>%
        transform_sample_counts(function(x) x / sum(x) * 100) %>% # conditional
        filter_tax_table(get(barplot_level) %in% out$most_ab_treat[[barplot_level]]) %>%
        microViz::comp_barplot(
          bar_width = 1,
          n_taxa = length(out$most_ab_treat[[barplot_level]]),
          tax_transform_for_plot = "identity",
          taxon_renamer = function(x) stringr::str_replace_all(x, "_", " "),
          label = plot_x,
          # x = plot_x,
          tax_level = barplot_level,
          merge_other = FALSE
        ) +
        ylab("Proportion - %") +
        theme_linedraw() +
        # theme(axis.ticks = element_blank(), axis.text.x = element_blank())
        theme(axis.ticks = element_blank())
      
      
      # Create the nested plot with a boxplot-style visualization
      p <- ggnested::ggnested(
        bar_plot$data,
        aes_string(main_group = boxplot_main_group, sub_group = "unique", x = plot_x, y = "Abundance")
      ) +
        scale_y_continuous(expand = c(0, 0)) +
        geom_bar(position = "stack", stat = "identity", color="grey5", linewidth = 0.1) +
        # theme_light() +  # Or any other preferred theme function like theme_minimal()
        # ylab("Proportion - %") +
        theme_linedraw() +
        theme(axis.ticks = element_blank(), axis.text.x = element_blank()) +
        facet_wrap(as.formula(facet_formula), scales = "free_x", drop = TRUE) +
        ylab("Proportion - %")
      
      # Extract legend for the nested plot
      out$nested_legend <- ggpubr::get_legend(p) %>%
        ggpubr::as_ggplot()
      
      # Remove legend from the main nested plot
      out$p <- p + theme(legend.position = "none")
      
      # Customize `bar_plot` with color scales from `p` and facet adjustments
      tax_pal <- p$data %>%
        distinct(unique, subgroup_colour) %>%
        pull(subgroup_colour, unique)
      
      out$bar_plot <- bar_plot +
        facet_wrap(as.formula(facet_formula), scales = "free_x", drop = TRUE) +
        # scale_color_manual(values = tax_pal) +
        scale_fill_manual(values = tax_pal) +
        theme(legend.position = "none")
      
      
      # Customize nested_plot per Subject
      
      if (isTRUE(per_sub_plot)){
        
        
        p$data %>% 
          select(-subgroup_colour, -group_subgroup, -group_colour) %>% 
          ggnested::ggnested(
            .,
            aes_string(main_group = "Class", sub_group = "unique", x = "Time", y = "Abundance")
          ) +
          scale_y_continuous(
            expand = expansion(add = c(0, 0.1))) + # axis starts exactly at 0
          geom_bar(position = "stack", stat = "identity", color="grey5", linewidth = 0.1) +
          theme_linedraw() +
          theme(axis.ticks = element_blank(), axis.text.x = element_blank()) +
          ylab("Proportion - %") + xlab(NULL) + facet_grid(rows = vars(Sample_Type), 
                                                           cols = vars(Subject ),
                                                           labeller = as_labeller(~ paste("", .))) + theme(legend.position = "none") -> out$ps_sub
        
      }
      
      # Return a list of plots and relevant components
      # out <- list(
      #   "heat_all" = heat_all,
      #   "legend" = legend,
      #   "tax" = out$most_ab_treat,
      #   "bar_plot" = bar_plot,
      #   "nested_legend" = nested_legend,
      #   "p" = p
      # )
      
      return(out)
      
    })
  })
}

# g <- ggplot_build(bar_plot)
# data.frame(colours = unique(g$data[[1]]["fill"]),
#            label = g$plot$scales$scales[[1]]$get_labels()) %>% 
#   pull(fill, label) -> tax_pal_mic
# 
# tax_pal_mic %>% 
#   length()
# 
# tax_pal_mic


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
#' 
#' @import tidyverse, phyloseq, vegan
#' 
#' @examples
#' # Example usage:
#' results <- compute_plot_beta(ps_up = physeq_obj, beta = beta_distances, color_group = "Group")
#' 
#' @export
compute_plot_beta <- function(ps_up = ps_up,
                              beta = beta,
                              m = "PCoA",
                              color_group = "Time",
                              shape_group = "Sample",
                              alpha = NULL,
                              col_pal = time_pal,
                              fill_pal = time_pal,
                              path_group = "interaction(Sample,Subject)",
                              facet_formula = "Sample ~ .",
                              axis1 = 1,
                              axis2 = 2,
                              seed = 123,
                              permanova_terms = c("Time", "cluster_Dtp2"),
                              metadata_dist_boxplot = NULL,
                              # terms_margins = "terms",
                              strata = "none",
                              perm = 999)
{
  ## ------------------------------------------------------------------------
  require(ggpubr); require(tidyverse);require(speedyseq);require(ampvis2);require(microViz);require(rstatix)
  
  source("https://raw.githubusercontent.com/fconstancias/DivComAnalyses/master/R/phyloseq_beta.R")
  
  # cat(paste0('\n##',"You are using tidyverse version ", packageVersion('tidyverse'),'\n\n'))
  # cat(paste0('\n##',"You are using phyloseq version ", packageVersion('phyloseq'),'\n\n'))
  
  
  
  suppressMessages(suppressWarnings({ 
    
    ####-------- 
    
    out <- NULL
    
    ####-------- compute beta div and arrange data
    
    # ps_up %>%
    #   # phyloseq_density_normalize(value_idx = "Dens") %>% 
    #   phyloseq_compute_bdiv() -> beta
    # 
    # ps_up %>% 
    #   microViz::dist_calc(., dist = "robust.aitchison") %>% 
    #   microViz::dist_get() -> beta$rAitchison
    # 
    # beta$sorensen = NULL
    # beta$bray = NULL
    
    # plot3D(beta$wjaccard, beta$rAitchison, beta$bjaccard)
    
    ####-------- Ordination
    
    ps_up %>% 
      phyloseq_plot_bdiv(dlist = beta, # list of distance computed from a phyloseq object
                         ps_rare = ., # phyloseq object
                         m = m, # PCoA or NMDS
                         seed = seed, # for reproducibility
                         axis1 = axis1, # axis to plot
                         axis2 = axis2) -> out$plot_list
    
    out$plot_list %>%
      phyloseq_plot_ordinations_facet(color_group = color_group,
                                      shape_group = shape_group,
                                      alpha = alpha)  + scale_color_manual(name = "", values = col_pal,
                                                                           na.value = "black") +
      scale_fill_manual(name = "", values = fill_pal,
                        na.value = "black") + theme_linedraw()   + theme(legend.position = "right") + facet_null()  + facet_wrap(distance ~., scales = "free", nrow = 3) -> out$pcoas
    
    
    # out$PCOA %>% 
    #   ggpubr::get_legend(.) %>% 
    #   ggpubr::as_ggplot(.) ->  out$PCOA_leg
    # 
    # out$PCOA + theme(legend.position = "none") ->    out$PCOA
    # pcoas
    
    out$plot_list %>%
      phyloseq_ordinations_expl_var() -> out$expl_var
    
    ####-------- Ordination detailed
    
    out$plot_list$rAitchison$layers[[1]] = NULL;  out$plot_list$rAitchison$layers[[1]] = NULL
    out$plot_list$rAitchison$layers[[2]] = NULL;  out$plot_list$rAitchison$layers[[1]] = NULL
    
    # plots_hall_humans$aichinson$layers[[1]] = NULL;plots_hall_humans$aichinson$layers[[1]] = NULL
    # plots_hall_humans$aichinson$layers[[2]] = NULL;plots_hall_humans$aichinson$layers[[2]] = NULL
    
    out$plot_list$rAitchison + geom_point(size = 3,
                                          aes_string(colour = color_group ,
                                                     shape = shape_group,
                                                     alpha = alpha)) +
      geom_path(data =  out$plot_list$rAitchison$data %>%
                  arrange(Subject) ,
                # aes(colour = Treatment, group = interaction(Model, Model2, Antibiotic, Treatment, Fermentation, Reactor,Antibiotic_mg.L)),
                aes_string(group = path_group),
                
                arrow = arrow(
                  angle = 30, length = unit(0.15, "inches"),
                  ends = "last", type = "open"
                ), linetype = "longdash", size = 0.1) +
      theme_light() +
      scale_color_manual(name = "", values = col_pal,
                         na.value = "black") +
      scale_fill_manual(name = "", values = fill_pal,
                        na.value = "black") +
      # scale_shape_manual(name = "" ,values = c(15,16,18,19), na.value =  17) +
      theme(legend.position = "right") +
      facet_grid(as.formula(facet_formula), scales = "free_y", space = "fixed", switch = "y") +   theme_linedraw() + theme(strip.placement = "outside")  -> out$PCOA
    
    
    out$PCOA %>% 
      ggpubr::get_legend(.) %>% 
      ggpubr::as_ggplot(.) ->  out$PCOA_leg
    
    out$PCOA + theme(legend.position = "none") ->    out$PCOA
    
    
    
    
    ####-------- envfit
    
    out$PCOA + 
      scale_fill_manual(values = c("transparent")) + 
      scale_color_manual(values = c(rep("transparent", length(col_pal)))) + 
      theme(panel.border = element_blank(), panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) -> empty_plot_tmp
    
    # ps_up %>%
    #   phyloseq_add_taxa_vector_fix(phyloseq = ., perm = perm,
    #                                dist = beta$rAitchison,
    #                                # tax_rank_plot = "Genus",
    #                                taxrank_glom = "Species",
    #                                figure_ord = empty_plot_tmp,
    #                                adj_method = "fdr",
    #                                fact = 0.8, pval_cutoff = 0.05,
    #                                top_r = 10) -> out$envfit
    
    # ps_up %>%
    #   phyloseq_add_metadata_vector(dist = beta$rAitchison,
    #                                phyloseq = .,
    #                                figure_ord = empty_plot_tmp,
    #                                m = "PCoA",
    #                                pval_cutoff = 0.05,
    #                                top_r = 12,
    #                                metadata_sel = c("mean_plaque", "mean_bleeding", "age"),
    #                                fact = 0.5,
    #                                seed = 123,
    #                                perm = 999,
    #                                norm_method = "center_scale",
    #                                color = "green",
    #                                linetype = "dashed",
    #                                na.rm = TRUE)
    
    
    ####-------- permnova
    
    tmp_out1 = NULL; tmp_out2 = NULL
    for (terms_margins in c("terms", "margin")) {
      
      form = paste0(permanova_terms, collapse=" * ")
      
      
      ps_up %>% 
        lapply(
          beta,
          FUN = phyloseq_adonis2,
          physeq = .,
          formula = form,
          nrep = perm,
          terms_margins = terms_margins
          # strata = "Subject"
        )  %>%
        bind_rows(.id = "Distance") %>% 
        mutate("Group" = (as.vector(form)),
               "terms_margin" = terms_margins) %>%   bind_rows(.,tmp_out1)  -> tmp_out1#perm1
      
      form = paste0(permanova_terms, collapse=" + ")
      ## TODO: Add all combination of order for the + one 
      
      
      ps_up %>% 
        lapply(
          beta,
          FUN = phyloseq_adonis2,
          physeq = .,
          formula = form,
          nrep = perm,
          terms_margins = terms_margins,
          strata = strata
        )  %>%
        bind_rows(.id = "Distance") %>% 
        mutate("Group" = (as.vector(form)),
               "terms_margin" = terms_margins) %>%   bind_rows(.,tmp_out2) -> tmp_out2
      
    }
    
    bind_rows(tmp_out1, tmp_out2) -> out$perm
    
    ####-------- PW permnova
    
    # tmp_pw_perm = NULL; tmp_tw_perm = NULL
    
    for (compare_header in permanova_terms) {
      
      # print(compare_header)
      
      ps_up %>% 
        lapply(
          beta,
          FUN = physeq_pairwise_permanovas_adonis2,
          physeq = .,
          compare_header = compare_header,
          n_perm = perm,
          strata = strata,
          terms_margins = "terms"
        )  %>%
        bind_rows(.id = "Distance") %>% 
        mutate("Group" = (as.vector(compare_header))) %>% 
        bind_rows(.,out$pw_perm) -> out$pw_perm #out[[paste0("pw_", compare_header)]] 
      
      
      # phyloseq_generate_pcoa_per_variables
      
      # adonis_OmegaSq
      # physeq_betadisper
      
      lapply(
        beta,
        FUN = phyloseq_TW,
        physeq = ps_up,
        variable = compare_header,
        nrep = perm) %>% #,
        # strata = ifelse(strata == "none", NULL, strata))  %>%
        bind_rows(.id = "Distance") %>%
        # ddd %>% 
        mutate("Group" = (as.vector(compare_header)))  %>% 
        bind_rows(., out$tw_perm) -> out$tw_perm
      
    }
    
    
    ####-------- permnova
    for (compare_header in permanova_terms) {
      
      phyloseq_generate_pcoa_per_variables(tmp = ps_up,
                                           group = compare_header,
                                           m = "PCoA",
                                           dist = beta,
                                           color_group = color_group,
                                           shape_group = shape_group,
                                           alpha = alpha,
                                           col_pal = col_pal,
                                           fill_pal =  fill_pal) -> out[[compare_header]]
      
      
    }
    
    
    ####-------- distance boxplot
    
    if( !is.null(metadata_dist_boxplot)){
      
      lapply(
        beta,
        FUN = phyloseq_distance_boxplot,
        p = ps_up,
        d = "Sample") -> dist_bx
      
      ps_up %>%
        sample_data()%>%
        data.frame() %>%
        dplyr::select(any_of(metadata_dist_boxplot)) %>%
        rownames_to_column("Var") %>% 
        drop_na() -> meta_sel
      
      meta_sel_2 <- meta_sel
      
      colnames(meta_sel) <- paste0(  colnames(meta_sel) , "_1")
      colnames(meta_sel_2) <- paste0(  colnames(meta_sel_2) , "_2")
      
      dist_df = NULL
      for (d in names(dist_bx))
      {
        dist_bx[[d]]$matrix %>% 
          mutate("Distance" = d) %>% 
          bind_rows(.,dist_df) -> dist_df
      }
      
      dist_df %>% 
        left_join(meta_sel,
                  by = c("Var1" = "Var_1")) %>% 
        left_join(meta_sel_2,
                  by =  c("Var2" = "Var_2")) -> dist_df
      
      dist_df %>% #dplyr::filter((Time_1 == "TP1" & Subject_1 == Subject_2 & Sample_1 == Sample_2 ))
        # mutate( 
        #   value = case_when(Time_1 == Time_2 & Subject_1 == Subject_2 & Sample_1 == Sample_2 ~  0,
        #                     .default = value))  %>% 
        dplyr::filter(Sample_1 == Sample_2,
                      # Time_1 != Time_2,
                      Time_1 == "TP1",# or add columns at TP1 distance = 0 ?
                      # Time_2 != "TP1",
                      Subject_1 == Subject_2) %>% # arrange(value)
        arrange(Distance, Subject_2, Time_2) %>% 
        # filter(dist_1 == "wjaccard") %>%
        ggplot(data = ., aes_string(x="Time_2", y="value")) +
        geom_boxplot(outlier.colour = NA, alpha=0.7, aes_string(fill = "Time_2")) +
        # ggbeeswarm::geom_beeswarm(size=1, alpha=0.2,
        #                           position=pd) +
        geom_jitter(size=1, position = pd, aes_string(color = "Subject_1")) +
        # aes_string(shape = "cluster_Dtp2_1")) + 
        geom_line(aes_string(group="Subject_1"), position = pd, linetype = "dashed", color = "grey50", linewidth = 0.08) +
        
        facet_grid(as.formula(paste0("Distance ~ Sample_1")), scales = "free_y", space = "fixed", switch = "y") +
        scale_color_manual(name = "", values = sub_pal,
                           na.value = "black") +
        scale_fill_manual(name = "", values = time_pal,
                          na.value = "black") +
        theme_light() + ylab("Distance to Baseline") + xlab(NULL) + theme(
          axis.text.x = element_blank()) +  theme_linedraw() + theme(strip.placement = "outside") -> dist_box
      
      dist_box %>% 
        ggpubr::get_legend(.) %>% 
        ggpubr::as_ggplot(.) -> dist_box_leg
      
      dist_box + theme(legend.position = "none") -> dist_box
      
      
      out$dist_df <- dist_df
      out$dist_box <- dist_box
      out$dist_box_leg <- dist_box_leg
      
    }
    return(out)
  }))
}

phyloseq_explore_beta <- function(ps_up = ps_up,
                                  beta = beta,
                                  m = "PCoA",
                                  color_group = "Time",
                                  shape_group = "Sample",
                                  alpha = NULL,
                                  col_pal = time_pal,
                                  fill_pal = time_pal,
                                  path_group = "interaction(Sample,Subject)",
                                  facet_formula = "Sample ~ .",
                                  axis1 = 1,
                                  axis2 = 2,
                                  seed = 123,
                                  permanova_terms = c("Time", "cluster_Dtp2"),
                                  metadata_dist_boxplot = NULL,
                                  strata = "none",
                                  run_phyloseq_TW = TRUE,
                                  perm = 999)
{
  suppressMessages({
    suppressWarnings({
      
      ## Load necessary packages
      # require(tidyverse); require(phyloseq); require(vegan)
      # cat(paste0('\n##',"Using tidyverse version ", packageVersion('tidyverse'),'\n'))
      # cat(paste0('\n##',"Using phyloseq version ", packageVersion('phyloseq'),'\n'))
      # 
      out <- list() # Initialize output list
      
      #### 1. Compute beta diversity and organize data
      
      # Compute beta diversity using multiple distance metrics
      # Uncomment and modify as needed to calculate additional beta diversity metrics
      # ps_up %>%
      #   microViz::dist_calc(dist = "robust.aitchison") %>% 
      #   microViz::dist_get() -> beta$rAitchison
      
      #### 2. Ordination
      ps_up %>% 
        phyloseq_plot_bdiv(dlist = beta, ps_rare = ., m = m, seed = seed, axis1 = axis1, axis2 = axis2) -> out$plot_list
      
      # Generate PCoA ordination plots with specified color and shape groups
      out$plot_list %>%
        phyloseq_plot_ordinations_facet(color_group = color_group, shape_group = shape_group, alpha = alpha) +
        scale_color_manual(name = "", values = col_pal, na.value = "black") +
        scale_fill_manual(name = "", values = fill_pal, na.value = "black") + 
        theme_linedraw() + theme(legend.position = "right") +
        facet_wrap(distance ~ ., scales = "free", nrow = 3) -> out$pcoas
      
      # Extract explained variance from ordination plot
      out$plot_list %>%
        phyloseq_ordinations_expl_var() -> out$expl_var
      
      #### 3. Detailed Ordination for Robust Aitchison Distance
      out$plot_list$rAitchison + geom_point(size = 3, aes_string(colour = color_group, shape = shape_group, alpha = alpha)) +
        geom_path(data = out$plot_list$rAitchison$data %>% arrange(Subject), aes_string(group = path_group),
                  arrow = arrow(angle = 30, length = unit(0.15, "inches"), ends = "last", type = "open"), linetype = "longdash", size = 0.1) +
        theme_light() + scale_color_manual(name = "", values = col_pal, na.value = "black") +
        scale_fill_manual(name = "", values = fill_pal, na.value = "black") +
        facet_grid(as.formula(facet_formula), scales = "free_y", space = "fixed", switch = "y") + theme_linedraw() + theme(strip.placement = "outside") -> out$PCOA
      
      # Extract and hide legend for the ordination plot
      out$PCOA %>% ggpubr::get_legend() %>% ggpubr::as_ggplot() -> out$PCOA_leg
      out$PCOA + theme(legend.position = "none") -> out$PCOA
      
      #### 4. Environmental Fit (optional)
      out$PCOA + scale_fill_manual(values = c("transparent")) +
        scale_color_manual(values = rep("transparent", length(col_pal))) + theme(panel.border = element_blank(),
                                                                                 panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) -> empty_plot_tmp
      
      #### 5. PERMANOVA Analysis
      tmp_out1 = NULL; tmp_out2 = NULL
      for (terms_margins in c("terms", "margin")) {
        form = paste0(permanova_terms, collapse = " * ")
        ps_up %>% lapply(beta, FUN = phyloseq_adonis2, physeq = ., formula = form, nrep = perm, terms_margins = terms_margins) %>%
          bind_rows(.id = "Distance") %>% mutate("Group" = as.vector(form), "terms_margin" = terms_margins) %>%
          bind_rows(., tmp_out1) -> tmp_out1
        form = paste0(permanova_terms, collapse = " + ")
        ps_up %>% lapply(beta, FUN = phyloseq_adonis2, physeq = ., formula = form, nrep = perm, terms_margins = terms_margins, strata = strata) %>%
          bind_rows(.id = "Distance") %>% mutate("Group" = as.vector(form), "terms_margin" = terms_margins) %>%
          bind_rows(., tmp_out2) -> tmp_out2
      }
      bind_rows(tmp_out1, tmp_out2) -> out$perm
      
      #### 6. Pairwise PERMANOVA
      for (compare_header in permanova_terms) {
        ps_up %>% lapply(beta, FUN = physeq_pairwise_permanovas_adonis2, physeq = ., compare_header = compare_header, n_perm = perm, strata = strata, terms_margins = "terms") %>%
          bind_rows(.id = "Distance") %>% mutate("Group" = as.vector(compare_header)) %>% bind_rows(., out$pw_perm) -> out$pw_perm
        if(run_phyloseq_TW){
        lapply(beta, FUN = phyloseq_TW, physeq = ps_up, variable = compare_header, nrep = perm) %>%
          bind_rows(.id = "Distance") %>% mutate("Group" = as.vector(compare_header)) %>% bind_rows(., out$tw_perm) -> out$tw_perm
        }
      }
      
      #### 7. Generate PCoA plots for each variable in permanova_terms
      for (compare_header in permanova_terms) {
        phyloseq_generate_pcoa_per_variables(tmp = ps_up, group = compare_header, m = "PCoA", dist = beta,
                                             color_group = color_group, shape_group = shape_group, alpha = alpha,
                                             col_pal = col_pal, fill_pal = fill_pal) -> out[[compare_header]]
      }
      
      ####-------- distance boxplot
      
      if( !is.null(metadata_dist_boxplot)){
        
        lapply(
          beta,
          FUN = phyloseq_distance_boxplot,
          p = ps_up,
          d = "Sample") -> dist_bx
        
        ps_up %>%
          sample_data()%>%
          data.frame() %>%
          dplyr::select(any_of(metadata_dist_boxplot)) %>%
          rownames_to_column("Var") %>% 
          drop_na() -> meta_sel
        
        meta_sel_2 <- meta_sel
        
        colnames(meta_sel) <- paste0(  colnames(meta_sel) , "_1")
        colnames(meta_sel_2) <- paste0(  colnames(meta_sel_2) , "_2")
        
        dist_df = NULL
        for (d in names(dist_bx))
        {
          dist_bx[[d]]$matrix %>% 
            mutate("Distance" = d) %>% 
            bind_rows(.,dist_df) -> dist_df
        }
        
        dist_df %>% 
          left_join(meta_sel,
                    by = c("Var1" = "Var_1")) %>% 
          left_join(meta_sel_2,
                    by =  c("Var2" = "Var_2")) -> dist_df
        
        dist_df %>% #dplyr::filter((Time_1 == "TP1" & Subject_1 == Subject_2 & Sample_1 == Sample_2 ))
          # mutate( 
          #   value = case_when(Time_1 == Time_2 & Subject_1 == Subject_2 & Sample_1 == Sample_2 ~  0,
          #                     .default = value))  %>% 
          dplyr::filter(Sample_1 == Sample_2,
                        # Time_1 != Time_2,
                        Time_1 == "TP1",# or add columns at TP1 distance = 0 ?
                        # Time_2 != "TP1",
                        Subject_1 == Subject_2) %>% # arrange(value)
          arrange(Distance, Subject_2, Time_2) %>% 
          # filter(dist_1 == "wjaccard") %>%
          ggplot(data = ., aes_string(x="Time_2", y="value")) +
          geom_boxplot(outlier.colour = NA, alpha=0.7, aes_string(fill = "Time_2")) +
          # ggbeeswarm::geom_beeswarm(size=1, alpha=0.2,
          #                           position=pd) +
          geom_jitter(size=1, position = pd, aes_string(color = "Subject_1")) +
          # aes_string(shape = "cluster_Dtp2_1")) + 
          geom_line(aes_string(group="Subject_1"), position = pd, linetype = "dashed", color = "grey50", linewidth = 0.08) +
          
          facet_grid(as.formula(paste0("Distance ~ Sample_1")), scales = "free_y", space = "fixed", switch = "y") +
          scale_color_manual(name = "", values = sub_pal,
                             na.value = "black") +
          scale_fill_manual(name = "", values = time_pal,
                            na.value = "black") +
          theme_light() + ylab("Distance to Baseline") + xlab(NULL) + theme(
            axis.text.x = element_blank()) +  theme_linedraw() + theme(strip.placement = "outside") -> dist_box
        
        dist_box %>% 
          ggpubr::get_legend(.) %>% 
          ggpubr::as_ggplot(.) -> dist_box_leg
        
        dist_box + theme(legend.position = "none") -> dist_box
        
        
        out$dist_df <- dist_df
        out$dist_box <- dist_box
        out$dist_box_leg <- dist_box_leg
        
      }
      
      return(out)
    })
  })
}



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
#' 
#' 
phyloseq_explore_alpha <- function(alpha_long_df,
                                   alpha_measures = c("observed", "diversity_shannon", "diversity_inverse_simpson", "diversity_coverage", "evenness_pielou"),
                                   pd = position_dodge(0.3),
                                   x = "Time",
                                   y = "value",
                                   color_point = "Subject",
                                   shape_point = "Sample",
                                   group_point = "interaction(Sample,Subject)",
                                   group_boxplot = "interaction(Sample,Time)",
                                   group_line = "Subject",
                                   facet_formula = "alphadiversiy ~ Sample + cluster_Dtp2",
                                   col_pal = sub_pal,
                                   fill_pal = sub_pal,
                                   stat_formula = "value ~ Time",
                                   anova_test_formula = "value ~ Time*Sample + Error(Subject/(Time*Sample))",
                                   group_by_stats = c("alphadiversiy", "Sample", "cluster_Dtp2"),
                                   padjust_method = "fdr",
                                   stat_paired = FALSE,
                                   ref_group_stat = NULL) {
  
  
  source("https://raw.githubusercontent.com/fconstancias/DivComAnalyses/master/R/phyloseq_alpha.R")
  require(ggpubr); require(tidyverse);require(speedyseq);require(ampvis2);require(microViz);require(rstatix)
  
  
  # Suppress messages and warnings for all subsequent steps
  suppressMessages(suppressWarnings({
    
    ####-------- Compute Summary Statistics
    alpha_summary <- alpha_long_df %>%
      group_by(alphadiversiy, Time, Sample) %>%
      rstatix::get_summary_stats("value")
    
    ####-------- Plotting
    alpha_plot <- alpha_long_df %>%
      ggplot(aes_string(x = x, y = y)) +
      geom_boxplot(aes_string(group = group_boxplot), outlier.shape = NA, alpha = 0.4) +
      geom_point(aes_string(color = color_point, shape = shape_point, group = group_point), position = pd) +
      geom_line(aes_string(group = group_line), position = pd, linetype = "dashed", color = "grey50", linewidth = 0.08) +
      facet_grid(as.formula(facet_formula), scales = "free", switch = "y") +
      theme_linedraw() + theme(strip.placement = "outside") +
      xlab(NULL) + ylab(NULL) +
      scale_color_manual(values = col_pal) + scale_fill_manual(values = fill_pal)
    
    # Separate legend for reuse or customization
    alpha_legend <- ggpubr::get_legend(alpha_plot) %>% ggpubr::as_ggplot()
    alpha_plot <- alpha_plot + theme(legend.position = "none")
    
    ####-------- Statistical Tests
    # Pairwise Wilcoxon Test
    wilcox_test <- ggpubr::compare_means(
      formula = as.formula(stat_formula),
      group.by = group_by_stats,
      paired = stat_paired,
      ref.group = ref_group_stat,
      data = alpha_long_df,
      method = "wilcox.test",
      p.adjust.method = padjust_method
    ) %>%
      select(-.y., -p.format, -p.signif) %>%
      arrange(p) %>%
      rstatix::add_significance(p.col = "p.adj")
    
    # Pairwise T-Test
    t_test <- ggpubr::compare_means(
      formula = as.formula(stat_formula),
      group.by = group_by_stats,
      paired = stat_paired,
      ref.group = ref_group_stat,
      data = alpha_long_df,
      method = "t.test",
      p.adjust.method = padjust_method
    ) %>%
      select(-.y., -p.format, -p.signif) %>%
      arrange(p) %>%
      rstatix::add_significance(p.col = "p.adj")
    
    # ANOVA
    anova_result <- alpha_long_df %>%
      group_by(alphadiversiy, Time, Subject) %>%
      filter(n() > 1) %>%
      ungroup() %>%
      group_by(alphadiversiy, Time, Sample) %>%
      filter(n() > 39) %>%
      ungroup() %>%
      group_by(alphadiversiy) %>%
      rstatix::anova_test(formula = as.formula(anova_test_formula))
    
    anova_table <- rstatix::get_anova_table(anova_result)
    
    # Kruskal-Wallis Test
    kruskal_test <- ggpubr::compare_means(
      formula = as.formula(stat_formula),
      group.by = group_by_stats,
      paired = stat_paired,
      ref.group = NULL,
      data = alpha_long_df,
      method = "kruskal.test",
      p.adjust.method = padjust_method
    ) %>%
      select(-.y., -p.format, -p.signif) %>%
      arrange(p) %>%
      rstatix::add_significance(p.col = "p.adj")
    
    ####-------- Output
    list(
      alpha_summary = alpha_summary,
      alpha_plot = alpha_plot,
      alpha_legend = alpha_legend,
      anova_result = anova_result,
      anova_table = anova_table,
      kruskal_test = kruskal_test,
      wilcox_test = wilcox_test,
      t_test = t_test
    )
    
  }))
}

compute_plot_alpha <- function(alpha_long_df = alpha_long_df,
                               alpha_measures = c("observed", "diversity_shannon", "diversity_inverse_simpson", "diversity_coverage", "evenness_pielou"),
                               pd = position_dodge(0.3),
                               x = "Time",
                               y = "value",
                               color_point = "Subject",
                               shape_point = "Sample",
                               group_point = "interaction(Sample,Subject)",
                               group_boxplot = "interaction(Sample,Time)",
                               group_line = "Subject",
                               facet_formula = "alphadiversiy ~ Sample + cluster_Dtp2",
                               col_pal = sub_pal,
                               fill_pal = sub_pal,
                               stat_formula = "value ~ Time",
                               anova_test_formula = "value ~ Time*Sample + Error(Subject/(Time*Sample))",
                               group_by_stats = c("alphadiversiy", "Sample", "cluster_Dtp2"),
                               padjust_method = "fdr",
                               stat_paired = FALSE,
                               ref_group_stat = NULL){
  out = NULL
  
  ####-------- compute alpha div and arrange data
  
  
  # ps_up %>%
  #   phyloseq_alphas() -> alpha
  
  # alpha %>%
  #   pivot_longer(cols = all_of(alpha_measures), values_to = "value", names_to = 'alphadiversiy', values_drop_na  = TRUE) %>%
  #   mutate(alphadiversiy = fct_relevel(alphadiversiy, alpha_measures)) -> alpha_long_df
  
  alpha_long_df %>%
    group_by(alphadiversiy, Time, Sample) %>%
    rstatix::get_summary_stats("value") -> alpha_sum
  
  ####-------- plot
  
  alpha_long_df %>%
    arrange(time) %>%
    # filter(!is.na(cluster_Dtp2)) %>%
    ggplot(aes_string(x = x,
                      y = y)) +
    geom_boxplot(aes_string(group= group_boxplot, fill = NULL, color = NULL), outlier.shape = NA, alpha = 0.4) +
    geom_point(aes_string(color = color_point, group=group_point,
                          shape = shape_point),  position = pd) +
    geom_line(aes_string(group=group_line), position = pd, linetype = "dashed", color = "grey50", linewidth = 0.08) +
    facet_grid(as.formula(facet_formula), scales = "free", switch = "y") +
    theme_linedraw() + theme(strip.placement = "outside") + xlab(NULL) + ylab(NULL) + scale_color_manual(values = col_pal) + scale_fill_manual(values = fill_pal) -> alpha_plot
  
  alpha_plot %>%
    ggpubr::get_legend(.) %>%
    ggpubr::as_ggplot(.) -> alpha_plot_leg
  
  alpha_plot + theme(legend.position = "none") -> alpha_plot
  
  
  ####-------- stats
  
  # would be cool to be able to write text -> code for group_by mainly
  
  # alpha_long_df %>%
  #   # filter(!is.na(cluster_Dtp2)) %>%
  #   # group_by(alphadiversiy, Sample, cluster_Dtp2) %>%
  #     # group_by(eval(parse(text="alphadiversiy, Sample"))) %>%
  #
  #   rstatix::pairwise_wilcox_test(formula = as.formula(stat_formula),
  #                                 p.adjust.method = "none", paired = stat_paired, ref.group = ref_group_stat,detailed = TRUE) %>%
  #   group_by(alphadiversiy, Sample, cluster_Dtp2) %>%
  #   rstatix::adjust_pvalue(method = padjust_method) %>%
  #   select(-'.y.', -method, -alternative) %>%
  #   arrange(p) %>%
  #   mutate_if(is.numeric, ~round(., 5)) -> alpha_pw_w_stats
  #
  # # !!! rlang::syms(crit1)
  #
  #
  
  ggpubr::compare_means(formula = as.formula(stat_formula),
                        group.by = group_by_stats,
                        paired = stat_paired,
                        ref.group = ref_group_stat,
                        data =alpha_long_df,
                        method = "wilcox.test",
                        p.adjust.method = padjust_method) %>%
    select(-.y., -p.format, -p.signif) %>%
    arrange(p) %>%
    # mutate(signif = ifelse(p.adj <= 0.05, 'SIGN', 'NS')) %>%
    rstatix::add_significance(p.col = "p.adj") -> wilcox.test_stat
  
  ggpubr::compare_means(formula = as.formula(stat_formula),
                        group.by = group_by_stats,
                        paired = stat_paired,
                        ref.group = ref_group_stat,
                        data =alpha_long_df,
                        method = "t.test",
                        p.adjust.method = padjust_method) %>%
    select(-.y., -p.format, -p.signif) %>%
    arrange(p) %>%
    # mutate(signif = ifelse(p.adj <= 0.05, 'SIGN', 'NS')) %>%
    rstatix::add_significance(p.col = "p.adj") -> t.test_stat
  
  ggpubr::compare_means(formula = as.formula(stat_formula),
                        group.by = group_by_stats,
                        paired = stat_paired,
                        ref.group = NULL,
                        data =alpha_long_df,
                        method = "anova",
                        p.adjust.method = padjust_method) %>%
    select(-.y., -p.format, -p.signif) %>%
    arrange(p) %>%
    # mutate(signif = ifelse(p.adj <= 0.05, 'SIGN', 'NS')) %>%
    rstatix::add_significance(p.col = "p.adj") -> anova_stat
  
  ggpubr::compare_means(formula = as.formula(stat_formula),
                        group.by = group_by_stats,
                        paired = stat_paired,
                        ref.group = NULL,
                        data = alpha_long_df,
                        method = "kruskal.test",
                        p.adjust.method = padjust_method) %>%
    select(-.y., -p.format, -p.signif) %>%
    arrange(p) %>%
    # mutate(signif = ifelse(p.adj <= 0.05, 'SIGN', 'NS')) %>%
    rstatix::add_significance(p.col = "p.adj") -> kruskal.test_stat
  
  
  alpha_long_df %>%
    group_by(alphadiversiy, Time, Subject) %>%
    add_count() %>%
    filter(n > 1) %>%
    ungroup() %>%
    group_by(alphadiversiy, Time, Sample) %>%
    add_count(name = "n2") %>%
    filter(n2 > 39) %>%
    ungroup() %>%
    group_by(alphadiversiy) %>%
    rstatix::anova_test(formula = as.formula(anova_test_formula)) -> anova
  
  anova %>%
    rstatix::get_anova_table(.) -> anova_table
  
  
  out = list("alpha" = alpha,
             "alpha_long_df" = alpha_long_df,
             "alpha_sum" = alpha_sum,
             "alpha_plot" = alpha_plot,
             "alpha_plot_legend" = alpha_plot_leg,
             "anova" = anova,
             "anova_table" = anova_table,
             "kruskal.test_stat" = kruskal.test_stat,
             "anova_stat" = anova_stat,
             "t.test_stat" = t.test_stat,
             "wilcox.test_stat" = wilcox.test_stat)
  
  return(out)
  
}
