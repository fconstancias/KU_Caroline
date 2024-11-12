
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
  require(tidyverse); require(phyloseq); require(vegan)
  cat(paste0('\n##',"You are using tidyverse version ", packageVersion('tidyverse'),'\n\n'))
  cat(paste0('\n##',"You are using phyloseq version ", packageVersion('phyloseq'),'\n\n'))
  
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
  
}

compute_plot_beta_rev <- function(ps_up = ps_up,
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
                              perm = 999)
{
  ## Load necessary packages
  require(tidyverse); require(phyloseq); require(vegan)
  cat(paste0('\n##',"Using tidyverse version ", packageVersion('tidyverse'),'\n'))
  cat(paste0('\n##',"Using phyloseq version ", packageVersion('phyloseq'),'\n'))
  
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
    lapply(beta, FUN = phyloseq_TW, physeq = ps_up, variable = compare_header, nrep = perm) %>%
      bind_rows(.id = "Distance") %>% mutate("Group" = as.vector(compare_header)) %>% bind_rows(., out$tw_perm) -> out$tw_perm
  }
  
  #### 7. Generate PCoA plots for each variable in permanova_terms
  for (compare_header in permanova_terms) {
    phyloseq_generate_pcoa_per_variables(tmp = ps_up, group = compare_header, m = "PCoA", dist = beta,
                                         color_group = color_group, shape_group = shape_group, alpha = alpha,
                                         col_pal = col_pal, fill_pal = fill_pal) -> out[[compare_header]]
  }
  
  #### 8. Distance Boxplot (optional)
  if (!is.null(metadata_dist_boxplot)) {
    lapply(beta, FUN = phyloseq_distance_boxplot, p = ps_up, d = "Sample") -> dist_bx
    ps_up %>% sample_data() %>% data.frame() %>% dplyr::select(any_of(metadata_dist_boxplot)) %>%
      rownames_to_column("Var") %>% drop_na() -> meta_sel
    meta_sel_2 <- meta_sel
    colnames(meta_sel) <- paste0(colnames(meta_sel), "_1")
    colnames(meta_sel_2) <- paste0(colnames(meta_sel_2), "_2")
    
    dist_df <- NULL
    for (d in names(dist_bx)) {
      dist_bx[[d]]$matrix %>% mutate("Distance" = d) %>% bind_rows(., dist_df) -> dist_df
    }
    
    dist_df %>% left_join(meta_sel, by = c("Var1" = "Var_1")) %>%
      left_join(meta_sel_2, by = c("Var2" = "Var_2")) -> dist_df
    
    dist_df %>% filter(Sample_1 == Sample_2, Time_1 == "TP1", Subject_1 == Subject_2) %>%
      arrange(Distance, Subject_2, Time_2) %>%
      ggplot(data = ., aes_string(x = "Time_2", y = "value", color = "Sample_2")) + geom_boxplot() + theme_minimal() -> out$dist_boxplot
  }
  
  return(out)
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
compute_plot_alpha_rev <- function(alpha_long_df,
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
