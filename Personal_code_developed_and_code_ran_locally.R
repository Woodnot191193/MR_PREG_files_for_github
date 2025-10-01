# Load environment variables from config file
readRenviron("config.env")

# Retrieve project directory from environment
projectdir <- Sys.getenv("MRPregdir")
print(projectdir)

# Use `projectdir` to navigate folders without hardcoding paths
setwd(paste0(projectdir, "/MR_PREG_files_for_github/"))
getwd()

# Return to main project directory
setwd(projectdir)
getwd()

exposure_dat <- format_data(
  exposure_list[["F-ALL"]],
  type = "exposure",
  snp_col = "RSID",
  beta_col = "beta",
  se_col = "se",
  pval_col = "pval",
  eaf_col = "EAF",
  effect_allele_col = "effect_allele",
  other_allele_col = "other_allele"
)
exposure_dat$phenotype <- "F-ALL"

outcome_data_raw <- outcome_list[["gdm_subsamp"]]

# Convert outcome
outcome_dat <- format_data(
  outcome_data_raw,
  type = "outcome",
  snp_col = "SNP",
  beta_col = "beta",
  se_col = "se",
  pval_col = "pval",
  eaf_col = "eaf",
  effect_allele_col = "effect_allele",
  other_allele_col = "other_allele",
  phenotype = "gdm_subsamp"
)

harmonised <- harmonise_data(exposure_dat, outcome_dat)
harmonised$exposure <- "F-ALL"
harmonised$outcome <- "gdm_subsamp"
harmonised_list[[paste0("F-ALL", "_", "gdm_subsamp")]] <- harmonised

if (nrow(harmonised) == 0) next

mr_result <- mr(harmonised)
mr_result$exposure <- "F-ALL"
mr_result$outcome <- "gdm_subsamp"

all_results[[paste0("F-ALL", "_", "gdm_subsamp")]] <- mr_result

####



outcome_list_sensitivity <- split(stu_final, stu_final$Phenotype)
harmonised_list_sensitivity <- list()
all_results_sensitivity <- list()

for (exp_name in names(exposure_list)) {
  
  # Convert exposure
  exposure_dat <- format_data(
    exposure_list[[exp_name]],
    type = "exposure",
    snp_col = "RSID",
    beta_col = "beta",
    se_col = "se",
    pval_col = "pval",
    eaf_col = "EAF",
    effect_allele_col = "effect_allele",
    other_allele_col = "other_allele"
  )
  exposure_dat$phenotype <- exp_name
  
  # Loop outcome
  for (outcome_name in names(outcome_list_sensitivity)) {
    
    outcome_data_raw <- outcome_list_sensitivity[[outcome_name]]
    
    # Convert outcome
    outcome_dat <- format_data(
      outcome_data_raw,
      type = "outcome",
      snp_col = "SNP",
      beta_col = "beta",
      se_col = "se",
      pval_col = "pval",
      eaf_col = "eaf",
      effect_allele_col = "effect_allele",
      other_allele_col = "other_allele",
      phenotype = outcome_name
    )
    
    # harmonise
    harmonised <- harmonise_data(exposure_dat, outcome_dat)
    harmonised$exposure <- exp_name
    harmonised$outcome <- outcome_name
    harmonised_list_sensitivity[[paste0(exp_name, "_", outcome_name)]] <- harmonised
    
    if (nrow(harmonised) == 0) next
    
    # MR 
    mr_result <- mr(harmonised)
    mr_result$exposure <- exp_name
    mr_result$outcome <- outcome_name
    
    all_results_sensitivity[[paste0(exp_name, "_", outcome_name)]] <- mr_result
  }
}

all_mr_summary_sensitivity <- dplyr::bind_rows(all_results)

# Check harmonise list


#

MR Plot
```{r}
# IVW method
method_results <- all_mr_summary_sensitivity %>%
  filter(method %in% c("Inverse variance weighted"))

ivw <- method_results %>%
  filter(method == "Inverse variance weighted") %>%
  mutate(
    OR = exp(b),
    lower = exp(b - 1.96 * se),
    upper = exp(b + 1.96 * se),
    A = case_when(
      outcome == "hdp_subsamp" ~ "HDP",
      outcome == "gh_subsamp" ~ "GH",
      outcome == "gdm_subsamp" ~ "GDM",
      outcome == "pe_subsamp" ~ "PE",
      outcome == "misc_subsamp" ~ "Miscarriage",
      outcome == "sb_subsamp" ~ "Stillbirth",
      outcome == "pretb_all" ~ "Preterm Birth",
      outcome == "vpretb_all" ~ "Very Preterm Birth",
      outcome == "sga" ~ "Small GA",
      outcome == "lga" ~ "Large GA",
      outcome == "lbw_all" ~ "Low Birth Weight",
      outcome == "hbw_all" ~ "High Birth Weight",
      TRUE ~ outcome
    )
  )

ivw$outcome <- factor(ivw$A, levels = unique(ivw$A))
ivw$label <- paste(ivw$exposure, ivw$outcome, sep = " - ")
# Volcano Plot
volcano_plot <- ggplot(
  ivw, aes(x = OR, y = -log10(pval), color = exposure)) +
  geom_point(size = 3) +
  geom_vline(xintercept = 1, linetype = "dashed") +  
  geom_hline(yintercept = -log10(0.05), linetype = "dotted", color = "red") +  
  geom_text_repel(
    data = subset(ivw, pval < 0.05),  
    aes(label = label),
    size = 3,
    max.overlaps = 20
  ) +
  labs(
    x = "Odds Ratio (OR)", 
    y = "-log10(p-value)",
    title = "Volcano Plot of Sensitivity Analysis",
    color = "Exposure Type"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5)
  )

# Save plot

ggsave("sen_overall_ivw_voplot.png", plot = volcano_plot, width = 8, height = 6, dpi = 300)

# Forest Plot
forest_plot <- ggplot(
  ivw, aes(x = OR, y = A, color = exposure)) +
  geom_point(position = position_dodge(width = 0.7), size = 2) +
  geom_errorbarh(
    aes(xmin = lower, xmax = upper),
    height = 0.2,
    position = position_dodge(width = 0.7)
  ) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "gray40") +
  coord_cartesian(xlim = c(0, 3)) +  
  labs(
    title = "Forest Plot of Sensitivity Analysis",
    x = "OR",
    y = "Pregnancy Outcome",
    color = "Exposure"
  ) +
  theme_bw()

# Save plot


ggsave("sen_overall_ivw_foplot_ld_clumped.png", plot = forest_plot, width = 8, height = 6, dpi = 300)



ivw_significant <- ivw %>%
  filter(pval < 0.05)


```

Calculate F statistics
```{r}
fall <- fall %>%
  mutate(F_stat = (beta^2) / (se^2))
fanat <- fanat %>%
  mutate(F_stat = (beta^2) / (se^2))
fanov <- fanov %>%
  mutate(F_stat = (beta^2) / (se^2))
fexcl <- fexcl %>%
  mutate(F_stat = (beta^2) / (se^2))
fincl <- fincl %>%
  mutate(F_stat = (beta^2) / (se^2))

mean_f <- mean(fall$F_stat, na.rm = TRUE)
mean_f1 <- mean(fanat$F_stat, na.rm = TRUE)
mean_f2 <- mean(fanov$F_stat, na.rm = TRUE)
mean_f3 <- mean(fexcl$F_stat, na.rm = TRUE)
mean_f4 <- mean(fincl$F_stat, na.rm = TRUE)

```

MR Egger & Weighted Median & Cochran

results_sensitivity <- mr(harmonised_df_sensitivity,
                          method_list = c("mr_ivw", "mr_egger_regression", "mr_weighted_median"))
# Pleiotropy

pleiotropy_test <- mr_pleiotropy_test(harmonised_df_sensitivity)
significant_pleiotropy <- pleiotropy_test %>%
  filter(pval < 0.05)

# Cochran's Q
heterogeneity_results <- mr_heterogeneity(harmonised_df_sensitivity)
QQ <-subset(heterogeneity_results, Q_pval < 0.05)

# Leave one out analysis
leave1out_results <- mr_leaveoneout(harmonised_df_sensitivity)
mr_leaveoneout_plot(leave1out_results)

## Loop all exposures
output_dir <- "C:/Users/eo25503/OneDrive - University of Bristol/Documents/MR Preg stuff/leave_one_out_graphs"
if (!dir.exists(output_dir)) dir.create(output_dir)

exposures <- unique(leave1out_results$exposure)

lapply(exposures, function(exp_name) {
  
  leave1out_subset <- leave1out_results %>%
    filter(exposure == exp_name)
  
  overall_estimates <- leave1out_subset %>%
    filter(SNP == "All") %>%
    select(outcome, overall_b = b)
  
  leave1out_subset <- leave1out_subset %>%
    left_join(overall_estimates, by = "outcome")
  
  p <- ggplot(leave1out_subset, aes(x = b, y = SNP)) +
    geom_point() +
    geom_errorbarh(aes(xmin = b - 1.96 * se, xmax = b + 1.96 * se), height = 0.2) +
    geom_vline(aes(xintercept = overall_b), linetype = "dashed", color = "red") +
    facet_wrap(~ outcome, scales = "free_y") +
    labs(
      title = paste("Leave-One-Out Analysis for", exp_name),
      x = "Beta After Removing Each SNP",
      y = "SNP"
    ) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5))
  
  ggsave(
    filename = paste0(output_dir, "/", gsub("[^A-Za-z0-9]", "_", exp_name), "_leave1out_plot.png"),
    plot = p,
    width = 10,
    height = 8
  )
})
```

#####


ivw <- all_mr_summary %>%
  filter(method == "Inverse variance weighted") %>%
  mutate(
    OR = exp(b),
    lower = exp(b - 1.96 * se),
    upper = exp(b + 1.96 * se),
    A = case_when(
      outcome == "hdp_subsamp" ~ "HDP",
      outcome == "gh_subsamp" ~ "GH",
      outcome == "gdm_subsamp" ~ "GDM",
      outcome == "pe_subsamp" ~ "PE",
      outcome == "misc_subsamp" ~ "Miscarriage",
      outcome == "sb_subsamp" ~ "Stillbirth",
      outcome == "pretb_all" ~ "Preterm Birth",
      outcome == "vpretb_all" ~ "Very Preterm Birth",
      outcome == "sga" ~ "Small GA",
      outcome == "lga" ~ "Large GA",
      outcome == "lbw_all" ~ "Low Birth Weight",
      outcome == "hbw_all" ~ "High Birth Weight",
      TRUE ~ outcome
    )
  )


anat2 <- anat %>%
  mutate(
    OR = exp(b),
    lower = exp(b - 1.96 * se),
    upper = exp(b + 1.96 * se),
  )

volcano_plot <- ggplot(
  ivw, aes(x = OR, y = -log10(pval), color = exposure)) +
  geom_point(size = 3) +
  geom_vline(xintercept = 1, linetype = "dashed") +  
  geom_hline(yintercept = -log10(0.05), linetype = "dotted", color = "red") +  
  geom_text_repel(
    data = subset(ivw, pval < 0.05),  
    aes(label = label),
    size = 3,
    max.overlaps = 20
  ) +
  labs(
    x = "Odds Ratio (OR)", 
    y = "-log10(p-value)",
    title = "Volcano Plot of MR Results",
    color = "Exposure Type"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5)
  )
volcano_plot


create_volcano_plot <- function(data, exposure_name) {
  data_filtered <- data %>%
    filter(exposure == exposure_name) %>%
    mutate(significant = ifelse(pval < 0.05, "Significant", "Not Significant"))
  
  plot <- ggplot(
    data_filtered, aes(x = OR, y = -log10(pval), color = significant)) +
    geom_point(size = 3) +
    geom_vline(xintercept = 1, linetype = "dashed") +
    geom_hline(yintercept = -log10(0.05), linetype = "dotted", color = "red") +
    geom_text_repel(
      data = subset(data_filtered, pval < 0.05),
      aes(label = label),
      size = 3,
      max.overlaps = 20
    ) +
    scale_color_manual(values = c("Significant" = "red", "Not Significant" = "grey")) +
    labs(
      x = "Odds Ratio (OR)",
      y = "-log10(p-value)",
      title = paste("Volcano Plot of MR Results for", exposure_name),
      color = "Significance"
    ) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5))
  
  return(plot)
}


if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

for (exposure in exposure_types) {
  p <- create_volcano_plot(ivw, exposure)
  volcano_plots[[exposure]] <- p
  
  filename <- paste0(gsub("-", "", exposure), "_ivw_voplot.png")
  ggsave(file.path(output_dir, filename), plot = p, width = 8, height = 6, dpi = 300)
}

####

output_dir <- "C:/Users/eo25503/OneDrive - University of Bristol/Documents/MR Preg stuff/Forest_Plots"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

create_forest_plot <- function(data, exposure_name) {
  data_filtered <- data %>% filter(exposure == exposure_name)
  
  plot <- ggplot(
    data_filtered, aes(x = OR, y = A, color = exposure)) +
    geom_point(position = position_dodge(width = 0.7), size = 2) +
    geom_errorbarh(
      aes(xmin = lower, xmax = upper),
      height = 0.2,
      position = position_dodge(width = 0.7)
    ) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "gray40") +
    labs(
      title = paste("ORs of Female", exposure_name, "on APPOs"),
      x = "OR",
      y = "Pregnancy Outcome",
      color = "Exposure"
    ) +
    theme_bw()
  
  return(plot)
}


exposure_types <- unique(ivw$exposure)

for (exposure in exposure_types) {
  plot <- create_forest_plot(ivw, exposure)
  filename <- paste0(gsub("-", "", exposure), "_ivw_foplot.png")
  ggsave(file.path(output_dir, filename), plot = plot, width = 8, height = 6, dpi = 300)
}

####

anat3 <- anat2 %>%
  mutate(significant = ifelse(pval < 0.05, "Significant", "Not Significant"))

volcano_plot <- ggplot(
  anat3, aes(x = OR, y = -log10(pval), color = significant)) +
  geom_point(size = 3) +
  geom_vline(xintercept = 1, linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dotted", color = "red") +
  geom_text_repel(
    data = subset(anat3, pval < 0.05),
    aes(label = label),
    size = 3,
    max.overlaps = 20
  ) +
  scale_color_manual(
    values = c("Significant" = "red", "Not Significant" = "grey")) +
  labs(
    x = "Odds Ratio (OR)",
    y = "-log10(p-value)",
    title = "Volcano Plot of MR Results for F-ANAT",
    color = "Significance"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5)
  )

ggsave(
  filename = "C:/Users/eo25503/OneDrive - University of Bristol/Documents/MR Preg stuff/Volcano_Plots/FANAT-WR-voplot.png",
  plot = volcano_plot,
  width = 8,
  height = 6,
  dpi = 300
)

forest_plot <-ggplot(
  anat3, aes(x = OR, y = A, color = exposure)) +
  geom_point(position = position_dodge(width = 0.7), size = 2) +
  geom_errorbarh(
    aes(xmin = lower, xmax = upper),
    height = 0.2,
    position = position_dodge(width = 0.7)
  ) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "gray40") +
  labs(
    title = "Forest plot of MR results for F-ANAT",
    x = "OR",
    y = "Pregnancy Outcome",
    color = "Exposure"
  ) +
  theme_bw()

ggsave(
  filename = "C:/Users/eo25503/OneDrive - University of Bristol/Documents/MR Preg stuff/Forest_Plots/FANAT-WR-foplot.png",
  plot = forest_plot,
  width = 8,
  height = 6,
  dpi = 300
)

anat3 <- anat3 %>% mutate(A = case_when(
  outcome == "hdp_subsamp" ~ "HDP",
  outcome == "gh_subsamp" ~ "GH",
  outcome == "gdm_subsamp" ~ "GDM",
  outcome == "pe_subsamp" ~ "PE",
  outcome == "misc_subsamp" ~ "Miscarriage",
  outcome == "sb_subsamp" ~ "Stillbirth",
  outcome == "pretb_all" ~ "Preterm Birth",
  outcome == "vpretb_all" ~ "Very Preterm Birth",
  outcome == "sga" ~ "Small GA",
  outcome == "lga" ~ "Large GA",
  outcome == "lbw_all" ~ "Low Birth Weight",
  outcome == "hbw_all" ~ "High Birth Weight",
  TRUE ~ outcome
))

#####

write.table(egger, 
            file = file.path("C:/Users/eo25503/OneDrive - University of Bristol/Documents/MR Preg stuff/Updated Sensitivity analysis", "sen_egger.txt"),
            sep = "\t", 
            row.names = FALSE, 
            quote = FALSE)

write_xlsx(egger, "C:/Users/eo25503/OneDrive - University of Bristol/Documents/MR Preg stuff/Updated Sensitivity analysis/sen_egger.xlsx")


forest_plot1 <- ggplot(
  wm, aes(x = OR, y = outcome, color = exposure)) +
  geom_point(position = position_dodge(width = 0.7), size = 2) +
  geom_errorbarh(
    aes(xmin = lower, xmax = upper),
    height = 0.2,
    position = position_dodge(width = 0.7)
  ) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "gray40") +
  labs(
    title = "Forest plot of Weighted Median",
    x = "OR",
    y = "Pregnancy Outcome",
    color = "Exposure"
  ) +
  theme_bw() +
  coord_cartesian(xlim = c(0, 3))

for (this_exp in all_exposures) {
  
  outcomes <- unique(allmr$outcome[allmr$exposure == this_exp])
  
  plot_list <- list()
  
  for (this_out in outcomes) {
    mr_sub <- subset(allmr, exposure == this_exp & outcome == this_out)
    dat_sub <- subset(harmon, exposure == this_exp & outcome == this_out)
    
    mr_sub <- subset(mr_sub, method %in% c("Inverse variance weighted", "MR Egger", "Weighted median"))
    
    if (nrow(dat_sub) >= 3 & nrow(mr_sub) > 0) {
      p <- mr_scatter_plot(mr_sub, dat_sub)[[1]] +
        ggtitle(paste0(this_exp, " → ", this_out))
      plot_list[[this_out]] <- p
    }
  }
  
  if (length(plot_list) > 0) {
    combined_plot <- wrap_plots(plot_list, ncol = 4, nrow = 3) + 
      plot_annotation(title = paste0("MR Scatter Plots for ", this_exp))
    
    file_name <- paste0("scatter_", gsub("[^A-Za-z0-9]", "_", this_exp), "_combined.png")
    
    ggsave(filename = file.path(output_dir, file_name),
           plot = combined_plot,
           width = 4 * 4,   
           height = 4 * 3,  
           dpi = 300)
  }
}

for (i in seq_along(plots)) {
  plot_obj <- plots[[i]]
  
  this_exp <- unique(plot_obj$data$exposure)
  
  plot_obj <- plot_obj + ggtitle(this_exp)
  
  file_name <- paste0(gsub("[^A-Za-z0-9]", "_", this_exp), "_funnel_plot.png")
  
  ggsave(filename = file.path(output_dir, file_name),
         plot = plot_obj,
         width = 6, height = 5, dpi = 300)
}

lapply(exposures, function(exp_name) {
  
  leave1out_subset <- leave1out_results %>%
    filter(exposure == exp_name)
  
  overall_estimates <- leave1out_subset %>%
    filter(SNP == "All") %>%
    select(outcome, overall_b = b)
  
  leave1out_subset <- leave1out_subset %>%
    left_join(overall_estimates, by = "outcome")
  
  p <- ggplot(leave1out_subset, aes(x = b, y = SNP)) +
    geom_point() +
    geom_errorbarh(aes(xmin = b - 1.96 * se, xmax = b + 1.96 * se), height = 0.2) +
    geom_vline(aes(xintercept = overall_b), linetype = "dashed", color = "red") +
    facet_wrap(~ outcome, scales = "free_y") +
    labs(
      title = paste("Leave-One-Out Analysis for", exp_name),
      x = "Beta After Removing Each SNP",
      y = "SNP"
    ) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5))
  
  ggsave(
    filename = paste0(output_dir, "/", gsub("[^A-Za-z0-9]", "_", exp_name), "_leave1out_plot.png"),
    plot = p,
    width = 10,
    height = 8
  )
})

####

for (exp_name in names(exposure_list)) {
  
  # Convert exposure
  exposure_dat <- format_data(
    exposure_list[[exp_name]],
    type = "exposure",
    snp_col = "RSID",
    beta_col = "beta",
    se_col = "se",
    pval_col = "pval",
    eaf_col = "EAF",
    effect_allele_col = "effect_allele",
    other_allele_col = "other_allele"
  )
  exposure_dat$phenotype <- exp_name
  
  # Loop outcome
  for (outcome_name in names(outcome_list)) {
    
    outcome_data_raw <- outcome_list[[outcome_name]]
    
    for (study_name in unique(outcome_data_raw$study)) {
      
    outcome_study_specific <- outcome_data_raw[outcome_data_raw$study == study_name,]
    
    # Convert outcome
    outcome_dat <- format_data(
      outcome_study_specific,
      type = "outcome",
      snp_col = "SNP",
      beta_col = "beta",
      se_col = "se",
      pval_col = "pval",
      eaf_col = "eaf",
      effect_allele_col = "effect_allele",
      other_allele_col = "other_allele",
      phenotype = outcome_name
    )
    
    # harmonise
    harmonised <- harmonise_data(exposure_dat, outcome_dat)
    harmonised$exposure <- exp_name
    harmonised$outcome <- outcome_name
    harmonised_list[[paste0(exp_name, "_", outcome_name,"_",study_name)]] <- harmonised
    
    if (nrow(harmonised) == 0) next
    
    # MR 
    mr_result <- mr(harmonised)
    mr_result$exposure <- exp_name
    mr_result$outcome <- outcome_name
    mr_result$study <- study_name
    
    all_results[[paste0(exp_name, "_", outcome_name,"_",study_name)]] <- mr_result
    }
  }
}

forest_plot <- ggplot(
  ivw, aes(x = OR, y = label2, color = study)) +
  geom_point(position = position_dodge(width = 0.7), size = 2) +
  geom_errorbarh(
    aes(xmin = lower, xmax = upper),
    height = 0.2,
    position = position_dodge(width = 0.7)
  ) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "gray40") +
  coord_cartesian(xlim = c(0, 3)) +  
  labs(
    title = "Forest Plot of Sensitivity Analysis",
    x = "OR",
    y = "Pregnancy Exposure and Outcome",
    color = "Study"
  ) +
  theme_bw()

leave1out_results <- leave1out_results %>% mutate(A = case_when(
  id.outcome =="Lv2reN"~"F-ALL_gdm_subsamp_FinnGen",
  id.outcome =="VB11Gs"~"F-ALL_gdm_subsamp_MOBA",
  id.outcome =="B5AIqG"~"F-ALL_gdm_subsamp_Public",
  id.outcome =="B3kG6k"~"F-ALL_gh_subsamp_ALSPAC",
  id.outcome =="0YgLFr"~"F-ALL_gh_subsamp_BIB-SA",
  id.outcome =="qHTrhQ"~"F-ALL_gh_subsamp_BIB-WE",
  id.outcome =="h6OLWJ"~"F-ALL_gh_subsamp_FinnGen",
  id.outcome =="igyBky"~"F-ALL_gh_subsamp_MOBA",
  id.outcome =="P7y1tI"~"F-ALL_gh_subsamp_UKB",
  id.outcome =="2M6KbJ"~"F-ALL_hbw_all_ALSPAC",
  id.outcome =="Yg0pY1"~"F-ALL_hbw_all_BIB-SA",
  id.outcome =="4XYAOQ"~"F-ALL_hbw_all_BIB-WE",
  id.outcome =="XxUIk1"~"F-ALL_hbw_all_MOBA",
  id.outcome =="66Uola"~"F-ALL_hbw_all_UKB",
  id.outcome =="lnAeEd"~"F-ALL_hdp_subsamp_ALSPAC",
  id.outcome =="XbTj9q"~"F-ALL_hdp_subsamp_BIB-SA",
  id.outcome =="7vvFzj"~"F-ALL_hdp_subsamp_BIB-WE",
  id.outcome =="Sclfxb"~"F-ALL_hdp_subsamp_FinnGen",
  id.outcome =="u3oCpj"~"F-ALL_hdp_subsamp_MOBA",
  id.outcome =="09Z93s"~"F-ALL_hdp_subsamp_UKB",
  id.outcome =="FdsbYu"~"F-ALL_lbw_all_ALSPAC",
  id.outcome =="41F90d"~"F-ALL_lbw_all_BIB-SA",
  id.outcome =="CEbdWF"~"F-ALL_lbw_all_BIB-WE",
  id.outcome =="SjZICV"~"F-ALL_lbw_all_MOBA",
  id.outcome =="IPZGTg"~"F-ALL_lbw_all_UKB",
  id.outcome =="utkLVT"~"F-ALL_lga_ALSPAC",
  id.outcome =="WWzocQ"~"F-ALL_lga_BIB-SA",
  id.outcome =="J34vlo"~"F-ALL_lga_BIB-WE",
  id.outcome =="C2JbT0"~"F-ALL_lga_MOBA",
  id.outcome =="fDKvMb"~"F-ALL_lga_UKB",
  id.outcome =="Fn6fVm"~"F-ALL_misc_subsamp_ALSPAC",
  id.outcome =="gasrnT"~"F-ALL_misc_subsamp_BIB-SA",
  id.outcome =="WCPGBi"~"F-ALL_misc_subsamp_BIB-WE",
  id.outcome =="h38D0x"~"F-ALL_misc_subsamp_FinnGen",
  id.outcome =="XSVpi4"~"F-ALL_misc_subsamp_MOBA",
  id.outcome =="IYNtgi"~"F-ALL_misc_subsamp_UKB",
  id.outcome =="IbNs97"~"F-ALL_pe_subsamp_BIB-SA",
  id.outcome =="Ht2ltk"~"F-ALL_pe_subsamp_BIB-WE",
  id.outcome =="vuh4LV"~"F-ALL_pe_subsamp_FinnGen",
  id.outcome =="ST0smw"~"F-ALL_pe_subsamp_MOBA",
  id.outcome =="3LTlpR"~"F-ALL_pe_subsamp_UKB",
  id.outcome =="HFJYiF"~"F-ALL_pe_subsamp_Public",
  id.outcome =="SDXV7k"~"F-ALL_pretb_all_BIB-SA",
  id.outcome =="VzKgyl"~"F-ALL_pretb_all_BIB-WE",
  id.outcome =="ZR4hFV"~"F-ALL_pretb_all_FinnGen",
  id.outcome =="Bs2xZG"~"F-ALL_pretb_all_MOBA",
  id.outcome =="BaYO3u"~"F-ALL_pretb_all_UKB",
  id.outcome =="UBER6n"~"F-ALL_pretb_all_Public",
  id.outcome =="07UXxO"~"F-ALL_sb_subsamp_ALSPAC",
  id.outcome =="ohsOMW"~"F-ALL_sb_subsamp_BIB-SA",
  id.outcome =="lOTEf9"~"F-ALL_sb_subsamp_BIB-WE",
  id.outcome =="d9BeYE"~"F-ALL_sb_subsamp_MOBA",
  id.outcome =="7HeLXa"~"F-ALL_sb_subsamp_UKB",
  id.outcome =="E4ByVT"~"F-ALL_sga_ALSPAC",
  id.outcome =="TRXWDv"~"F-ALL_sga_BIB-SA",
  id.outcome =="1H3pMB"~"F-ALL_sga_BIB-WE",
  id.outcome =="HPlG9Z"~"F-ALL_sga_MOBA",
  id.outcome =="7nm5Cu"~"F-ALL_sga_UKB",
  id.outcome =="9BjQ8H"~"F-ALL_vpretb_all_ALSPAC",
  id.outcome =="vxKocQ"~"F-ALL_vpretb_all_BIB-SA",
  id.outcome =="XnwoeG"~"F-ALL_vpretb_all_BIB-WE",
  id.outcome =="hogAzv"~"F-ALL_vpretb_all_MOBA",
  id.outcome =="kvFVqa"~"F-ALL_vpretb_all_UKB",
  id.outcome =="Fqt33u"~"F-INCL_gdm_subsamp_FinnGen",
  id.outcome =="7fedQu"~"F-INCL_gdm_subsamp_MOBA",
  id.outcome =="MyEnmI"~"F-INCL_gdm_subsamp_Public",
  id.outcome =="ve2g1z"~"F-INCL_gh_subsamp_ALSPAC",
  id.outcome =="gF30IE"~"F-INCL_gh_subsamp_BIB-SA",
  id.outcome =="7fq22t"~"F-INCL_gh_subsamp_BIB-WE",
  id.outcome =="iRV6fr"~"F-INCL_gh_subsamp_FinnGen",
  id.outcome =="eDQJmV"~"F-INCL_gh_subsamp_MOBA",
  id.outcome =="1tiko4"~"F-INCL_gh_subsamp_UKB",
  id.outcome =="rh2oIz"~"F-INCL_hbw_all_ALSPAC",
  id.outcome =="8QDKPa"~"F-INCL_hbw_all_BIB-SA",
  id.outcome =="P3VsPq"~"F-INCL_hbw_all_BIB-WE",
  id.outcome =="86UINE"~"F-INCL_hbw_all_MOBA",
  id.outcome =="DV10s6"~"F-INCL_hbw_all_UKB",
  id.outcome =="MdigN9"~"F-INCL_hdp_subsamp_ALSPAC",
  id.outcome =="X3AJbt"~"F-INCL_hdp_subsamp_BIB-SA",
  id.outcome =="1K0Dnv"~"F-INCL_hdp_subsamp_BIB-WE",
  id.outcome =="Yw7l3h"~"F-INCL_hdp_subsamp_FinnGen",
  id.outcome =="Zo23YL"~"F-INCL_hdp_subsamp_MOBA",
  id.outcome =="uDCH1E"~"F-INCL_hdp_subsamp_UKB",
  id.outcome =="lBMAQD"~"F-INCL_lbw_all_ALSPAC",
  id.outcome =="UTMzLJ"~"F-INCL_lbw_all_BIB-SA",
  id.outcome =="sZqLUG"~"F-INCL_lbw_all_BIB-WE",
  id.outcome =="Ihntub"~"F-INCL_lbw_all_MOBA",
  id.outcome =="tViVB2"~"F-INCL_lbw_all_UKB",
  id.outcome =="kvNQIV"~"F-INCL_lga_ALSPAC",
  id.outcome =="yb6m1Y"~"F-INCL_lga_BIB-SA",
  id.outcome =="VjN6Wa"~"F-INCL_lga_BIB-WE",
  id.outcome =="QE82f6"~"F-INCL_lga_MOBA",
  id.outcome =="cYSHvC"~"F-INCL_lga_UKB",
  id.outcome =="oVaZ4v"~"F-INCL_misc_subsamp_ALSPAC",
  id.outcome =="iEEgSL"~"F-INCL_misc_subsamp_BIB-SA",
  id.outcome =="ajvPUa"~"F-INCL_misc_subsamp_BIB-WE",
  id.outcome =="kequkS"~"F-INCL_misc_subsamp_FinnGen",
  id.outcome =="8PrJTH"~"F-INCL_misc_subsamp_MOBA",
  id.outcome =="sahQPh"~"F-INCL_misc_subsamp_UKB",
  id.outcome =="5H5i5u"~"F-INCL_pe_subsamp_BIB-SA",
  id.outcome =="0l7bXd"~"F-INCL_pe_subsamp_BIB-WE",
  id.outcome =="1b9hQB"~"F-INCL_pe_subsamp_FinnGen",
  id.outcome =="xOJYzH"~"F-INCL_pe_subsamp_MOBA",
  id.outcome =="wugaU8"~"F-INCL_pe_subsamp_UKB",
  id.outcome =="P7rGkk"~"F-INCL_pe_subsamp_Public",
  id.outcome =="mKT3h5"~"F-INCL_pretb_all_BIB-SA",
  id.outcome =="2aVkvp"~"F-INCL_pretb_all_BIB-WE",
  id.outcome =="RdVUsq"~"F-INCL_pretb_all_FinnGen",
  id.outcome =="qkzXuA"~"F-INCL_pretb_all_MOBA",
  id.outcome =="EgZtDw"~"F-INCL_pretb_all_UKB",
  id.outcome =="AnFPU3"~"F-INCL_pretb_all_Public",
  id.outcome =="yDq6tU"~"F-INCL_sb_subsamp_ALSPAC",
  id.outcome =="IfL5e7"~"F-INCL_sb_subsamp_BIB-SA",
  id.outcome =="BRYKBs"~"F-INCL_sb_subsamp_BIB-WE",
  id.outcome =="mhgcaZ"~"F-INCL_sb_subsamp_MOBA",
  id.outcome =="35rR3g"~"F-INCL_sb_subsamp_UKB",
  id.outcome =="YxFYKG"~"F-INCL_sga_ALSPAC",
  id.outcome =="HxaPkP"~"F-INCL_sga_BIB-SA",
  id.outcome =="7PzujE"~"F-INCL_sga_BIB-WE",
  id.outcome =="1cbj2D"~"F-INCL_sga_MOBA",
  id.outcome =="FK2WC8"~"F-INCL_sga_UKB",
  id.outcome =="xkTgPb"~"F-INCL_vpretb_all_ALSPAC",
  id.outcome =="1TSlHM"~"F-INCL_vpretb_all_BIB-SA",
  id.outcome =="wbqbTB"~"F-INCL_vpretb_all_BIB-WE",
  id.outcome =="N1gX3x"~"F-INCL_vpretb_all_MOBA",
  id.outcome =="Vr9CO5"~"F-INCL_vpretb_all_UKB",
  id.outcome =="KcSw5K"~"F-ANOV_gdm_subsamp_FinnGen",
  id.outcome =="e59pEt"~"F-ANOV_gdm_subsamp_MOBA",
  id.outcome =="ZHiwgG"~"F-ANOV_gdm_subsamp_Public",
  id.outcome =="99SkGk"~"F-ANOV_gh_subsamp_ALSPAC",
  id.outcome =="ywYdjQ"~"F-ANOV_gh_subsamp_BIB-SA",
  id.outcome =="HszVnO"~"F-ANOV_gh_subsamp_BIB-WE",
  id.outcome =="udKrlp"~"F-ANOV_gh_subsamp_FinnGen",
  id.outcome =="D1yjwi"~"F-ANOV_gh_subsamp_MOBA",
  id.outcome =="4wNTVE"~"F-ANOV_gh_subsamp_UKB",
  id.outcome =="Ldp2GJ"~"F-ANOV_hbw_all_ALSPAC",
  id.outcome =="CqcC76"~"F-ANOV_hbw_all_BIB-SA",
  id.outcome =="L5RBKn"~"F-ANOV_hbw_all_BIB-WE",
  id.outcome =="Ftgwg0"~"F-ANOV_hbw_all_MOBA",
  id.outcome =="UGUi6l"~"F-ANOV_hbw_all_UKB",
  id.outcome =="tyaXa4"~"F-ANOV_hdp_subsamp_ALSPAC",
  id.outcome =="PuWpue"~"F-ANOV_hdp_subsamp_BIB-SA",
  id.outcome =="MW0doX"~"F-ANOV_hdp_subsamp_BIB-WE",
  id.outcome =="6ktTje"~"F-ANOV_hdp_subsamp_FinnGen",
  id.outcome =="AMWKJg"~"F-ANOV_hdp_subsamp_MOBA",
  id.outcome =="uu5syd"~"F-ANOV_hdp_subsamp_UKB",
  id.outcome =="wbh7HT"~"F-ANOV_lbw_all_ALSPAC",
  id.outcome =="cNgOny"~"F-ANOV_lbw_all_BIB-SA",
  id.outcome =="XtD4OK"~"F-ANOV_lbw_all_BIB-WE",
  id.outcome =="wgyrIs"~"F-ANOV_lbw_all_MOBA",
  id.outcome =="QeFVM2"~"F-ANOV_lbw_all_UKB",
  id.outcome =="H1Ugxq"~"F-ANOV_lga_ALSPAC",
  id.outcome =="McuXyg"~"F-ANOV_lga_BIB-SA",
  id.outcome =="jj15tg"~"F-ANOV_lga_BIB-WE",
  id.outcome =="cjfac1"~"F-ANOV_lga_MOBA",
  id.outcome =="h6rK0A"~"F-ANOV_lga_UKB",
  id.outcome =="owY8EY"~"F-ANOV_misc_subsamp_ALSPAC",
  id.outcome =="eaQADV"~"F-ANOV_misc_subsamp_BIB-SA",
  id.outcome =="vX4J3X"~"F-ANOV_misc_subsamp_BIB-WE",
  id.outcome =="EIKQQM"~"F-ANOV_misc_subsamp_FinnGen",
  id.outcome =="gBVcN7"~"F-ANOV_misc_subsamp_MOBA",
  id.outcome =="xCyFbs"~"F-ANOV_misc_subsamp_UKB",
  id.outcome =="Pkm7Rk"~"F-ANOV_pe_subsamp_BIB-SA",
  id.outcome =="SEWaSy"~"F-ANOV_pe_subsamp_BIB-WE",
  id.outcome =="eJxjUU"~"F-ANOV_pe_subsamp_FinnGen",
  id.outcome =="2HSLAV"~"F-ANOV_pe_subsamp_MOBA",
  id.outcome =="kUgJKL"~"F-ANOV_pe_subsamp_UKB",
  id.outcome =="ZiaBxU"~"F-ANOV_pe_subsamp_Public",
  id.outcome =="MgtjN1"~"F-ANOV_pretb_all_BIB-SA",
  id.outcome =="YBZJiK"~"F-ANOV_pretb_all_BIB-WE",
  id.outcome =="YIUqXY"~"F-ANOV_pretb_all_FinnGen",
  id.outcome =="s3nOhZ"~"F-ANOV_pretb_all_MOBA",
  id.outcome =="r0CsZI"~"F-ANOV_pretb_all_UKB",
  id.outcome =="PhZHGC"~"F-ANOV_pretb_all_Public",
  id.outcome =="cwy1js"~"F-ANOV_sb_subsamp_ALSPAC",
  id.outcome =="LIY08W"~"F-ANOV_sb_subsamp_BIB-SA",
  id.outcome =="nm86Gw"~"F-ANOV_sb_subsamp_BIB-WE",
  id.outcome =="6Hi4FP"~"F-ANOV_sb_subsamp_MOBA",
  id.outcome =="m5S5cp"~"F-ANOV_sb_subsamp_UKB",
  id.outcome =="6EBYug"~"F-ANOV_sga_ALSPAC",
  id.outcome =="ccn2qa"~"F-ANOV_sga_BIB-SA",
  id.outcome =="2782mj"~"F-ANOV_sga_BIB-WE",
  id.outcome =="Mvmtqu"~"F-ANOV_sga_MOBA",
  id.outcome =="J2ar7L"~"F-ANOV_sga_UKB",
  id.outcome =="U6opw6"~"F-ANOV_vpretb_all_ALSPAC",
  id.outcome =="e1lOU8"~"F-ANOV_vpretb_all_BIB-SA",
  id.outcome =="0hBXHr"~"F-ANOV_vpretb_all_BIB-WE",
  id.outcome =="BE4NGE"~"F-ANOV_vpretb_all_MOBA",
  id.outcome =="XyjBFh"~"F-ANOV_vpretb_all_UKB",
  id.outcome =="k1v4kR"~"F-EXCL_gdm_subsamp_FinnGen",
  id.outcome =="cPL94O"~"F-EXCL_gdm_subsamp_MOBA",
  id.outcome =="kEDG88"~"F-EXCL_gdm_subsamp_Public",
  id.outcome =="xMqSTr"~"F-EXCL_gh_subsamp_ALSPAC",
  id.outcome =="woLTPd"~"F-EXCL_gh_subsamp_BIB-SA",
  id.outcome =="GqLdD9"~"F-EXCL_gh_subsamp_BIB-WE",
  id.outcome =="USI1kA"~"F-EXCL_gh_subsamp_FinnGen",
  id.outcome =="WZ7VyQ"~"F-EXCL_gh_subsamp_MOBA",
  id.outcome =="lPvr3P"~"F-EXCL_gh_subsamp_UKB",
  id.outcome =="yMd3vh"~"F-EXCL_hbw_all_ALSPAC",
  id.outcome =="CeGg5r"~"F-EXCL_hbw_all_BIB-SA",
  id.outcome =="G4cTlM"~"F-EXCL_hbw_all_BIB-WE",
  id.outcome =="XPvkZt"~"F-EXCL_hbw_all_MOBA",
  id.outcome =="KG2ZJ6"~"F-EXCL_hbw_all_UKB",
  id.outcome =="RTpftL"~"F-EXCL_hdp_subsamp_ALSPAC",
  id.outcome =="3FUoVh"~"F-EXCL_hdp_subsamp_BIB-SA",
  id.outcome =="tjArAx"~"F-EXCL_hdp_subsamp_BIB-WE",
  id.outcome =="t2rK0L"~"F-EXCL_hdp_subsamp_FinnGen",
  id.outcome =="JQPkVb"~"F-EXCL_hdp_subsamp_MOBA",
  id.outcome =="v4X5aO"~"F-EXCL_hdp_subsamp_UKB",
  id.outcome =="kBuVau"~"F-EXCL_lbw_all_ALSPAC",
  id.outcome =="ExhEnR"~"F-EXCL_lbw_all_BIB-SA",
  id.outcome =="IlHvGU"~"F-EXCL_lbw_all_BIB-WE",
  id.outcome =="rOSoZy"~"F-EXCL_lbw_all_MOBA",
  id.outcome =="IrUnS4"~"F-EXCL_lbw_all_UKB",
  id.outcome =="c21eWl"~"F-EXCL_lga_ALSPAC",
  id.outcome =="9tvIGg"~"F-EXCL_lga_BIB-SA",
  id.outcome =="rZXj6b"~"F-EXCL_lga_BIB-WE",
  id.outcome =="zTsv21"~"F-EXCL_lga_MOBA",
  id.outcome =="UwbWrn"~"F-EXCL_lga_UKB",
  id.outcome =="P102C0"~"F-EXCL_misc_subsamp_ALSPAC",
  id.outcome =="ZfEJDS"~"F-EXCL_misc_subsamp_BIB-SA",
  id.outcome =="QytqCX"~"F-EXCL_misc_subsamp_BIB-WE",
  id.outcome =="7PW1cv"~"F-EXCL_misc_subsamp_FinnGen",
  id.outcome =="PfFUvm"~"F-EXCL_misc_subsamp_MOBA",
  id.outcome =="wN0Lvk"~"F-EXCL_misc_subsamp_UKB",
  id.outcome =="BVGjqA"~"F-EXCL_pe_subsamp_BIB-SA",
  id.outcome =="xd5vzD"~"F-EXCL_pe_subsamp_BIB-WE",
  id.outcome =="Gjvf5a"~"F-EXCL_pe_subsamp_FinnGen",
  id.outcome =="o7Co7N"~"F-EXCL_pe_subsamp_MOBA",
  id.outcome =="bFDueJ"~"F-EXCL_pe_subsamp_UKB",
  id.outcome =="6CUkz5"~"F-EXCL_pe_subsamp_Public",
  id.outcome =="AJyDoX"~"F-EXCL_pretb_all_BIB-SA",
  id.outcome =="bLfpak"~"F-EXCL_pretb_all_BIB-WE",
  id.outcome =="AhFHCJ"~"F-EXCL_pretb_all_FinnGen",
  id.outcome =="W4CLo0"~"F-EXCL_pretb_all_MOBA",
  id.outcome =="h51h24"~"F-EXCL_pretb_all_UKB",
  id.outcome =="IzLFCH"~"F-EXCL_pretb_all_Public",
  id.outcome =="9hRQZ3"~"F-EXCL_sb_subsamp_ALSPAC",
  id.outcome =="n0YgTC"~"F-EXCL_sb_subsamp_BIB-SA",
  id.outcome =="c3jCD5"~"F-EXCL_sb_subsamp_BIB-WE",
  id.outcome =="gVmIAV"~"F-EXCL_sb_subsamp_MOBA",
  id.outcome =="58QZo6"~"F-EXCL_sb_subsamp_UKB",
  id.outcome =="kwLkyn"~"F-EXCL_sga_ALSPAC",
  id.outcome =="9vhDta"~"F-EXCL_sga_BIB-SA",
  id.outcome =="zxNV39"~"F-EXCL_sga_BIB-WE",
  id.outcome =="avZ0ov"~"F-EXCL_sga_MOBA",
  id.outcome =="H2DTU4"~"F-EXCL_sga_UKB",
  id.outcome =="LE0l6k"~"F-EXCL_vpretb_all_ALSPAC",
  id.outcome =="7imcd6"~"F-EXCL_vpretb_all_BIB-SA",
  id.outcome =="Ru788n"~"F-EXCL_vpretb_all_BIB-WE",
  id.outcome =="VEyZh4"~"F-EXCL_vpretb_all_MOBA",
  id.outcome =="WBtA9I"~"F-EXCL_vpretb_all_UKB",
  id.outcome =="G9xFd0"~"F-ANAT_gdm_subsamp_FinnGen",
  id.outcome =="Qf1tWo"~"F-ANAT_gdm_subsamp_MOBA",
  id.outcome =="LSwC7B"~"F-ANAT_gdm_subsamp_Public",
  id.outcome =="aDorF8"~"F-ANAT_gh_subsamp_ALSPAC",
  id.outcome =="VHOLbB"~"F-ANAT_gh_subsamp_BIB-SA",
  id.outcome =="tuKV4o"~"F-ANAT_gh_subsamp_BIB-WE",
  id.outcome =="7JaT9W"~"F-ANAT_gh_subsamp_FinnGen",
  id.outcome =="0zefbA"~"F-ANAT_gh_subsamp_MOBA",
  id.outcome =="lZDgM1"~"F-ANAT_gh_subsamp_UKB",
  id.outcome =="99Zdx6"~"F-ANAT_hbw_all_ALSPAC",
  id.outcome =="xMlX3c"~"F-ANAT_hbw_all_BIB-SA",
  id.outcome =="J4HoCq"~"F-ANAT_hbw_all_BIB-WE",
  id.outcome =="7vaPtn"~"F-ANAT_hbw_all_MOBA",
  id.outcome =="RmIHr6"~"F-ANAT_hbw_all_UKB",
  id.outcome =="Pp1NgF"~"F-ANAT_hdp_subsamp_ALSPAC",
  id.outcome =="e0XIU0"~"F-ANAT_hdp_subsamp_BIB-SA",
  id.outcome =="fmWkeK"~"F-ANAT_hdp_subsamp_BIB-WE",
  id.outcome =="jXWWld"~"F-ANAT_hdp_subsamp_FinnGen",
  id.outcome =="L2ALuW"~"F-ANAT_hdp_subsamp_MOBA",
  id.outcome =="o25e7x"~"F-ANAT_hdp_subsamp_UKB",
  id.outcome =="latmgX"~"F-ANAT_lbw_all_ALSPAC",
  id.outcome =="6KwWPK"~"F-ANAT_lbw_all_BIB-SA",
  id.outcome =="nnoKQw"~"F-ANAT_lbw_all_BIB-WE",
  id.outcome =="HoAZ3X"~"F-ANAT_lbw_all_MOBA",
  id.outcome =="RvBpsj"~"F-ANAT_lbw_all_UKB",
  id.outcome =="siWtjC"~"F-ANAT_lga_ALSPAC",
  id.outcome =="YKizr7"~"F-ANAT_lga_BIB-SA",
  id.outcome =="skrrzQ"~"F-ANAT_lga_BIB-WE",
  id.outcome =="EvHkqU"~"F-ANAT_lga_MOBA",
  id.outcome =="yjJW25"~"F-ANAT_lga_UKB",
  id.outcome =="KJfFMT"~"F-ANAT_misc_subsamp_ALSPAC",
  id.outcome =="Fa9ck9"~"F-ANAT_misc_subsamp_BIB-SA",
  id.outcome =="yvcvYW"~"F-ANAT_misc_subsamp_BIB-WE",
  id.outcome =="Kkc8QH"~"F-ANAT_misc_subsamp_FinnGen",
  id.outcome =="YPXbda"~"F-ANAT_misc_subsamp_MOBA",
  id.outcome =="JwjfQo"~"F-ANAT_misc_subsamp_UKB",
  id.outcome =="Bx4Shy"~"F-ANAT_pe_subsamp_BIB-SA",
  id.outcome =="lDfMPV"~"F-ANAT_pe_subsamp_BIB-WE",
  id.outcome =="1znWK5"~"F-ANAT_pe_subsamp_FinnGen",
  id.outcome =="cPUomq"~"F-ANAT_pe_subsamp_MOBA",
  id.outcome =="5bEk2X"~"F-ANAT_pe_subsamp_UKB",
  id.outcome =="VjxJb3"~"F-ANAT_pe_subsamp_Public",
  id.outcome =="n0TgJ2"~"F-ANAT_pretb_all_BIB-SA",
  id.outcome =="oQ1rFq"~"F-ANAT_pretb_all_BIB-WE",
  id.outcome =="YDhQoj"~"F-ANAT_pretb_all_FinnGen",
  id.outcome =="uGxYOK"~"F-ANAT_pretb_all_MOBA",
  id.outcome =="RtyqAB"~"F-ANAT_pretb_all_UKB",
  id.outcome =="U5QKfm"~"F-ANAT_pretb_all_Public",
  id.outcome =="5xrLFK"~"F-ANAT_sb_subsamp_ALSPAC",
  id.outcome =="P8DMPd"~"F-ANAT_sb_subsamp_BIB-SA",
  id.outcome =="yV9lsC"~"F-ANAT_sb_subsamp_BIB-WE",
  id.outcome =="HZnSfF"~"F-ANAT_sb_subsamp_MOBA",
  id.outcome =="hSaAto"~"F-ANAT_sb_subsamp_UKB",
  id.outcome =="7Zh4su"~"F-ANAT_sga_ALSPAC",
  id.outcome =="wmQq6W"~"F-ANAT_sga_BIB-SA",
  id.outcome =="p2KCye"~"F-ANAT_sga_BIB-WE",
  id.outcome =="BUKz8F"~"F-ANAT_sga_MOBA",
  id.outcome =="evJNTR"~"F-ANAT_sga_UKB",
  id.outcome =="9iAJS3"~"F-ANAT_vpretb_all_ALSPAC",
  id.outcome =="YzzfX4"~"F-ANAT_vpretb_all_BIB-SA",
  id.outcome =="iJVNEK"~"F-ANAT_vpretb_all_BIB-WE",
  id.outcome =="c5byzs"~"F-ANAT_vpretb_all_MOBA",
  id.outcome =="5gaFtW"~"F-ANAT_vpretb_all_UKB")
)

lapply(exposures, function(exp_name) {
  
  leave1out_subset <- leave1out_results %>%
    filter(exposure == exp_name)
  
  overall_estimates <- leave1out_subset %>%
    filter(SNP == "All") %>%
    select(A, overall_b = b)
  
  leave1out_subset <- leave1out_subset %>%
    left_join(overall_estimates, by = "A")
  
  p <- ggplot(leave1out_subset, aes(x = b, y = SNP)) +
    geom_point() +
    geom_errorbarh(aes(xmin = b - 1.96 * se, xmax = b + 1.96 * se), height = 0.2) +
    geom_vline(aes(xintercept = overall_b), linetype = "dashed", color = "red") +
    facet_wrap(~ A, scales = "free_y") +
    labs(
      title = paste("Leave-One-Out Analysis for", exp_name),
      x = "Beta After Removing Each SNP",
      y = "SNP"
    ) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5))
  
  ggsave(
    filename = paste0(output_dir, "/", gsub("[^A-Za-z0-9]", "_", exp_name), "_leave1out_plot.png"),
    plot = p,
    width = 25,
    height = 20
  )
})
#####

lapply(exposures, function(exp_name) {
  
  ivw2 <- ivw %>%
    filter(exposure == exp_name)

forest_plot <- ggplot(
  ivw2, aes(x = OR, y = A, color = study)) +
  geom_point(position = position_dodge(width = 0.7), size = 2) +
  geom_errorbarh(
    aes(xmin = lower, xmax = upper),
    height = 0.2,
    position = position_dodge(width = 0.7)
  ) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "gray40") +
  coord_cartesian(xlim = c(0, 3)) +  
  labs(
    title = paste(exp_name, "Forest Plot of Sensitivity Analysis"),
    x = paste("OR for", exp_name),
    y = "Pregnancy Outcome",
    color = "Study"
  ) + geom_stripped_rows() +
  theme_bw()

ggsave(
  filename = paste0(output_dir, "/", gsub("[^A-Za-z0-9]", "_", exp_name), "_forest_plot.png"),
  plot = forest_plot,
  width = 25,
  height = 20
)
})

#####

fall <- fall %>%
  mutate(R_sqr = 2*(beta^2)*EAF*(1-EAF))
fanat <- fanat %>%
  mutate(R_sqr = 2*(beta^2)*EAF*(1-EAF))
fanov <- fanov %>%
  mutate(R_sqr = 2*(beta^2)*EAF*(1-EAF))
fexcl <- fexcl %>%
  mutate(R_sqr = 2*(beta^2)*EAF*(1-EAF))
fincl <- fincl %>%
  mutate(R_sqr = 2*(beta^2)*EAF*(1-EAF))

Rsqr_f <- sum(fall$R_sqr, na.rm = TRUE)
Rsqr_f1 <- sum(fanat$R_sqr, na.rm = TRUE)
Rsqr_f2 <- sum(fanov$R_sqr, na.rm = TRUE)
Rsqr_f3 <- sum(fexcl$R_sqr, na.rm = TRUE)
Rsqr_f4 <- sum(fincl$R_sqr, na.rm = TRUE)

#####

for (i in seq_along(plots)) {
  plot_obj <- plots[[i]]
  
  this_exp <- unique(plot_obj$data$exposure)
  
  this_out <- unique(plot_obj$data$outcome)
  
  plot_obj <- plot_obj + ggtitle(paste(this_exp,"and",this_out,sep = " "))
  
  file_name <- paste0(gsub("[^A-Za-z0-9]", "_", this_exp),"_", this_out,"_funnel_plot.png")
  
  ggsave(filename = file.path(output_dir, file_name),
         plot = plot_obj,
         width = 6, height = 5, dpi = 300)
}

#####

prep_func <- function(results_file) {
  
  #rename outcomes
  results_file <- results_file %>% mutate(A = case_when(
    outcome == "hdp_subsamp" ~ "HDP",
    outcome == "gh_subsamp" ~ "GH",
    outcome == "gdm_subsamp" ~ "GDM",
    outcome == "pe_subsamp" ~ "PE",
    outcome == "misc_subsamp" ~ "Miscarriage",
    outcome == "sb_subsamp" ~ "Stillbirth",
    outcome == "pretb_all" ~ "Preterm Birth",
    outcome == "vpretb_all" ~ "Very Preterm Birth",
    outcome == "sga" ~ "Small GA",
    outcome == "lga" ~ "Large GA",
    outcome == "lbw_all" ~ "Low Birth Weight",
    outcome == "hbw_all" ~ "High Birth Weight",
    TRUE ~ outcome
  ))
  
  #reorder studies
  as.factor(results_file$study)
  results_file$study <- factor(results_file$study, levels = c("Overall", "Public" ,"UKB", "MOBA", 
                                                              "FinnGen", "BIB-WE", "BIB-SA", "ALSPAC"))
  return(results_file)
  
}


myforestplot <- function(df, log_T = TRUE, xlab, limits)
{
  x <- forestplot(
    df = df,
    estimate = b,    # b and se have already been multiplied by 10 
    se = se,
    pvalue = pval,
    name = A,
    logodds = log_T,
    colour = study,
    #title = title,
    xlab = xlab,
    xlim= limits, 
    shape = study
  ) 
  colours_BP <- c("black", "gold2", "deepskyblue1", "coral", "springgreen3", "purple", "magenta", "turquoise")
  shapes_BP <- c(23L, 21L, 21L, 21L, 21L, 21L, 21L, 21L)
  x <- x + scale_color_manual(values=colours_BP)
  #x <- x + scale_x_continuous(breaks = c(0.5, 1, 1.5, 2, 2.5, 3, 5))
  x <- x + ggplot2::scale_shape_manual(values = shapes_BP)
  print(x)
}

SBP_DBP_func <- function(A, B, C, D){
  plot.new()
  par(mar=c(1,1,1,1), mgp=c(3,1,0))
  x <- ggarrange(A, B, C, D, labels=c('A', 'B','C','D'),
                 ncol = 2, nrow = 2, common.legend = T, hjust = -3, legend = "bottom")
  print(x)
}

######

fall <- read.csv("exposure_FALL_clumped.csv")
fanat <- read.csv("exposure_FANAT_clumped.csv")
fanov <- read.csv("exposure_FANOV_clumped.csv")
fexcl <- read.csv("exposure_FEXCL_clumped.csv")
fincl <- read.csv("exposure_FINCL_clumped.csv")

for (exp_name in names(exposure_list)) {
  
  # Convert exposure
  exposure_dat <- format_data(
    exposure_list[[exp_name]],
    type = "exposure",
    snp_col = "rsid",
    beta_col = "beta",
    se_col = "standard_error",
    pval_col = "pval",
    eaf_col = "effect_allele_frequency",
    effect_allele_col = "effect_allele",
    other_allele_col = "other_allele"
  )
  exposure_dat$phenotype <- exp_name
  
  # Loop outcome
  for (outcome_name in names(outcome_list)) {
    
    outcome_data_raw <- outcome_list[[outcome_name]]
    
    # Convert outcome
    outcome_dat <- format_data(
      outcome_data_raw,
      type = "outcome",
      snp_col = "SNP",
      beta_col = "beta",
      se_col = "se",
      pval_col = "pval",
      eaf_col = "eaf",
      effect_allele_col = "effect_allele",
      other_allele_col = "other_allele",
      phenotype = outcome_name
    )
    
    # harmonise
    harmonised <- harmonise_data(exposure_dat, outcome_dat)
    harmonised$exposure <- exp_name
    harmonised$outcome <- outcome_name
    harmonised_list[[paste0(exp_name, "_", outcome_name)]] <- harmonised
    
    if (nrow(harmonised) == 0) next
    
    # MR 
    mr_result <- mr(harmonised)
    mr_result$exposure <- exp_name
    mr_result$outcome <- outcome_name
    
    all_results[[paste0(exp_name, "_", outcome_name)]] <- mr_result
  }
}

ivw <- all_mr_summary %>%
  filter(method == "Inverse variance weighted") %>%
  mutate(
    OR = exp(b),
    lower = exp(b - 1.96 * se),
    upper = exp(b + 1.96 * se),
    A = case_when(
      outcome == "hdp_subsamp" ~ "HDP",
      outcome == "gh_subsamp" ~ "GH",
      outcome == "gdm_subsamp" ~ "GDM",
      outcome == "pe_subsamp" ~ "PE",
      outcome == "misc_subsamp" ~ "Miscarriage",
      outcome == "sb_subsamp" ~ "Stillbirth",
      outcome == "pretb_all" ~ "Preterm Birth",
      outcome == "vpretb_all" ~ "Very Preterm Birth",
      outcome == "sga" ~ "Small GA",
      outcome == "lga" ~ "Large GA",
      outcome == "lbw_all" ~ "Low Birth Weight",
      outcome == "hbw_all" ~ "High Birth Weight",
      TRUE ~ outcome
    )
  )

for (exposure in exposure_types) {
  p <- create_volcano_plot(ivw, exposure)
  volcano_plots[[exposure]] <- p
  
  filename <- paste0(gsub("-", "", exposure), "_ivw_voplot_LD_clumped.png")
  ggsave(file.path(output_dir, filename), plot = p, width = 8, height = 6, dpi = 300)
}

for (exposure in exposure_types) {
  plot <- create_forest_plot(ivw, exposure)
  filename <- paste0(gsub("-", "", exposure), "_ivw_foplot_LD_clumped.png")
  ggsave(file.path(output_dir, filename), plot = plot, width = 8, height = 6, dpi = 300)
}

#####Doing 4-way Plot 

myforestplot <- function(df, log_T = TRUE, xlab, limits)
{
  x <- print(plot_name)

    df = df,
    estimate = b,     
    se = se,
    pvalue = pval,
    name = outcome,
    logodds = log_T,
    #title = title,
    xlab = xlab,
    xlim= limits
  ) + theme_forest(base_size = 10) + ggforce::facet_col(
    facets = ~group,
    scales = "free_y",
    space = "free"
  )
  colours_BP <- c("black", "gold2", "deepskyblue1", "coral", "springgreen3", "purple", "magenta", "turquoise")
  shapes_BP <- c(23L, 21L, 21L, 21L, 21L, 21L, 21L, 21L)
  x <- x + scale_color_manual(values=colours_BP)
  #x <- x + scale_x_continuous(breaks = c(0.5, 1, 1.5, 2, 2.5, 3, 5))
  x <- x + ggplot2::scale_shape_manual(values = shapes_BP)
  print(x)
}

plot_merger_func2 <- function(A, B, C, D){
  plot.new()
  par(mar=c(1,1,1,1), mgp=c(3,1,0))
  x <- ggarrange(A, B,C,D, labels=c('A', 'B','C','D'),
                 ncol = 2, nrow = 2, common.legend = T, hjust = -3, legend = "bottom")
  print(x)
}

save_func <- function(file_name, plot_name, height)
{
  png(file_name, res=330, height=height, width=6000)
  print(plot_name)
  dev.off()
}

myforestplot2 <- function(df, log_T = TRUE, xlab, limits)
{
  x <- forestplot(
    df = df,
    estimate = b,     
    se = se,
    pvalue = pval,
    name = outcome,
    logodds = log_T,
    #title = title,
    xlab = xlab,
    xlim= limits,
    colour = method,
    shape = method
  ) + theme_forest(base_size = 10) + ggforce::facet_col(
    facets = ~group,
    scales = "free_y",
    space = "free"
  )
  colours_BP <- c("black", "gold2", "deepskyblue1", "coral", "springgreen3", "purple", "magenta", "turquoise")
  shapes_BP <- c(23L, 21L, 21L, 21L, 21L, 21L, 21L, 21L)
  x <- x + scale_color_manual(values=colours_BP)
  #x <- x + scale_x_continuous(breaks = c(0.5, 1, 1.5, 2, 2.5, 3, 5))
  x <- x + ggplot2::scale_shape_manual(values = shapes_BP)
  print(x)
}

####

myforestplot4 <- function(df, log_T = TRUE, xlab, limits)
{
  x <- forestplot(
    df = df,
    estimate = b,    # b and se have already been multiplied by 10 
    se = se,
    pvalue = pval,
    name = A,
    logodds = log_T,
    colour = adjustment,
    #title = title,
    xlab = xlab,
    xlim= limits, 
    shape = adjustment
  ) + theme_forest(base_size = 10) +
    ggforce::facet_col(
      facets = ~group,
      scales = "free_y",
      space = "free"
    )
  colours_BP <- c("black", "gold2", "deepskyblue1", "coral", "springgreen3", "purple", "magenta", "turquoise")
  shapes_BP <- c(23L, 21L, 21L, 21L, 21L, 21L, 21L, 21L)
  x <- x + scale_color_manual(values=colours_BP)
  #x <- x + scale_x_continuous(breaks = c(0.5, 1, 1.5, 2, 2.5, 3, 5))
  x <- x + ggplot2::scale_shape_manual(values = shapes_BP)
  print(x)
}


####

for (exp_name in names(exposure_list)) {
  
  # Convert exposure
  exposure_dat <- format_data(
    exposure_list[[exp_name]],
    type = "exposure",
    snp_col = "RSID",
    beta_col = "beta",
    se_col = "se",
    pval_col = "pval",
    eaf_col = "EAF",
    effect_allele_col = "effect_allele",
    other_allele_col = "other_allele"
  )
  exposure_dat$phenotype <- exp_name
  
  # Loop outcome
  for (outcome_name in names(outcome_list_trios)) {
    
    outcome_data_raw <- outcome_list_trios[[outcome_name]]
    
    # Convert outcome
    outcome_dat <- format_data(
      outcome_data_raw,
      type = "outcome",
      snp_col = "SNP",
      beta_col = "beta_mat",
      se_col = "se_mat",
      pval_col = "p_mat",
      eaf_col = "eaf_mat",
      effect_allele_col = "effect_allele",
      other_allele_col = "other_allele",
      phenotype = outcome_name
    )
    
    # harmonise
    harmonised <- harmonise_data(exposure_dat, outcome_dat)
    harmonised$exposure <- exp_name
    harmonised$outcome <- outcome_name
    harmonised_list_mat_trios[[paste0(exp_name, "_", outcome_name)]] <- harmonised
    
    if (nrow(harmonised) == 0) next
    
    # MR 
    mr_result <- mr(harmonised)
    mr_result$exposure <- exp_name
    mr_result$outcome <- outcome_name
    
    all_results_mat_trios[[paste0(exp_name, "_", outcome_name)]] <- mr_result
  }
}


method_results_mat_donut_trios_clean <- method_results_mat_donut_trios %>%
  mutate(
    b = round(b, 2),
    se = round(se, 2),
    pval = round(pval, 4),
    outcome_1 = case_when(
      outcome == "hdp_subsamp" ~ "HDP",
      outcome == "gh_subsamp" ~ "GH",
      outcome == "gdm_subsamp" ~ "GDM",
      outcome == "pe_subsamp" ~ "PE",
      outcome == "misc_subsamp" ~ "Miscarriage",
      outcome == "sb_subsamp" ~ "Stillbirth",
      outcome == "pretb_all" ~ "Preterm Birth",
      outcome == "vpretb_all" ~ "Very Preterm Birth",
      outcome == "sga" ~ "Small GA",
      outcome == "lga" ~ "Large GA",
      outcome == "lbw_all" ~ "Low Birth Weight",
      outcome == "hbw_all" ~ "High Birth Weight",
      TRUE ~ outcome)
  ) %>%
  select(-id.exposure, -id.outcome, -outcome)

####

pleiotropy_test2 <- pleiotropy_test %>% mutate(study = case_when(
  id.outcome =="Lv2reN"~"FinnGen",
  id.outcome =="VB11Gs"~"MOBA",
  id.outcome =="B5AIqG"~"Public",
  id.outcome =="B3kG6k"~"ALSPAC",
  id.outcome =="0YgLFr"~"BIB-SA",
  id.outcome =="qHTrhQ"~"BIB-WE",
  id.outcome =="h6OLWJ"~"FinnGen",
  id.outcome =="igyBky"~"MOBA",
  id.outcome =="P7y1tI"~"UKB",
  id.outcome =="2M6KbJ"~"ALSPAC",
  id.outcome =="Yg0pY1"~"BIB-SA",
  id.outcome =="4XYAOQ"~"BIB-WE",
  id.outcome =="XxUIk1"~"MOBA",
  id.outcome =="66Uola"~"UKB",
  id.outcome =="lnAeEd"~"ALSPAC",
  id.outcome =="XbTj9q"~"BIB-SA",
  id.outcome =="7vvFzj"~"BIB-WE",
  id.outcome =="Sclfxb"~"FinnGen",
  id.outcome =="u3oCpj"~"MOBA",
  id.outcome =="09Z93s"~"UKB",
  id.outcome =="FdsbYu"~"ALSPAC",
  id.outcome =="41F90d"~"BIB-SA",
  id.outcome =="CEbdWF"~"BIB-WE",
  id.outcome =="SjZICV"~"MOBA",
  id.outcome =="IPZGTg"~"UKB",
  id.outcome =="utkLVT"~"ALSPAC",
  id.outcome =="WWzocQ"~"BIB-SA",
  id.outcome =="J34vlo"~"BIB-WE",
  id.outcome =="C2JbT0"~"MOBA",
  id.outcome =="fDKvMb"~"UKB",
  id.outcome =="Fn6fVm"~"ALSPAC",
  id.outcome =="gasrnT"~"BIB-SA",
  id.outcome =="WCPGBi"~"BIB-WE",
  id.outcome =="h38D0x"~"FinnGen",
  id.outcome =="XSVpi4"~"MOBA",
  id.outcome =="IYNtgi"~"UKB",
  id.outcome =="IbNs97"~"BIB-SA",
  id.outcome =="Ht2ltk"~"BIB-WE",
  id.outcome =="vuh4LV"~"FinnGen",
  id.outcome =="ST0smw"~"MOBA",
  id.outcome =="3LTlpR"~"UKB",
  id.outcome =="HFJYiF"~"Public",
  id.outcome =="SDXV7k"~"BIB-SA",
  id.outcome =="VzKgyl"~"BIB-WE",
  id.outcome =="ZR4hFV"~"FinnGen",
  id.outcome =="Bs2xZG"~"MOBA",
  id.outcome =="BaYO3u"~"UKB",
  id.outcome =="UBER6n"~"Public",
  id.outcome =="07UXxO"~"ALSPAC",
  id.outcome =="ohsOMW"~"BIB-SA",
  id.outcome =="lOTEf9"~"BIB-WE",
  id.outcome =="d9BeYE"~"MOBA",
  id.outcome =="7HeLXa"~"UKB",
  id.outcome =="E4ByVT"~"ALSPAC",
  id.outcome =="TRXWDv"~"BIB-SA",
  id.outcome =="1H3pMB"~"BIB-WE",
  id.outcome =="HPlG9Z"~"MOBA",
  id.outcome =="7nm5Cu"~"UKB",
  id.outcome =="9BjQ8H"~"ALSPAC",
  id.outcome =="vxKocQ"~"BIB-SA",
  id.outcome =="XnwoeG"~"BIB-WE",
  id.outcome =="hogAzv"~"MOBA",
  id.outcome =="kvFVqa"~"UKB",
  id.outcome =="Fqt33u"~"FinnGen",
  id.outcome =="7fedQu"~"MOBA",
  id.outcome =="MyEnmI"~"Public",
  id.outcome =="ve2g1z"~"ALSPAC",
  id.outcome =="gF30IE"~"BIB-SA",
  id.outcome =="7fq22t"~"BIB-WE",
  id.outcome =="iRV6fr"~"FinnGen",
  id.outcome =="eDQJmV"~"MOBA",
  id.outcome =="1tiko4"~"UKB",
  id.outcome =="rh2oIz"~"ALSPAC",
  id.outcome =="8QDKPa"~"BIB-SA",
  id.outcome =="P3VsPq"~"BIB-WE",
  id.outcome =="86UINE"~"MOBA",
  id.outcome =="DV10s6"~"UKB",
  id.outcome =="MdigN9"~"ALSPAC",
  id.outcome =="X3AJbt"~"BIB-SA",
  id.outcome =="1K0Dnv"~"BIB-WE",
  id.outcome =="Yw7l3h"~"FinnGen",
  id.outcome =="Zo23YL"~"MOBA",
  id.outcome =="uDCH1E"~"UKB",
  id.outcome =="lBMAQD"~"ALSPAC",
  id.outcome =="UTMzLJ"~"BIB-SA",
  id.outcome =="sZqLUG"~"BIB-WE",
  id.outcome =="Ihntub"~"MOBA",
  id.outcome =="tViVB2"~"UKB",
  id.outcome =="kvNQIV"~"ALSPAC",
  id.outcome =="yb6m1Y"~"BIB-SA",
  id.outcome =="VjN6Wa"~"BIB-WE",
  id.outcome =="QE82f6"~"MOBA",
  id.outcome =="cYSHvC"~"UKB",
  id.outcome =="oVaZ4v"~"ALSPAC",
  id.outcome =="iEEgSL"~"BIB-SA",
  id.outcome =="ajvPUa"~"BIB-WE",
  id.outcome =="kequkS"~"FinnGen",
  id.outcome =="8PrJTH"~"MOBA",
  id.outcome =="sahQPh"~"UKB",
  id.outcome =="5H5i5u"~"BIB-SA",
  id.outcome =="0l7bXd"~"BIB-WE",
  id.outcome =="1b9hQB"~"FinnGen",
  id.outcome =="xOJYzH"~"MOBA",
  id.outcome =="wugaU8"~"UKB",
  id.outcome =="P7rGkk"~"Public",
  id.outcome =="mKT3h5"~"BIB-SA",
  id.outcome =="2aVkvp"~"BIB-WE",
  id.outcome =="RdVUsq"~"FinnGen",
  id.outcome =="qkzXuA"~"MOBA",
  id.outcome =="EgZtDw"~"UKB",
  id.outcome =="AnFPU3"~"Public",
  id.outcome =="yDq6tU"~"ALSPAC",
  id.outcome =="IfL5e7"~"BIB-SA",
  id.outcome =="BRYKBs"~"BIB-WE",
  id.outcome =="mhgcaZ"~"MOBA",
  id.outcome =="35rR3g"~"UKB",
  id.outcome =="YxFYKG"~"ALSPAC",
  id.outcome =="HxaPkP"~"BIB-SA",
  id.outcome =="7PzujE"~"BIB-WE",
  id.outcome =="1cbj2D"~"MOBA",
  id.outcome =="FK2WC8"~"UKB",
  id.outcome =="xkTgPb"~"ALSPAC",
  id.outcome =="1TSlHM"~"BIB-SA",
  id.outcome =="wbqbTB"~"BIB-WE",
  id.outcome =="N1gX3x"~"MOBA",
  id.outcome =="Vr9CO5"~"UKB",
  id.outcome =="KcSw5K"~"FinnGen",
  id.outcome =="e59pEt"~"MOBA",
  id.outcome =="ZHiwgG"~"Public",
  id.outcome =="99SkGk"~"ALSPAC",
  id.outcome =="ywYdjQ"~"BIB-SA",
  id.outcome =="HszVnO"~"BIB-WE",
  id.outcome =="udKrlp"~"FinnGen",
  id.outcome =="D1yjwi"~"MOBA",
  id.outcome =="4wNTVE"~"UKB",
  id.outcome =="Ldp2GJ"~"ALSPAC",
  id.outcome =="CqcC76"~"BIB-SA",
  id.outcome =="L5RBKn"~"BIB-WE",
  id.outcome =="Ftgwg0"~"MOBA",
  id.outcome =="UGUi6l"~"UKB",
  id.outcome =="tyaXa4"~"ALSPAC",
  id.outcome =="PuWpue"~"BIB-SA",
  id.outcome =="MW0doX"~"BIB-WE",
  id.outcome =="6ktTje"~"FinnGen",
  id.outcome =="AMWKJg"~"MOBA",
  id.outcome =="uu5syd"~"UKB",
  id.outcome =="wbh7HT"~"ALSPAC",
  id.outcome =="cNgOny"~"BIB-SA",
  id.outcome =="XtD4OK"~"BIB-WE",
  id.outcome =="wgyrIs"~"MOBA",
  id.outcome =="QeFVM2"~"UKB",
  id.outcome =="H1Ugxq"~"ALSPAC",
  id.outcome =="McuXyg"~"BIB-SA",
  id.outcome =="jj15tg"~"BIB-WE",
  id.outcome =="cjfac1"~"MOBA",
  id.outcome =="h6rK0A"~"UKB",
  id.outcome =="owY8EY"~"ALSPAC",
  id.outcome =="eaQADV"~"BIB-SA",
  id.outcome =="vX4J3X"~"BIB-WE",
  id.outcome =="EIKQQM"~"FinnGen",
  id.outcome =="gBVcN7"~"MOBA",
  id.outcome =="xCyFbs"~"UKB",
  id.outcome =="Pkm7Rk"~"BIB-SA",
  id.outcome =="SEWaSy"~"BIB-WE",
  id.outcome =="eJxjUU"~"FinnGen",
  id.outcome =="2HSLAV"~"MOBA",
  id.outcome =="kUgJKL"~"UKB",
  id.outcome =="ZiaBxU"~"Public",
  id.outcome =="MgtjN1"~"BIB-SA",
  id.outcome =="YBZJiK"~"BIB-WE",
  id.outcome =="YIUqXY"~"FinnGen",
  id.outcome =="s3nOhZ"~"MOBA",
  id.outcome =="r0CsZI"~"UKB",
  id.outcome =="PhZHGC"~"Public",
  id.outcome =="cwy1js"~"ALSPAC",
  id.outcome =="LIY08W"~"BIB-SA",
  id.outcome =="nm86Gw"~"BIB-WE",
  id.outcome =="6Hi4FP"~"MOBA",
  id.outcome =="m5S5cp"~"UKB",
  id.outcome =="6EBYug"~"ALSPAC",
  id.outcome =="ccn2qa"~"BIB-SA",
  id.outcome =="2782mj"~"BIB-WE",
  id.outcome =="Mvmtqu"~"MOBA",
  id.outcome =="J2ar7L"~"UKB",
  id.outcome =="U6opw6"~"ALSPAC",
  id.outcome =="e1lOU8"~"BIB-SA",
  id.outcome =="0hBXHr"~"BIB-WE",
  id.outcome =="BE4NGE"~"MOBA",
  id.outcome =="XyjBFh"~"UKB",
  id.outcome =="k1v4kR"~"FinnGen",
  id.outcome =="cPL94O"~"MOBA",
  id.outcome =="kEDG88"~"Public",
  id.outcome =="xMqSTr"~"ALSPAC",
  id.outcome =="woLTPd"~"BIB-SA",
  id.outcome =="GqLdD9"~"BIB-WE",
  id.outcome =="USI1kA"~"FinnGen",
  id.outcome =="WZ7VyQ"~"MOBA",
  id.outcome =="lPvr3P"~"UKB",
  id.outcome =="yMd3vh"~"ALSPAC",
  id.outcome =="CeGg5r"~"BIB-SA",
  id.outcome =="G4cTlM"~"BIB-WE",
  id.outcome =="XPvkZt"~"MOBA",
  id.outcome =="KG2ZJ6"~"UKB",
  id.outcome =="RTpftL"~"ALSPAC",
  id.outcome =="3FUoVh"~"BIB-SA",
  id.outcome =="tjArAx"~"BIB-WE",
  id.outcome =="t2rK0L"~"FinnGen",
  id.outcome =="JQPkVb"~"MOBA",
  id.outcome =="v4X5aO"~"UKB",
  id.outcome =="kBuVau"~"ALSPAC",
  id.outcome =="ExhEnR"~"BIB-SA",
  id.outcome =="IlHvGU"~"BIB-WE",
  id.outcome =="rOSoZy"~"MOBA",
  id.outcome =="IrUnS4"~"UKB",
  id.outcome =="c21eWl"~"ALSPAC",
  id.outcome =="9tvIGg"~"BIB-SA",
  id.outcome =="rZXj6b"~"BIB-WE",
  id.outcome =="zTsv21"~"MOBA",
  id.outcome =="UwbWrn"~"UKB",
  id.outcome =="P102C0"~"ALSPAC",
  id.outcome =="ZfEJDS"~"BIB-SA",
  id.outcome =="QytqCX"~"BIB-WE",
  id.outcome =="7PW1cv"~"FinnGen",
  id.outcome =="PfFUvm"~"MOBA",
  id.outcome =="wN0Lvk"~"UKB",
  id.outcome =="BVGjqA"~"BIB-SA",
  id.outcome =="xd5vzD"~"BIB-WE",
  id.outcome =="Gjvf5a"~"FinnGen",
  id.outcome =="o7Co7N"~"MOBA",
  id.outcome =="bFDueJ"~"UKB",
  id.outcome =="6CUkz5"~"Public",
  id.outcome =="AJyDoX"~"BIB-SA",
  id.outcome =="bLfpak"~"BIB-WE",
  id.outcome =="AhFHCJ"~"FinnGen",
  id.outcome =="W4CLo0"~"MOBA",
  id.outcome =="h51h24"~"UKB",
  id.outcome =="IzLFCH"~"Public",
  id.outcome =="9hRQZ3"~"ALSPAC",
  id.outcome =="n0YgTC"~"BIB-SA",
  id.outcome =="c3jCD5"~"BIB-WE",
  id.outcome =="gVmIAV"~"MOBA",
  id.outcome =="58QZo6"~"UKB",
  id.outcome =="kwLkyn"~"ALSPAC",
  id.outcome =="9vhDta"~"BIB-SA",
  id.outcome =="zxNV39"~"BIB-WE",
  id.outcome =="avZ0ov"~"MOBA",
  id.outcome =="H2DTU4"~"UKB",
  id.outcome =="LE0l6k"~"ALSPAC",
  id.outcome =="7imcd6"~"BIB-SA",
  id.outcome =="Ru788n"~"BIB-WE",
  id.outcome =="VEyZh4"~"MOBA",
  id.outcome =="WBtA9I"~"UKB")
)

####

for (type in infertility_types) {
  subset_data <- subset(df2, infertility_type == type)
  
  filename <- paste0("exposure_2025_All_", gsub("[^A-Za-z0-9]", "_", type), ".txt")
  
  write.table(subset_data, file = filename, sep = "\t", quote = FALSE, row.names = FALSE)
}

for (exp_name in names(exposure_list)) {
  
  # Convert exposure
  exposure_dat <- format_data(
    exposure_list[[exp_name]],
    type = "exposure",
    snp_col = "RSID",
    beta_col = "beta",
    se_col = "se",
    pval_col = "pval",
    eaf_col = "EAF",
    effect_allele_col = "effect_allele",
    other_allele_col = "other_allele"
  )
  exposure_dat$phenotype <- exp_name
  
  # Loop outcome
  for (outcome_name in names(outcome_list)) {
    
    outcome_data_raw <- outcome_list[[outcome_name]]
    
    # Convert outcome
    outcome_dat <- format_data(
      outcome_data_raw,
      type = "outcome",
      snp_col = "SNP",
      beta_col = "beta",
      se_col = "se",
      pval_col = "pval",
      eaf_col = "eaf",
      effect_allele_col = "effect_allele",
      other_allele_col = "other_allele",
      phenotype = outcome_name
    )
    
    # harmonise
    harmonised <- harmonise_data(exposure_dat, outcome_dat)
    harmonised$exposure <- exp_name
    harmonised$outcome <- outcome_name
    harmonised_list[[paste0(exp_name, "_", outcome_name)]] <- harmonised
    
    if (nrow(harmonised) == 0) next
    
    # MR 
    mr_result <- mr(harmonised)
    mr_result$exposure <- exp_name
    mr_result$outcome <- outcome_name
    
    all_results[[paste0(exp_name, "_", outcome_name)]] <- mr_result
  }
}

ivw <- method_results %>%
  filter(method %in% c("Inverse variance weighted","Wald ratio")) %>%
  mutate(
    OR = exp(b),
    lower = exp(b - 1.96 * se),
    upper = exp(b + 1.96 * se),
    A = case_when(
      outcome == "hdp_subsamp" ~ "HDP",
      outcome == "gh_subsamp" ~ "GH",
      outcome == "gdm_subsamp" ~ "GDM",
      outcome == "pe_subsamp" ~ "PE",
      outcome == "misc_subsamp" ~ "Miscarriage",
      outcome == "sb_subsamp" ~ "Stillbirth",
      outcome == "pretb_all" ~ "Preterm Birth",
      outcome == "vpretb_all" ~ "Very Preterm Birth",
      outcome == "sga" ~ "Small GA",
      outcome == "lga" ~ "Large GA",
      outcome == "lbw_all" ~ "Low Birth Weight",
      outcome == "hbw_all" ~ "High Birth Weight",
      TRUE ~ outcome
    )
  )


myforestplot <- function(df, log_T = TRUE, xlab, limits)
{
  x <- forestplot(
    df = df,
    estimate = b,     
    se = se,
    pvalue = pval,
    name = outcome,
    logodds = log_T,
    #title = title,
    xlab = xlab,
    xlim= limits
  ) + theme_forest(base_size = 10) + ggforce::facet_col(
    facets = ~group,
    scales = "free_y",
    space = "free"
  )
  colours_BP <- c("black", "gold2", "deepskyblue1", "coral", "springgreen3", "purple", "magenta", "turquoise")
  shapes_BP <- c(23L, 21L, 21L, 21L, 21L, 21L, 21L, 21L)
  x <- x + scale_color_manual(values=colours_BP)
  #x <- x + scale_x_continuous(breaks = c(0.5, 1, 1.5, 2, 2.5, 3, 5))
  x <- x + ggplot2::scale_shape_manual(values = shapes_BP)
  print(x)
}

plot_merger_func <- function(A, B, C, D,E){
  plot.new()
  par(mar=c(1,1,1,1), mgp=c(3,1,0))
  x <- ggarrange(A, B,C,D,E, labels=c('A', 'B','C','D','E'),
                 ncol = 2, nrow = 3, common.legend = T, hjust = -3, legend = "bottom")
  print(x)
}

save_func <- function(file_name, plot_name, height)
{
  png(file_name, res=330, height=height, width=6000)
  print(plot_name)
  dev.off()
}

plot_merger_func2 <- function(A, B, C, D){
  plot.new()
  par(mar=c(1,1,1,1), mgp=c(3,1,0))
  x <- ggarrange(A, B,C,D, labels=c('A', 'B','C','D'),
                 ncol = 2, nrow = 2, common.legend = T, hjust = -3, legend = "bottom")
  print(x)
}


for (this_exp in all_exposures) {
  
  outcomes <- unique(all_mr_summary$outcome[all_mr_summary$exposure == this_exp])
  
  plot_list <- list()
  
  for (this_out in outcomes) {
    mr_sub <- subset(all_mr_summary, exposure == this_exp & outcome == this_out)
    dat_sub <- subset(harmonised_df, exposure == this_exp & outcome == this_out)
    
    mr_sub <- subset(mr_sub, method %in% c("Inverse variance weighted", "MR Egger", "Weighted median"))
    
    if (nrow(dat_sub) >= 3 & nrow(mr_sub) > 0) {
      p <- mr_scatter_plot(mr_sub, dat_sub)[[1]] +
        ggtitle(paste0(this_exp, " → ", this_out))
      plot_list[[this_out]] <- p
    }
  }
  
  if (length(plot_list) > 0) {
    combined_plot <- wrap_plots(plot_list, ncol = 4, nrow = 3) + 
      plot_annotation(title = paste0("MR Scatter Plots for ", this_exp))
    
    file_name <- paste0("scatter_", gsub("[^A-Za-z0-9]", "_", this_exp), "_combined.png")
    
    ggsave(filename = file.path(output_dir, file_name),
           plot = combined_plot,
           width = 4 * 4,   
           height = 4 * 3,  
           dpi = 300)
  }
}


for (exp_name in names(exposure_list)) {
  
  # Convert exposure
  exposure_dat <- format_data(
    exposure_list[[exp_name]],
    type = "exposure",
    snp_col = "RSID",
    beta_col = "beta",
    se_col = "se",
    pval_col = "pval",
    eaf_col = "EAF",
    effect_allele_col = "effect_allele",
    other_allele_col = "other_allele"
  )
  exposure_dat$phenotype <- exp_name
  
  # Loop outcome
  for (outcome_name in names(outcome_list2)) {
    
    outcome_data_raw <- outcome_list2[[outcome_name]]
    
    for (study_name in unique(outcome_data_raw$study)) {
      
      outcome_study_specific <- outcome_data_raw[outcome_data_raw$study == study_name,]
      
      # Convert outcome
      outcome_dat <- format_data(
        outcome_study_specific,
        type = "outcome",
        snp_col = "SNP",
        beta_col = "beta",
        se_col = "se",
        pval_col = "pval",
        eaf_col = "eaf",
        effect_allele_col = "effect_allele",
        other_allele_col = "other_allele",
        phenotype = outcome_name
      )
      
      # harmonise
      harmonised <- harmonise_data(exposure_dat, outcome_dat)
      harmonised$exposure <- exp_name
      harmonised$outcome <- outcome_name
      harmonised_list2[[paste0(exp_name, "_", outcome_name,"_",study_name)]] <- harmonised
      
      if (nrow(harmonised) == 0) next
      
      # MR 
      mr_result <- mr(harmonised)
      mr_result$exposure <- exp_name
      mr_result$outcome <- outcome_name
      mr_result$study <- study_name
      
      all_results2[[paste0(exp_name, "_", outcome_name,"_",study_name)]] <- mr_result
    }
  }
}

myforestplot3 <- function(df, log_T = TRUE, xlab, limits)
{
  x <- forestplot(
    df = df,
    estimate = b,     
    se = se,
    pvalue = pval,
    name = A,
    logodds = log_T,
    #title = title,
    xlab = xlab,
    xlim= limits,
    colour = study,
    shape = study
  ) + theme_forest(base_size = 10) + ggforce::facet_col(
    facets = ~group,
    scales = "free_y",
    space = "free"
  )
  colours_BP <- c("black", "gold2", "deepskyblue1", "coral", "springgreen3", "purple", "magenta", "turquoise")
  shapes_BP <- c(23L, 21L, 21L, 21L, 21L, 21L, 21L, 21L)
  x <- x + scale_color_manual(values=colours_BP)
  #x <- x + scale_x_continuous(breaks = c(0.5, 1, 1.5, 2, 2.5, 3, 5))
  x <- x + ggplot2::scale_shape_manual(values = shapes_BP)
  print(x)
}

leave1out_results2 <- leave1out_results2 %>% mutate(A = case_when(
  id.outcome == "5S2uBD"~"F-ALL_gdm_subsamp_FinnGen",
  id.outcome == "FbdyLa"~"F-ALL_gdm_subsamp_MOBA",
  id.outcome == "pnJEcL"~"F-ALL_gdm_subsamp_Public",
  id.outcome == "qBL2hf"~"F-ALL_gh_subsamp_ALSPAC",
  id.outcome == "neKQmh"~"F-ALL_gh_subsamp_BIB-SA",
  id.outcome == "8nnehI"~"F-ALL_gh_subsamp_BIB-WE",
  id.outcome == "v1EoqF"~"F-ALL_gh_subsamp_FinnGen",
  id.outcome == "meqnqp"~"F-ALL_gh_subsamp_MOBA",
  id.outcome == "Cjmesq"~"F-ALL_gh_subsamp_UKB",
  id.outcome == "l82H54"~"F-ALL_hbw_all_ALSPAC",
  id.outcome == "TRAEHU"~"F-ALL_hbw_all_BIB-SA",
  id.outcome == "msaiUC"~"F-ALL_hbw_all_BIB-WE",
  id.outcome == "YtXgpV"~"F-ALL_hbw_all_MOBA",
  id.outcome == "rh2xrL"~"F-ALL_hbw_all_UKB",
  id.outcome == "9d793e"~"F-ALL_hdp_subsamp_ALSPAC",
  id.outcome == "1SpfCW"~"F-ALL_hdp_subsamp_BIB-SA",
  id.outcome == "6UCOQb"~"F-ALL_hdp_subsamp_BIB-WE",
  id.outcome == "lIwwgE"~"F-ALL_hdp_subsamp_FinnGen",
  id.outcome == "Vsr2sF"~"F-ALL_hdp_subsamp_MOBA",
  id.outcome == "V3wmW0"~"F-ALL_hdp_subsamp_UKB",
  id.outcome == "P87HTA"~"F-ALL_lbw_all_ALSPAC",
  id.outcome == "1Q3k9E"~"F-ALL_lbw_all_BIB-SA",
  id.outcome == "oZ3UoH"~"F-ALL_lbw_all_BIB-WE",
  id.outcome == "dIxonP"~"F-ALL_lbw_all_MOBA",
  id.outcome == "VJK524"~"F-ALL_lbw_all_UKB",
  id.outcome == "He1EXm"~"F-ALL_lga_ALSPAC",
  id.outcome == "NdeU8I"~"F-ALL_lga_BIB-SA",
  id.outcome == "5tgpO6"~"F-ALL_lga_BIB-WE",
  id.outcome == "wsTRbo"~"F-ALL_lga_MOBA",
  id.outcome == "3jFInA"~"F-ALL_lga_UKB",
  id.outcome == "n7knSb"~"F-ALL_misc_subsamp_ALSPAC",
  id.outcome == "Its0wc"~"F-ALL_misc_subsamp_BIB-SA",
  id.outcome == "gT5j1L"~"F-ALL_misc_subsamp_BIB-WE",
  id.outcome == "7KjhG7"~"F-ALL_misc_subsamp_FinnGen",
  id.outcome == "IWbgkl"~"F-ALL_misc_subsamp_MOBA",
  id.outcome == "22mCEf"~"F-ALL_misc_subsamp_UKB",
  id.outcome == "RnHVtA"~"F-ALL_pe_subsamp_BIB-SA",
  id.outcome == "btlrKB"~"F-ALL_pe_subsamp_BIB-WE",
  id.outcome == "KjXgAa"~"F-ALL_pe_subsamp_FinnGen",
  id.outcome == "oT398E"~"F-ALL_pe_subsamp_MOBA",
  id.outcome == "ZzEy0G"~"F-ALL_pe_subsamp_UKB",
  id.outcome == "zN7W3I"~"F-ALL_pe_subsamp_Public",
  id.outcome == "joW66G"~"F-ALL_pretb_all_BIB-SA",
  id.outcome == "1xcl5N"~"F-ALL_pretb_all_BIB-WE",
  id.outcome == "nmfDES"~"F-ALL_pretb_all_FinnGen",
  id.outcome == "ZQ44hw"~"F-ALL_pretb_all_MOBA",
  id.outcome == "HK2mFV"~"F-ALL_pretb_all_UKB",
  id.outcome == "iQi4bb"~"F-ALL_pretb_all_Public",
  id.outcome == "VSY02P"~"F-ALL_sb_subsamp_ALSPAC",
  id.outcome == "hUegd2"~"F-ALL_sb_subsamp_BIB-SA",
  id.outcome == "hHFTA2"~"F-ALL_sb_subsamp_BIB-WE",
  id.outcome == "5ksrAj"~"F-ALL_sb_subsamp_MOBA",
  id.outcome == "peITPn"~"F-ALL_sb_subsamp_UKB",
  id.outcome == "S1xMVE"~"F-ALL_sga_ALSPAC",
  id.outcome == "7AFzlm"~"F-ALL_sga_BIB-SA",
  id.outcome == "ZiT8mz"~"F-ALL_sga_BIB-WE",
  id.outcome == "ygJY5j"~"F-ALL_sga_MOBA",
  id.outcome == "sKxpyr"~"F-ALL_sga_UKB",
  id.outcome == "bL4bKv"~"F-ALL_vpretb_all_ALSPAC",
  id.outcome == "v35LGW"~"F-ALL_vpretb_all_BIB-SA",
  id.outcome == "kWedKq"~"F-ALL_vpretb_all_BIB-WE",
  id.outcome == "9e5nIy"~"F-ALL_vpretb_all_MOBA",
  id.outcome == "fQPe9Y"~"F-ALL_vpretb_all_UKB",
  id.outcome == "p2zSZ7"~"F-INCL_gdm_subsamp_FinnGen",
  id.outcome == "F0AqBC"~"F-INCL_gdm_subsamp_MOBA",
  id.outcome == "Pazs6j"~"F-INCL_gdm_subsamp_Public",
  id.outcome == "SSuDR1"~"F-INCL_gh_subsamp_ALSPAC",
  id.outcome == "dX9lPf"~"F-INCL_gh_subsamp_BIB-SA",
  id.outcome == "te95tW"~"F-INCL_gh_subsamp_BIB-WE",
  id.outcome == "B6AZif"~"F-INCL_gh_subsamp_FinnGen",
  id.outcome == "kaoyli"~"F-INCL_gh_subsamp_MOBA",
  id.outcome == "VNeHM7"~"F-INCL_gh_subsamp_UKB",
  id.outcome == "xrwG7H"~"F-INCL_hbw_all_ALSPAC",
  id.outcome == "3U5X8P"~"F-INCL_hbw_all_BIB-SA",
  id.outcome == "iZuCrC"~"F-INCL_hbw_all_BIB-WE",
  id.outcome == "dMJYBB"~"F-INCL_hbw_all_MOBA",
  id.outcome == "FZM1Hy"~"F-INCL_hbw_all_UKB",
  id.outcome == "vFp21M"~"F-INCL_hdp_subsamp_ALSPAC",
  id.outcome == "hY10vn"~"F-INCL_hdp_subsamp_BIB-SA",
  id.outcome == "DpTvyt"~"F-INCL_hdp_subsamp_BIB-WE",
  id.outcome == "Juptg3"~"F-INCL_hdp_subsamp_FinnGen",
  id.outcome == "Vi56UN"~"F-INCL_hdp_subsamp_MOBA",
  id.outcome == "rDgLTE"~"F-INCL_hdp_subsamp_UKB",
  id.outcome == "rrVb73"~"F-INCL_lbw_all_ALSPAC",
  id.outcome == "71CDWY"~"F-INCL_lbw_all_BIB-SA",
  id.outcome == "s00kcH"~"F-INCL_lbw_all_BIB-WE",
  id.outcome == "EsQR1J"~"F-INCL_lbw_all_MOBA",
  id.outcome == "FVPFlr"~"F-INCL_lbw_all_UKB",
  id.outcome == "2Tq2Rh"~"F-INCL_lga_ALSPAC",
  id.outcome == "MBku59"~"F-INCL_lga_BIB-SA",
  id.outcome == "bU1pGG"~"F-INCL_lga_BIB-WE",
  id.outcome == "EFlIKQ"~"F-INCL_lga_MOBA",
  id.outcome == "n9wDBE"~"F-INCL_lga_UKB",
  id.outcome == "Vu1YcI"~"F-INCL_misc_subsamp_ALSPAC",
  id.outcome == "Jf8Tpe"~"F-INCL_misc_subsamp_BIB-SA",
  id.outcome == "gfRByE"~"F-INCL_misc_subsamp_BIB-WE",
  id.outcome == "DlLOj7"~"F-INCL_misc_subsamp_FinnGen",
  id.outcome == "ryjG5q"~"F-INCL_misc_subsamp_MOBA",
  id.outcome == "8u7Stz"~"F-INCL_misc_subsamp_UKB",
  id.outcome == "yAcq8h"~"F-INCL_pe_subsamp_BIB-SA",
  id.outcome == "AmhlcD"~"F-INCL_pe_subsamp_BIB-WE",
  id.outcome == "xRoJtt"~"F-INCL_pe_subsamp_FinnGen",
  id.outcome == "914EkG"~"F-INCL_pe_subsamp_MOBA",
  id.outcome == "1l9xeN"~"F-INCL_pe_subsamp_UKB",
  id.outcome == "taHFMy"~"F-INCL_pe_subsamp_Public",
  id.outcome == "xAR6eo"~"F-INCL_pretb_all_BIB-SA",
  id.outcome == "wAVqYA"~"F-INCL_pretb_all_BIB-WE",
  id.outcome == "Ogt9pm"~"F-INCL_pretb_all_FinnGen",
  id.outcome == "3n9z8K"~"F-INCL_pretb_all_MOBA",
  id.outcome == "2eErAt"~"F-INCL_pretb_all_UKB",
  id.outcome == "QjcN3P"~"F-INCL_pretb_all_Public",
  id.outcome == "iN50Xj"~"F-INCL_sb_subsamp_ALSPAC",
  id.outcome == "NYLjlK"~"F-INCL_sb_subsamp_BIB-SA",
  id.outcome == "6izIDs"~"F-INCL_sb_subsamp_BIB-WE",
  id.outcome == "tcVlzQ"~"F-INCL_sb_subsamp_MOBA",
  id.outcome == "11itMc"~"F-INCL_sb_subsamp_UKB",
  id.outcome == "RgsoTM"~"F-INCL_sga_ALSPAC",
  id.outcome == "ULNGg4"~"F-INCL_sga_BIB-SA",
  id.outcome == "GlP7q0"~"F-INCL_sga_BIB-WE",
  id.outcome == "NoHZSE"~"F-INCL_sga_MOBA",
  id.outcome == "LxV8OB"~"F-INCL_sga_UKB",
  id.outcome == "790Bsj"~"F-INCL_vpretb_all_ALSPAC",
  id.outcome == "KVENIk"~"F-INCL_vpretb_all_BIB-SA",
  id.outcome == "OcTPSE"~"F-INCL_vpretb_all_BIB-WE",
  id.outcome == "P8iuly"~"F-INCL_vpretb_all_MOBA",
  id.outcome == "Y4Z9ZQ"~"F-INCL_vpretb_all_UKB",
  id.outcome == "0kTw0P"~"F-ANOV_gdm_subsamp_FinnGen",
  id.outcome == "oRNyLx"~"F-ANOV_gdm_subsamp_MOBA",
  id.outcome == "GyApHL"~"F-ANOV_gdm_subsamp_Public",
  id.outcome == "3e8ofU"~"F-ANOV_gh_subsamp_ALSPAC",
  id.outcome == "aZie4Q"~"F-ANOV_gh_subsamp_BIB-SA",
  id.outcome == "rwNXda"~"F-ANOV_gh_subsamp_BIB-WE",
  id.outcome == "S1IMkg"~"F-ANOV_gh_subsamp_FinnGen",
  id.outcome == "L7JZr4"~"F-ANOV_gh_subsamp_MOBA",
  id.outcome == "EiZSdB"~"F-ANOV_gh_subsamp_UKB",
  id.outcome == "TZgxrp"~"F-ANOV_hbw_all_ALSPAC",
  id.outcome == "8yBXey"~"F-ANOV_hbw_all_BIB-SA",
  id.outcome == "d76SfP"~"F-ANOV_hbw_all_BIB-WE",
  id.outcome == "Dh1rLQ"~"F-ANOV_hbw_all_MOBA",
  id.outcome == "5i76E4"~"F-ANOV_hbw_all_UKB",
  id.outcome == "GfxBH7"~"F-ANOV_hdp_subsamp_ALSPAC",
  id.outcome == "oPPKq9"~"F-ANOV_hdp_subsamp_BIB-SA",
  id.outcome == "SFHsIa"~"F-ANOV_hdp_subsamp_BIB-WE",
  id.outcome == "OUpAhn"~"F-ANOV_hdp_subsamp_FinnGen",
  id.outcome == "NqNswn"~"F-ANOV_hdp_subsamp_MOBA",
  id.outcome == "fxdKwN"~"F-ANOV_hdp_subsamp_UKB",
  id.outcome == "YwAPy0"~"F-ANOV_lbw_all_ALSPAC",
  id.outcome == "5ckhuv"~"F-ANOV_lbw_all_BIB-SA",
  id.outcome == "X9T5oM"~"F-ANOV_lbw_all_BIB-WE",
  id.outcome == "3yvzUW"~"F-ANOV_lbw_all_MOBA",
  id.outcome == "9Tn9Qq"~"F-ANOV_lbw_all_UKB",
  id.outcome == "rd8u4m"~"F-ANOV_lga_ALSPAC",
  id.outcome == "gRoPZK"~"F-ANOV_lga_BIB-SA",
  id.outcome == "wqNE2I"~"F-ANOV_lga_BIB-WE",
  id.outcome == "9ftG0a"~"F-ANOV_lga_MOBA",
  id.outcome == "p2qk9n"~"F-ANOV_lga_UKB",
  id.outcome == "9WxvB5"~"F-ANOV_misc_subsamp_ALSPAC",
  id.outcome == "UE08Kh"~"F-ANOV_misc_subsamp_BIB-SA",
  id.outcome == "OLnZBL"~"F-ANOV_misc_subsamp_BIB-WE",
  id.outcome == "9BPpDf"~"F-ANOV_misc_subsamp_FinnGen",
  id.outcome == "4xRz5U"~"F-ANOV_misc_subsamp_MOBA",
  id.outcome == "Bu1cYe"~"F-ANOV_misc_subsamp_UKB",
  id.outcome == "ifnAMR"~"F-ANOV_pe_subsamp_BIB-SA",
  id.outcome == "dNArNP"~"F-ANOV_pe_subsamp_BIB-WE",
  id.outcome == "nytb2U"~"F-ANOV_pe_subsamp_FinnGen",
  id.outcome == "GpDYFW"~"F-ANOV_pe_subsamp_MOBA",
  id.outcome == "vC880p"~"F-ANOV_pe_subsamp_UKB",
  id.outcome == "PqxCaV"~"F-ANOV_pe_subsamp_Public",
  id.outcome == "S9migX"~"F-ANOV_pretb_all_BIB-SA",
  id.outcome == "2hV01E"~"F-ANOV_pretb_all_BIB-WE",
  id.outcome == "PfZ7Vb"~"F-ANOV_pretb_all_FinnGen",
  id.outcome == "EVFxuy"~"F-ANOV_pretb_all_MOBA",
  id.outcome == "j16pVi"~"F-ANOV_pretb_all_UKB",
  id.outcome == "LKxyBk"~"F-ANOV_pretb_all_Public",
  id.outcome == "b71TGV"~"F-ANOV_sb_subsamp_ALSPAC",
  id.outcome == "eNJn2U"~"F-ANOV_sb_subsamp_BIB-SA",
  id.outcome == "xayBrz"~"F-ANOV_sb_subsamp_BIB-WE",
  id.outcome == "EfTGhz"~"F-ANOV_sb_subsamp_MOBA",
  id.outcome == "b0JODB"~"F-ANOV_sb_subsamp_UKB",
  id.outcome == "1zgbEH"~"F-ANOV_sga_ALSPAC",
  id.outcome == "py8b3j"~"F-ANOV_sga_BIB-SA",
  id.outcome == "K78WS6"~"F-ANOV_sga_BIB-WE",
  id.outcome == "pt2RCQ"~"F-ANOV_sga_MOBA",
  id.outcome == "oCTwJ2"~"F-ANOV_sga_UKB",
  id.outcome == "DVWDo1"~"F-ANOV_vpretb_all_ALSPAC",
  id.outcome == "EFaJ7r"~"F-ANOV_vpretb_all_BIB-SA",
  id.outcome == "suCSr9"~"F-ANOV_vpretb_all_BIB-WE",
  id.outcome == "Voa8c5"~"F-ANOV_vpretb_all_MOBA",
  id.outcome == "XIC2Xn"~"F-ANOV_vpretb_all_UKB",
  id.outcome == "xLIwYt"~"F-EXCL_gdm_subsamp_FinnGen",
  id.outcome == "qXt9hF"~"F-EXCL_gdm_subsamp_MOBA",
  id.outcome == "ddkgoK"~"F-EXCL_gdm_subsamp_Public",
  id.outcome == "WrktGP"~"F-EXCL_gh_subsamp_ALSPAC",
  id.outcome == "YhPR3E"~"F-EXCL_gh_subsamp_BIB-SA",
  id.outcome == "cn3dVE"~"F-EXCL_gh_subsamp_BIB-WE",
  id.outcome == "aXKcgT"~"F-EXCL_gh_subsamp_FinnGen",
  id.outcome == "0F3M7R"~"F-EXCL_gh_subsamp_MOBA",
  id.outcome == "aVLhJR"~"F-EXCL_gh_subsamp_UKB",
  id.outcome == "ssaUwT"~"F-EXCL_hbw_all_ALSPAC",
  id.outcome == "yukcJG"~"F-EXCL_hbw_all_BIB-SA",
  id.outcome == "3KjSLw"~"F-EXCL_hbw_all_BIB-WE",
  id.outcome == "0OPlms"~"F-EXCL_hbw_all_MOBA",
  id.outcome == "va7BlO"~"F-EXCL_hbw_all_UKB",
  id.outcome == "s4fhDD"~"F-EXCL_hdp_subsamp_ALSPAC",
  id.outcome == "Pkurtg"~"F-EXCL_hdp_subsamp_BIB-SA",
  id.outcome == "a2Je0P"~"F-EXCL_hdp_subsamp_BIB-WE",
  id.outcome == "wIvnFd"~"F-EXCL_hdp_subsamp_FinnGen",
  id.outcome == "FnVDNm"~"F-EXCL_hdp_subsamp_MOBA",
  id.outcome == "HjCMlF"~"F-EXCL_hdp_subsamp_UKB",
  id.outcome == "fj2cZf"~"F-EXCL_lbw_all_ALSPAC",
  id.outcome == "wPksvV"~"F-EXCL_lbw_all_BIB-SA",
  id.outcome == "GfTXOT"~"F-EXCL_lbw_all_BIB-WE",
  id.outcome == "6Zw0vT"~"F-EXCL_lbw_all_MOBA",
  id.outcome == "WT0431"~"F-EXCL_lbw_all_UKB",
  id.outcome == "54YLnA"~"F-EXCL_lga_ALSPAC",
  id.outcome == "sCaSmm"~"F-EXCL_lga_BIB-SA",
  id.outcome == "bL9G2F"~"F-EXCL_lga_BIB-WE",
  id.outcome == "4r2Dk8"~"F-EXCL_lga_MOBA",
  id.outcome == "CZzi3f"~"F-EXCL_lga_UKB",
  id.outcome == "aJhmDL"~"F-EXCL_misc_subsamp_ALSPAC",
  id.outcome == "EZItFG"~"F-EXCL_misc_subsamp_BIB-SA",
  id.outcome == "739KLy"~"F-EXCL_misc_subsamp_BIB-WE",
  id.outcome == "Cvc19A"~"F-EXCL_misc_subsamp_FinnGen",
  id.outcome == "75LGDr"~"F-EXCL_misc_subsamp_MOBA",
  id.outcome == "CQYx1t"~"F-EXCL_misc_subsamp_UKB",
  id.outcome == "w4Tmei"~"F-EXCL_pe_subsamp_BIB-SA",
  id.outcome == "Z2vVgL"~"F-EXCL_pe_subsamp_BIB-WE",
  id.outcome == "Z4qFob"~"F-EXCL_pe_subsamp_FinnGen",
  id.outcome == "LM1QH6"~"F-EXCL_pe_subsamp_MOBA",
  id.outcome == "OWwVF5"~"F-EXCL_pe_subsamp_UKB",
  id.outcome == "zgfY4v"~"F-EXCL_pe_subsamp_Public",
  id.outcome == "4u4oxY"~"F-EXCL_pretb_all_BIB-SA",
  id.outcome == "2IEknM"~"F-EXCL_pretb_all_BIB-WE",
  id.outcome == "0yaFbY"~"F-EXCL_pretb_all_FinnGen",
  id.outcome == "ofMXs9"~"F-EXCL_pretb_all_MOBA",
  id.outcome == "kJHuHj"~"F-EXCL_pretb_all_UKB",
  id.outcome == "cS2I0O"~"F-EXCL_pretb_all_Public",
  id.outcome == "4J3WPv"~"F-EXCL_sb_subsamp_ALSPAC",
  id.outcome == "HL2VHE"~"F-EXCL_sb_subsamp_BIB-SA",
  id.outcome == "fui4OX"~"F-EXCL_sb_subsamp_BIB-WE",
  id.outcome == "4xkX0m"~"F-EXCL_sb_subsamp_MOBA",
  id.outcome == "u0HhLr"~"F-EXCL_sb_subsamp_UKB",
  id.outcome == "nh6Ra8"~"F-EXCL_sga_ALSPAC",
  id.outcome == "EZebN9"~"F-EXCL_sga_BIB-SA",
  id.outcome == "BJAqLi"~"F-EXCL_sga_BIB-WE",
  id.outcome == "j9EwhC"~"F-EXCL_sga_MOBA",
  id.outcome == "WThTLx"~"F-EXCL_sga_UKB",
  id.outcome == "hJPrmL"~"F-EXCL_vpretb_all_ALSPAC",
  id.outcome == "fEPQR1"~"F-EXCL_vpretb_all_BIB-SA",
  id.outcome == "Vcn6uI"~"F-EXCL_vpretb_all_BIB-WE",
  id.outcome == "tv03FD"~"F-EXCL_vpretb_all_MOBA",
  id.outcome == "Ujuu9g"~"F-EXCL_vpretb_all_UKB",
  id.outcome == "o05X0R"~"F-ANAT_gdm_subsamp_FinnGen",
  id.outcome == "zngol7"~"F-ANAT_gdm_subsamp_MOBA",
  id.outcome == "Pq2TW9"~"F-ANAT_gdm_subsamp_Public",
  id.outcome == "uQUX3G"~"F-ANAT_gh_subsamp_ALSPAC",
  id.outcome == "LAxy3Z"~"F-ANAT_gh_subsamp_BIB-SA",
  id.outcome == "GJbOsc"~"F-ANAT_gh_subsamp_BIB-WE",
  id.outcome == "pVA3SS"~"F-ANAT_gh_subsamp_FinnGen",
  id.outcome == "2hfISW"~"F-ANAT_gh_subsamp_MOBA",
  id.outcome == "kjLxyL"~"F-ANAT_gh_subsamp_UKB",
  id.outcome == "LnhM33"~"F-ANAT_hbw_all_ALSPAC",
  id.outcome == "9ml7E3"~"F-ANAT_hbw_all_BIB-SA",
  id.outcome == "ocKyBy"~"F-ANAT_hbw_all_BIB-WE",
  id.outcome == "SeJRZS"~"F-ANAT_hbw_all_MOBA",
  id.outcome == "sH9fZ8"~"F-ANAT_hbw_all_UKB",
  id.outcome == "W3j95p"~"F-ANAT_hdp_subsamp_ALSPAC",
  id.outcome == "wZ9L6t"~"F-ANAT_hdp_subsamp_BIB-SA",
  id.outcome == "RWBYi3"~"F-ANAT_hdp_subsamp_BIB-WE",
  id.outcome == "6samSN"~"F-ANAT_hdp_subsamp_FinnGen",
  id.outcome == "QiyonP"~"F-ANAT_hdp_subsamp_MOBA",
  id.outcome == "9UejSz"~"F-ANAT_hdp_subsamp_UKB",
  id.outcome == "xGj5JQ"~"F-ANAT_lbw_all_ALSPAC",
  id.outcome == "R5ZJGv"~"F-ANAT_lbw_all_BIB-SA",
  id.outcome == "Nos5Mj"~"F-ANAT_lbw_all_BIB-WE",
  id.outcome == "DhNJr2"~"F-ANAT_lbw_all_MOBA",
  id.outcome == "uOPHvY"~"F-ANAT_lbw_all_UKB",
  id.outcome == "OWVp8l"~"F-ANAT_lga_ALSPAC",
  id.outcome == "KCmuw2"~"F-ANAT_lga_BIB-SA",
  id.outcome == "uXtxtA"~"F-ANAT_lga_BIB-WE",
  id.outcome == "5Qbxxd"~"F-ANAT_lga_MOBA",
  id.outcome == "KkDCc1"~"F-ANAT_lga_UKB",
  id.outcome == "jTbMV1"~"F-ANAT_misc_subsamp_ALSPAC",
  id.outcome == "trSm8e"~"F-ANAT_misc_subsamp_BIB-SA",
  id.outcome == "VO8bWJ"~"F-ANAT_misc_subsamp_BIB-WE",
  id.outcome == "CX1qao"~"F-ANAT_misc_subsamp_FinnGen",
  id.outcome == "ARWdbP"~"F-ANAT_misc_subsamp_MOBA",
  id.outcome == "XIBCqh"~"F-ANAT_misc_subsamp_UKB",
  id.outcome == "ZJK3D0"~"F-ANAT_pe_subsamp_BIB-SA",
  id.outcome == "IndHdf"~"F-ANAT_pe_subsamp_BIB-WE",
  id.outcome == "BOJ9CG"~"F-ANAT_pe_subsamp_FinnGen",
  id.outcome == "1llAMa"~"F-ANAT_pe_subsamp_MOBA",
  id.outcome == "qZZLUH"~"F-ANAT_pe_subsamp_UKB",
  id.outcome == "4Y4O8N"~"F-ANAT_pe_subsamp_Public",
  id.outcome == "PyTTYI"~"F-ANAT_pretb_all_BIB-SA",
  id.outcome == "YkNEZ7"~"F-ANAT_pretb_all_BIB-WE",
  id.outcome == "H3C10l"~"F-ANAT_pretb_all_FinnGen",
  id.outcome == "kytaWO"~"F-ANAT_pretb_all_MOBA",
  id.outcome == "U4V7ra"~"F-ANAT_pretb_all_UKB",
  id.outcome == "Z2CFqS"~"F-ANAT_pretb_all_Public",
  id.outcome == "vEqPuo"~"F-ANAT_sb_subsamp_ALSPAC",
  id.outcome == "TkC5Lh"~"F-ANAT_sb_subsamp_BIB-SA",
  id.outcome == "sPWVID"~"F-ANAT_sb_subsamp_BIB-WE",
  id.outcome == "lZyeEo"~"F-ANAT_sb_subsamp_MOBA",
  id.outcome == "U4ccfc"~"F-ANAT_sb_subsamp_UKB",
  id.outcome == "052OWG"~"F-ANAT_sga_ALSPAC",
  id.outcome == "vSJgKF"~"F-ANAT_sga_BIB-SA",
  id.outcome == "k2OGPt"~"F-ANAT_sga_BIB-WE",
  id.outcome == "KIPaJF"~"F-ANAT_sga_MOBA",
  id.outcome == "NVEimD"~"F-ANAT_sga_UKB",
  id.outcome == "NYANUj"~"F-ANAT_vpretb_all_ALSPAC",
  id.outcome == "wXvsHu"~"F-ANAT_vpretb_all_BIB-SA",
  id.outcome == "M8eBZj"~"F-ANAT_vpretb_all_BIB-WE",
  id.outcome == "JYyqLm"~"F-ANAT_vpretb_all_MOBA",
  id.outcome == "IlNskX"~"F-ANAT_vpretb_all_UKB")
)

for (exp_name in names(exposure_list)) {
  
  # Convert exposure
  exposure_dat <- format_data(
    exposure_list[[exp_name]],
    type = "exposure",
    snp_col = "RSID",
    beta_col = "beta",
    se_col = "se",
    pval_col = "pval",
    eaf_col = "EAF",
    effect_allele_col = "effect_allele",
    other_allele_col = "other_allele"
  )
  exposure_dat$phenotype <- exp_name
  
  # Loop outcome
  for (outcome_name in names(outcome_list_trios)) {
    
    outcome_data_raw <- outcome_list_trios[[outcome_name]]
    
    # Convert outcome
    outcome_dat <- format_data(
      outcome_data_raw,
      type = "outcome",
      snp_col = "SNP",
      beta_col = "beta_mat_donuts",
      se_col = "se_mat_donuts",
      pval_col = "p_mat_donuts",
      eaf_col = "eaf_mat",
      effect_allele_col = "effect_allele",
      other_allele_col = "other_allele",
      phenotype = outcome_name
    )
    
    # harmonise
    harmonised <- harmonise_data(exposure_dat, outcome_dat)
    harmonised$exposure <- exp_name
    harmonised$outcome <- outcome_name
    harmonised_list_mat_donut_trios[[paste0(exp_name, "_", outcome_name)]] <- harmonised
    
    if (nrow(harmonised) == 0) next
    
    # MR 
    mr_result <- mr(harmonised)
    mr_result$exposure <- exp_name
    mr_result$outcome <- outcome_name
    
    all_results_mat_donut_trios[[paste0(exp_name, "_", outcome_name)]] <- mr_result
  }
}

#########

glyca_comb_hits <- rename(glyca_comb_hits, c( 
  "effect_allele"="ALLELE1",
  "other_allele"="ALLELE0",
  "beta"="BETA",
  "se"="SE",
  "pval"="P_BOLT_LMM",
  "eaf" = "A1FREQ",
  "chr" = "CHR",
  "pos" = "GENPOS"))


#####

out_func <- function(study_name)
{
  # Extract outcome SNPs matching the SNPs in the exposure dataset
  outcome_data <- outcome_variables %>% 
    filter(., study==study_name) %>%
    mutate("id" = paste(study, Phenotype, sep = "_")) %>% 
    data.frame()
  
  outcome_data <- format_data(outcome_data, snps = exp_dat$SNP, type = "outcome")
  
  return(outcome_data)
}

######

subset_df2 <- function(df, exp, meth = "Inverse variance weighted", 
                      stu = "Metanalysis",  typ = "binary", pri = "primary") {
  
  df <- filter(df, 
               id.exposure %in% exp,
               method %in% meth,
               study == stu,
               type == typ,
               priority == pri
  ) 
}

prep_func2 <- function(results_file) {
  
  #rename outcomes
  results_file <- rename(results_file, out_var = outcome) %>%
    merge(., out_info, by.x = "out_var", by.y = "Variable_name")
  
  #reorder methods
  as.factor(results_file$method)
  results_file$method <- factor(results_file$method, levels = c("Weighted mode", "Weighted median",
                                                                "MR Egger", "Inverse variance weighted"))
  
  #reorder outcome group
  as.factor(results_file$group)
  results_file$group <- factor(results_file$group, levels = c("Pregnancy loss outcomes", "Maternal morbidity outcomes", "Labour outcomes",
                                                              "Offspring birth outcomes", "Continuous outcomes"))
  
  # harmonise exposure names between Keaton and Ehret
  results_file <- mutate(results_file, id.exposure = case_when(
    id.exposure == "ebi-a-GCST90012005" ~ "IL6 (Folkerson et al) - IRNT",
    id.exposure == "ebi-a-GCST90012025" ~ "IL6-R (Folkerson et al) - IRNT",
    id.exposure == "ebi-a-GCST90029070" ~ "CRP (Said et al) - SD",
    id.exposure == "IL6" ~ "IL6 (Ahluwalia et al) - SD",
    id.exposure == "GlycA" ~ "GlycA (UKB) - SD",
    id.exposure == "GlycA-IRNT" ~ "GlycA (UKB) - IRNT",
    .default = exposure
  ))
  return(results_file)
  
}

#####

out_func3 <- function(study_name)
{
  # Extract outcome SNPs matching the SNPs in the exposure dataset
  outcome_data <- outcome_variables %>% 
    filter(., study==study_name) %>%
    mutate("id" = paste(study, Phenotype, sep = "_")) %>% 
    data.frame()
  
  outcome_data <- format_data(outcome_data, snps = cis_mr_dat$SNP, type = "outcome")
  
  return(outcome_data)
}

prep_func3 <- function(results_file) {
  
  #rename outcomes
  results_file <- rename(results_file, out_var = outcome) %>%
    merge(., out_info, by.x = "out_var", by.y = "Variable_name")
  
  #reorder methods
  as.factor(results_file$method)
  results_file$method <- factor(results_file$method, levels = c("Weighted mode", "Weighted median",
                                                                "MR Egger", "Inverse variance weighted"))
  
  #reorder outcome group
  as.factor(results_file$group)
  results_file$group <- factor(results_file$group, levels = c("Pregnancy loss outcomes", "Maternal morbidity outcomes", "Labour outcomes",
                                                              "Offspring birth outcomes", "Continuous outcomes"))
  
  # harmonise exposure names between Keaton and Ehret
  results_file <- mutate(results_file, id.exposure = case_when(
    id.exposure == "ebi-a-GCST90012025" ~ "IL6-R",
    id.exposure == "ebi-a-GCST90029070" ~ "CRP",
    .default = exposure
  ))
  return(results_file)
  
}

subset_df3 <- function(df, exp, meth, typ = "binary", pri = "primary") {
  
  df <- filter(df, 
               id.exposure %in% exp,
               method %in% meth,
               type == typ,
               priority == pri
  ) 
}


myforestplot2 <- function(dt, log_T = TRUE, exp_type, limits, effect = "Odds ratio", one_colour = T) {
  
  if("p.fdr" %in% names(dt)) {
    dt <- mutate(dt, p_val = p.fdr)
  } else {
    dt <- mutate(dt, p_val = pval)
  }
  
  dt <- arrange(dt, order) 
  
  xlabel <- paste(effect, "(95% CI) per unit higher", exp_type)
  
  if(one_colour == T) {
    p <- dt %>%  
      forestplot(
        df = .,
        estimate = b,    
        se = se,
        pvalue = p_val,
        name = Outcome,
        logodds = log_T,
        psignif = 0.05,
        colour = NULL,
        xlab = xlabel,
        xlim = limits
      ) 
  } else {  
    p <- dt %>%  
      forestplot(
        df = .,
        estimate = b,   
        se = se,
        pvalue = p_val,
        name = Outcome,
        logodds = log_T,
        psignif = 1,
        colour = var_colour,
        xlab = xlabel,
        xlim = limits
      ) +
      theme(legend.title = element_blank())
  }   
  
  p <- p + 
    theme(axis.title.x = element_text(color="black", size=10, face="bold")) +
    ggforce::facet_col(
      facets = ~group,
      scales = "free_y",
      space = "free"
    ) 
  
  return(p)
}

#####

out_func4 <- function(study_name)
{
  # Extract outcome SNPs matching the SNPs in the exposure dataset
  outcome_data <- outcome_variables %>% 
    filter(., study==study_name) %>%
    mutate("id" = paste(study, Phenotype, sep = "_")) %>% 
    data.frame()
  
  outcome_data <- format_data(outcome_data, snps = exp_dat_crp_2$SNP, type = "outcome")
  
  return(outcome_data)
}

out_func <- function(study_name)
{
  # Extract outcome SNPs matching the SNPs in the exposure dataset
  outcome_data <- outcome_variables %>% 
    filter(., study==study_name) %>%
    mutate("id" = paste(study, Phenotype, sep = "_")) %>% 
    data.frame()
  
  outcome_data <- format_data(outcome_data, snps = CRP_glucose_exp$SNP, type = "outcome")
  
  return(outcome_data)
}
CRP_glucose_gdm_mv <- mv_harmonise_data(CRP_glucose_exp,out_dat_gdm,harmonise_strictness = 2)

#######

myforestplot3 <- function(dt, log_T = TRUE, exp_type, limits, effect = "Odds ratio", one_colour = T) {
  
  if("p.fdr" %in% names(dt)) {
    dt <- mutate(dt, p_val = p.fdr)
  } else {
    dt <- mutate(dt, p_val = pval)
  }
  
  dt <- arrange(dt, order) 
  
  xlabel <- paste(effect, "(95% CI) per unit higher", exp_type)
  
  if(one_colour == T) {
    p <- dt %>%
      forestplot(
        estimate = b,    
        se = se,
        pvalue = p_val,
        name = Outcome,
        logodds = log_T,
        psignif = 1,
        colour = "constant",
        xlab = xlabel,
        xlim = limits
      )+scale_colour_manual(values = c("constant" = "#8F87FF"))+guides(colour = "none")
 } else {  
    p <- dt %>%  
      forestplot(
        estimate = b,    
        se = se,
        pvalue = p_val,
        name = Outcome,
        logodds = log_T,
        psignif = 1,
        colour = var_colour,
        xlab = xlabel,
        xlim = limits
      ) +
      theme(legend.title = element_blank())
  }   
  
  p <- p + 
    theme(axis.title.x = element_text(color="black", size=10, face="bold")) +
    ggforce::facet_col(
      facets = ~group,
      scales = "free_y",
      space = "free"
    ) 
  
  return(p)
}

####

out_func2 <- function(study_name)
{
  # Extract outcome SNPs matching the SNPs in the exposure dataset
  outcome_data <- outcome_variables %>% 
    filter(., study==study_name) %>%
    mutate("id" = paste(study, Phenotype, sep = "_")) %>% 
    data.frame()
  
  outcome_data <- format_data(outcome_data, snps = All_proxy_dat$rsid, type = "outcome")
  
  return(outcome_data)
}

#####

generate_outdat_with_proxies <- function(exposure_dat, outcome_dat, outcome_name,
                                         proxies){
  
  ### Set up
  i_outcome <- outcome_name
  # Filter specific outcome dataset to available exposure SNPs
  tmp_outcome_dat <- outcome_dat %>%
    filter(SNP %in% exposure_dat$SNP, outcome == i_outcome)
  
  ### Remove any duplicated SNPs
  tmp_outcome_dat <- tmp_outcome_dat[!duplicated(tmp_outcome_dat$SNP), ]
  
  ### Check whether outcome is available
  if(dim(tmp_outcome_dat)[1] == 0){
    print(paste("Outcome", outcome_name, "not available"))
    
    # Return data with 0 rows if outcome was not available
    tmp_outcome_dat
    
    ### Check whether proxies needed
  } else if (dim(tmp_outcome_dat)[1] == dim(exposure_dat)[1]){
    print(paste("No proxies needed for", outcome_name))
    
    # Return full outcome dataset if no proxies needed
    tmp_outcome_dat
    
  } else {
    ### Identify proxies
    
    # Specify all available SNPs
    outcome_snps <- expand.grid(
      SNP = c(exposure_dat$SNP), outcome = unique(outcome_dat$outcome))
    # Left join outcome_dat, NAs if missing
    outcome_snps_available <-
      left_join(outcome_snps, outcome_dat, by = c("SNP", "outcome")) %>%
      # Add indicator column for NAs
      mutate(missing = ifelse(is.na(effect_allele.outcome),
                              TRUE, FALSE))
    # Deduplicated list of all SNPs which we need a proxy for
    need_proxies <- outcome_snps_available %>%
      filter(missing == TRUE) %>%
      distinct(SNP)
    # Count of proxies needed by outcome
    table(outcome_snps_available$outcome[outcome_snps_available$missing == TRUE])
    # List of all proxies
    proxy_snps <- unique(proxies$rsid)
    # Proxies without query SNPs
    proxies_merge_diff <- subset(proxies, query_rsid!=rsid)
    proxies_merge_diff$SNP<-proxies_merge_diff$rsid
    
    outcome_dat_proxies_tmp <-subset(outcome_dat,
                                     outcome_dat$outcome== i_outcome)
    
    # All SNPs we need to proxy for this outcome
    outcome_snps_tmp <- subset(outcome_snps, outcome_snps$outcome==i_outcome)
    # Available SNPs in outcome data:
    outcome_snps_available <-
      left_join(outcome_snps_tmp, outcome_dat, by = c("SNP", "outcome")) %>%
      # Add indicator column for NAs
      mutate(missing = ifelse(is.na(effect_allele.outcome),
                              TRUE, FALSE))
    # SNPs not available for our outcome in outcome data, for which we need proxies:
    need_proxies <- outcome_snps_available %>%
      filter(missing == TRUE) %>%
      distinct(SNP)
    
    if(dim(need_proxies)[1] > 0){
      
      print(paste0("Searching for proxies for ", dim(need_proxies)[1], " SNPs."))
      
      # Append proxies information -
      # For every SNP needing a proxy ('SNP'),
      # create row with information about the proxy ('rsid' & 'SNP.y' columns)
      tmp1 <- merge(x=need_proxies, y=proxies_merge_diff, by.x="SNP", by.y="query_rsid", all.x=TRUE)
      
      # For each proxy SNP ('rsid'), append the extracted proxy-outcome GWAS data
      tmp2 <-merge(x=tmp1, y=outcome_dat_proxies_tmp, by.x="rsid", by.y="SNP", all.y=TRUE)
      # this leaves NAs where any proxy SNPs could not be extracted
      # since they were not available in outcome data
      
      # Exclude palindromic proxies:
      tmp2 <-subset(tmp2, !(effect_allele.outcome == "A" & other_allele.outcome == "T" |
                              effect_allele.outcome == "T" & other_allele.outcome == "A") )
      tmp2 <-subset(tmp2, !(effect_allele.outcome == "C" & other_allele.outcome == "G" |
                              effect_allele.outcome == "G" & other_allele.outcome == "C") )
      
      # Select the proxy SNP ('rsid') in highest LD with the query SNP ('SNP') needing a proxy -
      set.seed(67898) # to ensure slice_sample() is reproducible
      tmp3 <- tmp2 %>%
        group_by(SNP) %>%
        # Remove any proxies with missing outcome data
        drop_na() %>%
        # Select max R2 value, keep ties if several have same R2
        slice_max(R2, with_ties = TRUE) %>%
        # Select random proxy SNP if several have same R2
        slice_sample(n = 1)
      
    }
    
    # Stop searching for proxies if no SNP-outcome data available
    if( (dim(need_proxies)[1] > 0) & (dim(tmp3)[1] == 0) ){
      
      print("No proxy-outcome associations available for any SNPs.")
      
    } else {
      
      # Add column specifying proxy effect allele
      # checking proxy outcome data and using original SNP alleles
      # using query SNP alleles but checking they're the right way round using proxy outcome data
      #'Correlated alleles' format:
      # query effect allele = proxy effect allele, query other allele = proxy other allele
      tmp4 <- tmp3
      tmp4$effect_allele.proxy <- apply(tmp4, 1, function(row) {
        
        # If outcome effect allele same as proxy effect allele, treat as if using query effect allele -
        if (row["effect_allele.outcome"] == substr(row["Correlated_Alleles"], 3, 3)) {
          return(substr(row["Correlated_Alleles"], 1, 1))
          # if outcome effect allele same as proxy other allele, treat as if using query other allele -
        } else if (row["effect_allele.outcome"] == substr(row["Correlated_Alleles"], 7, 7)) {
          return(substr(row["Correlated_Alleles"], 5, 5))
        } else {
          return("default_value")
        }
      })
      # Same again but for other allele -
      tmp4$other_allele.proxy <- apply(tmp4, 1, function(row) {
        if (row["other_allele.outcome"] == substr(row["Correlated_Alleles"], 3, 3)) {
          return(substr(row["Correlated_Alleles"], 1, 1))
        } else if (row["other_allele.outcome"] == substr(row["Correlated_Alleles"], 7, 7)) {
          return(substr(row["Correlated_Alleles"], 5, 5))
        } else {
          return("default_value")
        }
      })
      
      # Treat proxy alleles as the outcome alleles
      tmp5 <- tmp4
      tmp5$effect_allele.outcome<-tmp5$effect_allele.proxy
      tmp5$other_allele.outcome<-tmp5$other_allele.proxy
      tmp5$data_source.outcome <- NA
      
      # Create final dataframe which has replaced missing SNPs with their proxies
      outcome_dat_cols <- colnames(outcome_dat)
      tmp5 <- tmp5 %>% select(all_of(outcome_dat_cols))
      
      # Return message with number of SNPs which could be proxied
      print(paste0(dim(tmp5)[1], " proxies identified."))
      
      tmp_outcome_dat <- rbind(tmp_outcome_dat, tmp5)
    }
    
  }
  
  # Return full outcome dataset with added proxies
  tmp_outcome_dat
  
}
#######

for (exp_name in names(exposure_list)) {
  
  # Convert exposure
  exposure_dat <- format_data(
    exposure_list[[exp_name]],
    type = "exposure",
    snp_col = "RSID",
    beta_col = "beta",
    se_col = "se",
    pval_col = "pval",
    eaf_col = "EAF",
    effect_allele_col = "effect_allele",
    other_allele_col = "other_allele"
  )
  exposure_dat$phenotype <- exp_name
  
  # Loop outcome
  for (outcome_name in names(outcome_list)) {
    
    outcome_data_raw <- outcome_list[[outcome_name]]
    
    outcome_data_raw <- as.data.frame(outcome_data_raw)
    
    # Convert outcome
    outcome_dat <- format_data(
      outcome_data_raw,
      type = "outcome",
      snp_col = "SNP",
      beta_col = "beta",
      se_col = "se",
      pval_col = "pval",
      eaf_col = "eaf",
      effect_allele_col = "effect_allele",
      other_allele_col = "other_allele"
    )
    
    outcome_dat$phenotype <- outcome_name
    
    outcome_dat <- generate_outdat_with_proxies(exposure_dat = exposure_dat, outcome_dat = outcome_dat, outcome_name = outcome_name,proxies = All_proxy_SNPs)
    
    # harmonise
    harmonised <- harmonise_data(exposure_dat, outcome_dat)
    harmonised$exposure <- exp_name
    harmonised$outcome <- outcome_name
    harmonised_list[[paste0(exp_name, "_", outcome_name)]] <- harmonised
    
    if (nrow(harmonised) == 0) next
    
    # MR 
    mr_result <- mr(harmonised)
    mr_result$exposure <- exp_name
    mr_result$outcome <- outcome_name
    
    all_results[[paste0(exp_name, "_", outcome_name)]] <- mr_result
  }
}

#######

for (this_exp in all_exposures) {
  
  outcomes <- unique(all_mr_summary$outcome[all_mr_summary$exposure == this_exp])
  
  plot_list <- list()
  
  for (this_out in outcomes) {
    mr_sub <- subset(all_mr_summary, exposure == this_exp & outcome == this_out)
    dat_sub <- subset(harmonised_df, exposure == this_exp & outcome == this_out)
    
    mr_sub <- subset(mr_sub, method %in% c("Inverse variance weighted", "MR Egger", "Weighted median"))
    
    if (nrow(dat_sub) >= 3 & nrow(mr_sub) > 0) {
      p <- mr_scatter_plot(mr_sub, dat_sub)[[1]] +
        ggtitle(paste0(this_exp, " → ", this_out))
      plot_list[[this_out]] <- p
    }
  }
  
  if (length(plot_list) > 0) {
    combined_plot <- wrap_plots(plot_list, ncol = 4, nrow = 3) + 
      plot_annotation(title = paste0("MR Scatter Plots for ", this_exp))
    
    file_name <- paste0("scatter_", gsub("[^A-Za-z0-9]", "_", this_exp), "_combined.png")
    
    ggsave(filename = file.path(output_dir, file_name),
           plot = combined_plot,
           width = 4 * 4,   
           height = 4 * 3,  
           dpi = 300)
  }
}

######

for (exp_name in names(exposure_list)) {
  
  # Convert exposure
  exposure_dat <- format_data(
    exposure_list[[exp_name]],
    type = "exposure",
    snp_col = "RSID",
    beta_col = "beta",
    se_col = "se",
    pval_col = "pval",
    eaf_col = "EAF",
    effect_allele_col = "effect_allele",
    other_allele_col = "other_allele"
  )
  exposure_dat$phenotype <- exp_name
  
  # Loop outcome
  for (outcome_name in names(outcome_list)) {
    
    outcome_data_raw <- outcome_list2[[outcome_name]]
    
    outcome_data_raw <- as.data.frame(outcome_data_raw)
    
    for (study_name in unique(outcome_data_raw$study)) {
      
      outcome_study_specific <- outcome_data_raw[outcome_data_raw$study == study_name,]
      
      # Convert outcome
      outcome_dat <- format_data(
        outcome_study_specific,
        type = "outcome",
        snp_col = "SNP",
        beta_col = "beta",
        se_col = "se",
        pval_col = "pval",
        eaf_col = "eaf",
        effect_allele_col = "effect_allele",
        other_allele_col = "other_allele"
      )
      
      outcome_dat$phenotype <- outcome_name
      
      outcome_dat <- generate_outdat_with_proxies(exposure_dat = exposure_dat, outcome_dat = outcome_dat, outcome_name = outcome_name,proxies = All_proxy_SNPs)
      
      # harmonise
      harmonised <- harmonise_data(exposure_dat, outcome_dat)
      harmonised$exposure <- exp_name
      harmonised$outcome <- outcome_name
      harmonised_list2[[paste0(exp_name, "_", outcome_name,"_",study_name)]] <- harmonised
      
      if (nrow(harmonised) == 0) next
      
      # MR 
      mr_result <- mr(harmonised)
      mr_result$exposure <- exp_name
      mr_result$outcome <- outcome_name
      mr_result$study <- study_name
      
      all_results2[[paste0(exp_name, "_", outcome_name,"_",study_name)]] <- mr_result
    }
  }
}


lapply(exposures2, function(exp_name) {
  
  leave1out_subset <- leave1out_results2 %>%
    filter(exposure == exp_name)
  
  overall_estimates <- leave1out_subset %>%
    filter(SNP == "All") %>%
    select(id.outcome, overall_b = b)
  
  leave1out_subset <- leave1out_subset %>%
    left_join(overall_estimates, by = "id.outcome")
  
  p <- ggplot(leave1out_subset, aes(x = b, y = SNP)) +
    geom_point() +
    geom_errorbarh(aes(xmin = b - 1.96 * se, xmax = b + 1.96 * se), height = 0.2) +
    geom_vline(aes(xintercept = overall_b), linetype = "dashed", color = "red") +
    facet_wrap(~ id.outcome, scales = "free_y") +
    labs(
      title = paste("Leave-One-Out Analysis for", exp_name),
      x = "Beta After Removing Each SNP",
      y = "SNP"
    ) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5))
  
  ggsave(
    filename = paste0(output_dir, "/", gsub("[^A-Za-z0-9]", "_", exp_name), "_leave1out_plot.png"),
    plot = p,
    width = 25,
    height = 20
  )
})

for (exp_name in names(exposure_list)) {
  
  # Convert exposure
  exposure_dat <- format_data(
    exposure_list[[exp_name]],
    type = "exposure",
    snp_col = "RSID",
    beta_col = "beta",
    se_col = "se",
    pval_col = "pval",
    eaf_col = "EAF",
    effect_allele_col = "effect_allele",
    other_allele_col = "other_allele"
  )
  exposure_dat$phenotype <- exp_name
  
  # Loop outcome
  for (outcome_name in names(outcome_list_trios)) {
    
    outcome_data_raw <- outcome_list_trios[[outcome_name]]
    
    outcome_data_raw <- as.data.frame(outcome_data_raw)
    
    # Convert outcome
    outcome_dat <- format_data(
      outcome_data_raw,
      type = "outcome",
      snp_col = "SNP",
      beta_col = "beta_mat_donuts",
      se_col = "se_mat_donuts",
      pval_col = "p_mat_donuts",
      eaf_col = "eaf_mat",
      effect_allele_col = "effect_allele",
      other_allele_col = "other_allele"
    )
    
    outcome_dat$phenotype <- outcome_name
    
    outcome_dat <- generate_outdat_with_proxies(exposure_dat = exposure_dat, outcome_dat = outcome_dat, outcome_name = outcome_name,proxies = All_proxy_SNPs)
    
    # harmonise
    harmonised <- harmonise_data(exposure_dat, outcome_dat)
    harmonised$exposure <- exp_name
    harmonised$outcome <- outcome_name
    harmonised_list_mat_donut_trios[[paste0(exp_name, "_", outcome_name)]] <- harmonised
    
    if (nrow(harmonised) == 0) next
    
    # MR 
    mr_result <- mr(harmonised)
    mr_result$exposure <- exp_name
    mr_result$outcome <- outcome_name
    
    all_results_mat_donut_trios[[paste0(exp_name, "_", outcome_name)]] <- mr_result
  }
}

n_out <- harmonised_df %>%
  group_by(outcome) %>%
  summarise(
    Ncase_min = min(ncase.outcome),
    N_min = min(samplesize.outcome),
    Ncase_median = median(ncase.outcome),
    N_median = median(samplesize.outcome),
    Ncase_max = max(ncase.outcome),
    N_max = max(samplesize.outcome)
  ) %>%
  ungroup %>%
  merge(., out_info, by.x = "outcome", by.y = "Variable_name") %>%
  arrange(priority, group, order)

n_out2 <- harmonised_df %>%
  group_by(outcome) %>%
  summarise(
    Ncase_min = min(ncase.outcome),
    Ncontrol_min = min(ncontrol.outcome),
    N_min = min(samplesize.outcome),
    Ncase_median = median(ncase.outcome),
    Ncontrol_median = median(ncontrol.outcome),
    N_median = median(samplesize.outcome),
    Ncase_max = max(ncase.outcome),
    Ncontrol_max = max(ncontrol.outcome),
    N_max = max(samplesize.outcome)
  ) %>%
  ungroup %>%
  merge(., out_info, by.x = "outcome", by.y = "Variable_name") %>%
  arrange(priority, group, order)
