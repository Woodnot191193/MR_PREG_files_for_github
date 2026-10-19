# Load environment variables from config file
readRenviron("config.env")

# Retrieve project directory from environment
projectdir <- Sys.getenv("MRPregdir")
print(projectdir)

# Use `projectdir` to navigate folders without hardcoding paths
setwd(paste0(projectdir, "/MR_PREG_files_for_github/Main_code_and_analysis_plan_for_infertility_paper/"))
getwd()

# Return to main project directory
setwd(projectdir)
getwd()

###Loading packages required for task
library(TwoSampleMR)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(ggforce)
library(ggforestplot)
library(ggpubr)
library(simex)
library(purrr)
library(stringi)

###Forest plot function for leave-one-out and leave-one-cohort-out analysis
myforestplot_5 <- function(df, log_T = TRUE, xlab, limits)
{
  x <- forestplot(
    df = df,
    estimate = b,     
    se = se,
    pvalue = p,
    name = SNP,
    logodds = log_T,
    #title = title,
    xlab = xlab,
    xlim= limits
  ) + theme(panel.grid.major.x = element_blank(),axis.text.y = element_text(size = 15)) + geom_point(size = 1) + geom_effect(ggplot2::aes(xmin = .data$.xmin, xmax = .data$.xmax), 
                                                                                                                             position = ggstance::position_dodgev(height = 0.5),
                                                                                                                             size=2)
  colours_BP <- c("black", "gold2", "deepskyblue1", "coral", "springgreen3", "purple", "magenta", "turquoise")
  shapes_BP <- c(23L, 21L, 21L, 21L, 21L, 21L, 21L, 21L)
  x <- x + scale_color_manual(values=colours_BP)
  #x <- x + scale_x_continuous(breaks = c(0.5, 1, 1.5, 2, 2.5, 3, 5))
  x <- x + ggplot2::scale_shape_manual(values = shapes_BP)
  print(x)
}

###Doing single SNP, leave-one-out and leave-one-cohort-out for PTB
##Selecting the required exposure and outcome for the analysis
dat_sub <- subset(harmonised_df3, exposure == "F-ALL" & outcome == "pretb_all")
##Renaming the exposure and outcome
dat_sub$exposure <- "Any female infertility"
dat_sub$outcome <- "PTB"
##Doing single SNP formatting
dat_sub2 <- mr_singlesnp(dat_sub)
##Doing leave-one-out analysis
p3 <- mr_leaveoneout(dat_sub)
##Preparing single SNP data frame for forest plot
#Adding 95 % confidence intervals
dat_sub2$LCI <- dat_sub2$b - (1.96*dat_sub2$se)
dat_sub2$UCI <- dat_sub2$b + (1.96*dat_sub2$se)
#Removing MR-Egger estimate
dat_sub3 <- dat_sub2[dat_sub2$SNP != "All - MR Egger",]
#Ordering the results for forest plot
dat_sub3 <- dat_sub3[order(dat_sub3$b,decreasing = TRUE),]
dat_sub3$order <- ifelse(dat_sub3$SNP == "All - Inverse variance weighted",2,1)
dat_sub3 <- dat_sub3[order(dat_sub3$order),]
##Preparing leave-one-out results for forest plot
#Ordering the results for forest plot
p3 <- p3[order(p3$b,decreasing = TRUE),]
p3$order <- ifelse(p3$SNP == "All",2,1)
p3 <- p3[order(p3$order),]

##Making figures for single SNP and leave-one-out analysis. Note, the figures were exported as PDFs using Rstudios in built function.
#Single SNP figure
p2_own_function <- dat_sub3 %>% myforestplot_5(df = ., xlab = "MR effect size for \n any female infertility on PTB", limits = c(0.15,10))
#Leave-one-out figure
p4_own_function <- p3 %>% myforestplot_5(df = ., xlab = "MR leave-one-out analysis for \n any female infertility on PTB", limits = c(0.15,10))

##Function developed by Vibert et al 2026 for Endomtriosis MR study for leave-one-cohort out analysis
run_loo_for_outcome <- function(one_outcome_df) {
  outcome_id   <- unique(one_outcome_df$id.outcome)[1]
  outcome_name <- unique(one_outcome_df$outcome)[1]
  cohorts      <- sort(unique(one_outcome_df$study))
  
  message(">> Outcome: ", outcome_id, " | ", outcome_name,
          " | Cohorts: ", length(cohorts))
  
  # Main meta-analysis across all cohorts
  main_meta <- meta_outcome_after_drop(one_outcome_df)
  main_outcome_fmt <- build_outcome_formatted(main_meta, outcome_id, outcome_name)
  if (is.null(main_outcome_fmt) || nrow(main_outcome_fmt) < 3) return(NULL)
  
  main_harm <- tryCatch(
    TwoSampleMR::harmonise_data(exposure_formatted, main_outcome_fmt),
    error = function(e) NULL
  )
  if (is.null(main_harm) || dplyr::n_distinct(main_harm$SNP) < 3) return(NULL)
  
  main_ivw <- TwoSampleMR::mr(main_harm, method_list = "mr_ivw")
  if (!nrow(main_ivw)) return(NULL)
  
  main_ivw$analysis_type       <- "Main analysis"
  main_ivw$nsnp                <- dplyr::n_distinct(main_harm$SNP)
  main_ivw$left_out_cohort     <- "None"
  main_ivw$n_snps_remaining    <- main_ivw$nsnp
  main_ivw$n_cohorts_remaining <- length(cohorts)
  main_ivw$outcome             <- outcome_name
  main_ivw$id.outcome          <- outcome_id
  
  # Leave-one-cohort-out analyses
  loo_list <- lapply(cohorts, function(left_out) {
    sub_df <- dplyr::filter(one_outcome_df, study != left_out)
    
    meta   <- meta_outcome_after_drop(sub_df)
    out_fmt <- build_outcome_formatted(meta, outcome_id, outcome_name)
    if (is.null(out_fmt) || nrow(out_fmt) < 3) return(NULL)
    
    harm <- tryCatch(
      TwoSampleMR::harmonise_data(exposure_formatted, out_fmt),
      error = function(e) NULL
    )
    if (is.null(harm) || dplyr::n_distinct(harm$SNP) < 3) return(NULL)
    
    res <- TwoSampleMR::mr(harm, method_list = "mr_ivw")
    if (!nrow(res)) return(NULL)
    
    res$analysis_type       <- paste0("Without ", left_out)
    res$left_out_cohort     <- left_out
    res$n_snps_remaining    <- dplyr::n_distinct(harm$SNP)
    res$n_cohorts_remaining <- dplyr::n_distinct(sub_df$study)
    res$outcome             <- outcome_name
    res$id.outcome          <- outcome_id
    
    res
  })
  
  dplyr::bind_rows(main_ivw, dplyr::bind_rows(loo_list))
}

##Preparing exposure data
exposure_formatted <- format_data(
  exposure_list2[["F-ALL"]],
  type = "exposure",
  snp_col = "rsid",
  beta_col = "beta",
  se_col = "standard_error",
  pval_col = "pval",
  eaf_col = "effect_allele_frequency",
  effect_allele_col = "effect_allele",
  other_allele_col = "other_allele"
)
exposure_formatted$phenotype <- "F-ALL"
exposure_formatted <- exposure_formatted[exposure_formatted$SNP != "rs10929759",]

##Preparing outcome data
#Selecting PTB data and Any female infertility SNPs
PTB_study_specific3 <- outcomes_list_stu[["pretb_all"]]
PTB_study_specific3 <- PTB_study_specific3[PTB_study_specific3$SNP %in% c("rs61768001","rs17803970","rs75045132","rs9643050"),]
#Manually harmonizing SNPs, as the SNPs were not harmonized the same across cohorts
PTB_study_specific3$True_EA <- NA
PTB_study_specific3$True_EA <- ifelse(PTB_study_specific3$SNP == "rs61768001","T",PTB_study_specific3$True_EA)
PTB_study_specific3$True_EA <- ifelse(PTB_study_specific3$SNP == "rs17803970","A",PTB_study_specific3$True_EA)
PTB_study_specific3$True_EA <- ifelse(PTB_study_specific3$SNP == "rs75045132","C",PTB_study_specific3$True_EA)
PTB_study_specific3$True_EA <- ifelse(PTB_study_specific3$SNP == "rs9643050","T",PTB_study_specific3$True_EA)
PTB_study_specific3$eaf.outcome <- ifelse(PTB_study_specific3$True_EA != PTB_study_specific3$effect_allele.outcome,1 - PTB_study_specific3$eaf.outcome,PTB_study_specific3$eaf.outcome)
PTB_study_specific3$beta.outcome <- ifelse(PTB_study_specific3$True_EA != PTB_study_specific3$effect_allele.outcome,PTB_study_specific3$beta.outcome*-1,PTB_study_specific3$beta.outcome)
PTB_study_specific3$other_allele.outcome <- ifelse(PTB_study_specific3$True_EA != PTB_study_specific3$effect_allele.outcome,PTB_study_specific3$effect_allele.outcome,PTB_study_specific3$other_allele.outcome)
PTB_study_specific3$effect_allele.outcome <- PTB_study_specific3$True_EA
#Removing duplicate cohorts or cohorts excluded from the main meta-analysis
PTB_study_specific4 <- PTB_study_specific3[!duplicated(PTB_study_specific3[,c("SNP","study")]),]
PTB_study_specific5 <- PTB_study_specific4[PTB_study_specific4$study != "FinnGen",]

##Running the leave-one-cohort out analysis
results_list_stu_PTB_corrected3 <- run_loo_for_outcome(PTB_study_specific5)
#Converting leave-one-cohort out results to data frame
PTB_corrected_results_raw_stu3 <- dplyr::bind_rows(results_list_stu_PTB_corrected3)

#Preparing file forest plot function
PTB_corrected_results_raw_stu3 <- rename(PTB_corrected_results_raw_stu3, "SNP" = "analysis_type")
PTB_corrected_results_raw_stu3 <- rename(PTB_corrected_results_raw_stu3, "p" = "pval")
PTB_corrected_results_raw_stu3$SNP <- ifelse(PTB_corrected_results_raw_stu3$SNP == "Without Public","Without Solé-Navais et al",PTB_corrected_results_raw_stu3$SNP)

#Making figure of leave-one-cohort-out analysis. The results were exported as a PDF using RStudios in built function.
PTB_corrected_results_raw_stu_plots <- PTB_corrected_results_raw_stu3 %>% myforestplot_5(df = ., xlab = "MR effect size for \n any female infertility on PTB", limits = c(0.10,10))

###Doing single SNP, leave-one-out and leave-one-cohort-out for very PTB
##Selecting the required exposure and outcome for the analysis
dat_sub_vpret <- subset(harmonised_df3, exposure == "F-ALL" & outcome == "vpretb_all")
##Renaming the exposure and outcome
dat_sub_vpret$exposure <- "Any female infertility"
dat_sub_vpret$outcome <- "Very PTB"
##Doing single SNP formatting
dat_sub_vpret2 <- mr_singlesnp(dat_sub_vpret)
##Doing leave-one-out analysis
p3_vpret <- mr_leaveoneout(dat_sub_vpret)
##Preparing single SNP data frame for forest plot
#Adding 95 % confidence intervals
dat_sub_vpret2$LCI <- dat_sub_vpret2$b - (1.96*dat_sub_vpret2$se)
dat_sub_vpret2$UCI <- dat_sub_vpret2$b + (1.96*dat_sub_vpret2$se)
#Removing MR-Egger estimate
dat_sub_vpret3 <- dat_sub_vpret2[dat_sub_vpret2$SNP != "All - MR Egger",]
#Ordering the results for forest plot
dat_sub_vpret3 <- dat_sub_vpret3[order(dat_sub_vpret3$b,decreasing = TRUE),]
dat_sub_vpret3$order <- ifelse(dat_sub_vpret3$SNP == "All - Inverse variance weighted",2,1)
dat_sub_vpret3 <- dat_sub3[order(dat_sub_vpret3$order),]
##Making figure for single SNP results
p2_own_function_vpret <- dat_sub_vpret3 %>% myforestplot_5(df = ., xlab = "MR effect size for any female infertility on Very PTB", limits = c(0.0125,10))
##Preparing leave-one-out data frame for forest plot
p3_vpret <- p3_vpret[order(p3_vpret$b,decreasing = TRUE),]
p3_vpret$order <- ifelse(p3_vpret$SNP == "All",2,1)
p3_vpret <- p3_vpret[order(p3_vpret$order),]
##Making figure for leave-one-out analysis
p4_own_function_vpret <- p3_vpret %>% myforestplot_5(df = ., xlab = "MR leave-one-out analysis for \n any female infertility on Very PTB", limits = c(0.05,10))

##Preparing outcome data for leave-one-cohort out analysis
#Selecting PTB data and Any female infertility SNPs
VPTB_study_specific3 <- outcomes_list_stu[["vpretb_all"]]
VPTB_study_specific3 <- VPTB_study_specific3[VPTB_study_specific3$SNP %in% c("rs61768001","rs17803970","rs75045132","rs9643050"),]
#Manually harmonizing SNPs, as the SNPs were not harmonized the same across cohorts
VPTB_study_specific3$True_EA <- NA
VPTB_study_specific3$True_EA <- ifelse(VPTB_study_specific3$SNP == "rs61768001","T",VPTB_study_specific3$True_EA)
VPTB_study_specific3$True_EA <- ifelse(VPTB_study_specific3$SNP == "rs17803970","A",VPTB_study_specific3$True_EA)
VPTB_study_specific3$True_EA <- ifelse(VPTB_study_specific3$SNP == "rs75045132","C",VPTB_study_specific3$True_EA)
VPTB_study_specific3$True_EA <- ifelse(VPTB_study_specific3$SNP == "rs9643050","T",VPTB_study_specific3$True_EA)
VPTB_study_specific3$eaf.outcome <- ifelse(VPTB_study_specific3$True_EA != VPTB_study_specific3$effect_allele.outcome,1 - VPTB_study_specific3$eaf.outcome,VPTB_study_specific3$eaf.outcome)
VPTB_study_specific3$beta.outcome <- ifelse(VPTB_study_specific3$True_EA != VPTB_study_specific3$effect_allele.outcome,VPTB_study_specific3$beta.outcome*-1,VPTB_study_specific3$beta.outcome)
VPTB_study_specific3$other_allele.outcome <- ifelse(VPTB_study_specific3$True_EA != VPTB_study_specific3$effect_allele.outcome,VPTB_study_specific3$effect_allele.outcome,VPTB_study_specific3$other_allele.outcome)
VPTB_study_specific3$effect_allele.outcome <- VPTB_study_specific3$True_EA
VPTB_study_specific3 <- VPTB_study_specific3[!duplicated(VPTB_study_specific3[,c("SNP","study")]),]
#Removing duplicate cohorts or cohorts excluded from the main meta-analysis
VPTB_study_specific3 <- VPTB_study_specific3[VPTB_study_specific3$study != "FinnGen",]
VPTB_study_specific3 <- VPTB_study_specific3[VPTB_study_specific3$study != "BIB-WE",]
VPTB_study_specific3 <- VPTB_study_specific3[VPTB_study_specific3$study != "BIB-SA",]
##Running the leave-one-cohort out analysis
results_list_stu_VPTB_corrected <- run_loo_for_outcome(VPTB_study_specific3)
#Converting leave-one-cohort out results to data frame
VPTB_corrected_results_raw_stu <- dplyr::bind_rows(results_list_stu_VPTB_corrected)
#Preparing file forest plot function
VPTB_corrected_results_raw_stu <- rename(VPTB_corrected_results_raw_stu, "SNP" = "analysis_type")
VPTB_corrected_results_raw_stu <- rename(VPTB_corrected_results_raw_stu, "p" = "pval")
#Making figure of leave-one-cohort-out analysis
VPTB_corrected_results_raw_stu_plots <- VPTB_corrected_results_raw_stu %>% myforestplot_5(df = ., xlab = "MR effect size for \n any female infertility on Very PTB", limits = c(0.05,10))