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

#Loading required packages for this task
library(TwoSampleMR)
library(dplyr)
library(tidyr)

##Estimating pseudo-F statististics
#Preparing data-frame for analysis. Note, a SNP for Any female infertility was removed, as it was palindromic and had a MAF > 0.42
exposure_all <- rbind(fall_EUR,fanat_EUR,fanov_EUR,fexcl_EUR,fincl_EUR,mall_EUR)
exposure_all <- exposure_all[exposure_all$rsid != "rs10929759",]
exposure_all <- exposure_all %>% mutate(exposure = case_when(
exposure == "F-ALL" ~ "Any female infertility",
exposure == "F-ANAT" ~ "Anatomical infertility",
exposure == "F-ANOV" ~ "Anovulatory infertility",
exposure == "F-INCL" ~ "Unknown causes of female infertility",
exposure == "male infertility" ~ "Male infertility",
TRUE ~ exposure)
)

#Estimating pseudo-F statistics for all exposure SNPs
exposure_all$F_stat <-(exposure_all$beta^2)/(exposure_all$standard_error^2)

#Estimating mean pseudo-F statistics for all exposure IVs
mean_f1 <- mean(exposure_all$F_stat[exposure_all$exposure == "Any female infertility"], na.rm = TRUE)
mean_f1
mean_f2 <- mean(exposure_all$F_stat[exposure_all$exposure == "Anatomical infertility"], na.rm = TRUE)
mean_f2
mean_f3 <- mean(exposure_all$F_stat[exposure_all$exposure == "Anovulatory infertility"], na.rm = TRUE)
mean_f3
mean_f4 <- mean(exposure_all$F_stat[exposure_all$exposure == "Unknown causes of female infertility"], na.rm = TRUE)
mean_f4
mean_f5 <- mean(exposure_all$F_stat[exposure_all$exposure == "Male infertility"], na.rm = TRUE)
mean_f5

##Getting the Cochran's Q heterogeneity results
harmonised_df4 <- harmonised_df3 %>% mutate(id.outcome = outcome)
heterogeneity_results <- mr_heterogeneity(harmonised_df4)
heterogeneity_results
##Getting the MR-Egger intercept results
pleiotropy_test <- mr_pleiotropy_test(dat = harmonised_df4)
pleiotropy_test

##Doing MR with SNP outcome associations adjusted for fetal and paternal genotypes

#Creating list for harmonised SNP associations
harmonised_list5 <- list()
#Creating list for MR results
all_results_mat_donut_trios2 <- list()

#Main code to perform MR for all exposures on all of the selected outcomes
for (exp_name in names(exposure_list2)) {
#Convert exposure
exposure_dat <- format_data(
exposure_list2[[exp_name]],
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
#Loop outcome
for (outcome_name in names(outcome_list_trios)) {
outcome_data_raw <- outcome_list_trios[[outcome_name]]
outcome_data_raw <- as.data.frame(outcome_data_raw)
#Convert outcome
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
outcome_dat <- generate_outdat_with_proxies2(exposure_dat = exposure_dat, outcome_dat = outcome_dat, outcome_name = outcome_name,proxies = All_proxy_SNPs)
#Harmonise
harmonised <- harmonise_data(exposure_dat, outcome_dat)
harmonised$exposure <- exp_name
harmonised$outcome <- outcome_name
harmonised_list5[[paste0(exp_name, "_", outcome_name)]] <- harmonised
if (nrow(harmonised) == 0) next
#MR
mr_result <- mr(harmonised)
mr_result$exposure <- exp_name
mr_result$outcome <- outcome_name
all_results_mat_donut_trios2[[paste0(exp_name, "_", outcome_name)]] <- mr_result
}
}

#Converting list of MR results into data frame
mr_summary_donut_trios2 <- dplyr::bind_rows(all_results_mat_donut_trios2)
#Exporting MR results for use in preparing file for figures
write.csv(mr_summary_donut_trios2, "mr_summary_results_donut_trios_EUR_2025_proxies_mall_updated.csv")












