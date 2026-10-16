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
library(stats)

#Creating list for harmonised SNP associations
harmonised_list3 <- list()
#Creating list for MR results
all_results3 <- list()
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
for (outcome_name in names(outcome_list2)) {
outcome_data_raw <- outcome_list2[[outcome_name]]
outcome_data_raw <- as.data.frame(outcome_data_raw)
#Convert outcome
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
#Replacing the SNP outcome associations for SNPs with missing outcome data with the SNP outcome associations of proxy SNPs
outcome_dat <- generate_outdat_with_proxies2(exposure_dat = exposure_dat, outcome_dat = outcome_dat, outcome_name = outcome_name,proxies = All_proxy_SNPs)
#Harmonise
harmonised <- harmonise_data(exposure_dat, outcome_dat)
harmonised$exposure <- exp_name
harmonised$outcome <- outcome_name
harmonised_list3[[paste0(exp_name, "_", outcome_name)]] <- harmonised
if (nrow(harmonised) == 0) next
#MR
mr_result <- mr(harmonised)
mr_result$exposure <- exp_name
mr_result$outcome <- outcome_name
all_results3[[paste0(exp_name, "_", outcome_name)]] <- mr_result
}
}

#Converting list of MR results into data frame
all_mr_summary3 <- dplyr::bind_rows(all_results3)

#Converting list of harmonised SNPs into data frame
harmonised_df3 <- bind_rows(harmonised_list3, .id = "harmonised_id")

#Exporting harmonised SNPs for Supplementary Table
write.csv(harmonised_df3, "harmonised_results_EUR_2025_proxies_mall_updated.csv", row.names = FALSE)

#Exporting all MR results for Supplementary Table (Odds ratios and 95% confidence intervals estimated separately)
write.csv(all_mr_summary3, "all_mr_summary_results_EUR_2025_proxies_mall_updated.csv", row.names = FALSE)

#Doing FDR corrected p-values for MR results
all_IVW_for_FDR3 <- all_mr_summary3 %>%
filter(method %in% c("Inverse variance weighted", "Wald ratio")) %>% filter(!exposure %in% c("M-ALL"))

all_IVW_for_FDR3$pval_fdr <- p.adjust(all_IVW_for_FDR3$pval,method = "BH")

#Exporting MR results with FDR corrected p-values for Supplementary Table (Odds ratios and 95% confidence intervals estimated separately)
write.csv(all_IVW_for_FDR3,"all_IVW_for_FRD3.csv")





