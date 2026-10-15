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

#Due to restrictions for server access, the SNP associations for the exposure and proxy SNPs with the MR-PREG outcomes were extracted on the server by Nancy McBride, and the results were moved to the shared drive to be downloaded.

#Loading the required package for the task
library(data.table)

#Defining the phenotypes used the the analysis, namely Hypertensive Disorders of Pregnancy (hdp_subsamp), Gestational Hypertension (gh_subsamp), Gestational Diabetes Mellitus (gdm_subsamp), Pre-Eclampsia (pe_subsamp), Stillbirth (sb_subsamp), Miscarriage (misc_subsamp),
#Very Preterm Birth (vpretb_all), Preterm Birth (pretb_all), Small for Gestational Age (SGA), Large for Gestational Age (LGA), Low Birth Weight (lbw_all) and High Birth Weight (hbw_all)
pheno <- c("hdp_subsamp", "gh_subsamp", "gdm_subsamp", "pe_subsamp",
"sb_subsamp", "misc_subsamp",
"vpretb_all", "pretb_all", "sga", "lga", "lbw_all", "hbw_all")

#Downloading main meta-analysis outcome SNP associations from the shared drive for the main analyses
ma_new <- fread("~/ma_out_dat.txt")
#Selecting the outcome phenotypes used for the analyses
ma_new_1 <- subset(ma_new, Phenotype %in% pheno)
rm(ma_new)
#Making the outcome list for the main MR
outcome_list2 <- split(ma_new_1, ma_new_1$Phenotype)

#Downloading the WLM SNP associations adjusted for fetal and paternal genotypes for sensitivity analyses (same steps as above)
trios <- fread("~/trios_out_dat.txt")
trios1 <- subset(trios, Phenotype %in% pheno)
rm(trios)
outcome_list_trios <- split(trios1, trios1$Phenotype)

#Downloading the SNP associations stratified by cohort for the post-hoc leave-one-cohort-out analyses of Preterm and Very Preterm Birth
stu <- fread("~/stu_out_dat.txt")
stu1 <- subset(stu, Phenotype %in% pheno)
rm(stu)
outcome_list_stu <- split(stu1, stu1$Phenotype)





