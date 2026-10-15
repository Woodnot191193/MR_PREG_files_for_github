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

#Loading the required package for the task
library(ieugwasr)

#For the exposures of Any female infertility (F_all), Unknown causes of infertility (F_incl) and Anovulatory infertility (F_anov), loading the Mixed Ancestry GWAS datasets onto R, selecting the SNPs associated at genome-wide p-value significance (p < 5e-8) and running remote LD clumping to get the main loci
F_all <- read.table("GCST90483469.tsv.gz", sep ="\t",header = TRUE)
F_all_5e8 <- F_all[F_all$p_value <= 5e-8,]
rm(F_all)
names(F_all_5e8)[9] <- "rsid"
names(F_all_5e8)[8] <- "pval"
fall <- ld_clump(F_all_5e8)
rm(F_all_5e8)
F_incl <- read.table("GCST90483473.tsv.gz", sep ="\t",header =TRUE)
F_incl_5e8 <- F_incl[F_incl$p_value <= 5e-8,]
rm(F_incl)
names(F_incl_5e8)[8] <- "pval"
names(F_incl_5e8)[9] <- "rsid"
fincl <- ld_clump(F_incl_5e8)
rm(F_incl_5e8)
F_anov <- read.table(file = "GCST90483471.tsv.gz",sep = "\t",header=TRUE)
F_anov_5e8 <- F_anov[F_anov$p_value <= 5e-8,]
rm(F_anov)
names(F_anov_5e8)[8] <- "pval"
names(F_anov_5e8)[9] <- "rsid"
fanov <- ld_clump(F_anov_5e8)
rm(F_anov_5e8)
#For Anatomical infertility (F_anat) and Male infertility (male_infertility), no LD clumping was needed
F_anat <- read.table(file = "GCST90483470.tsv.gz",sep = "\t",header=TRUE)
F_anat_5e8 <- F_anat[F_anat$p_value <= 5e-8,]
rm(F_anat)
F_anat_5e8
names(F_anat_5e8)[8] <- "pval"
names(F_anat_5e8)[9] <- "rsid"
fanat <- F_anat_5e8
rm(F_anat_5e8)
male_infertility <- read.table("GCST90483474.tsv.gz",sep = "\t",header = TRUE)
male_infertility_5e8 <- male_infertility[male_infertility$p_value <= 5e-8,]
rm(male_infertility)
names(male_infertility_5e8)[8] <- "pval"
names(male_infertility_5e8)[9] <- "rsid"
mall <- male_infertility_5e8
rm(male_infertility_5e8)

#For all the exposures, loading the European ancestry GWAS, and extracting the SNP associations for the SNPs identified in the Mixed ancestry GWAS
F_all_Eur <- read.table("GCST90483463.tsv.gz",sep = "\t",header = TRUE)
fall_EUR <- F_all_Eur[F_all_Eur$rs_id %in% fall$rsid,]
names(fall_EUR)[8] <- "pval"
names(fall_EUR)[9] <- "rsid"
rm(F_all_Eur)
F_anat_Eur <- read.table("GCST90483464.tsv.gz",sep = "\t",header = TRUE)
fanat_EUR <- F_anat_Eur[F_anat_Eur$rs_id %in% fanat$rsid,]
names(fanat_EUR)[8] <- "pval"
names(fanat_EUR)[9] <- "rsid"
rm(F_anat_Eur)
F_anov_Eur <- read.table("GCST90483465.tsv.gz",sep = "\t",header = TRUE)
fanov_EUR <- F_anov_Eur[F_anov_Eur$rs_id %in% fanov$rsid,]
names(fanov_EUR)[8] <- "pval"
names(fanov_EUR)[9] <- "rsid"
rm(F_anov_Eur)
F_incl_Eur <- read.table("GCST90483467.tsv.gz",sep = "\t",header = TRUE)
fincl_EUR <- F_incl_Eur[F_incl_Eur$rs_id %in% fincl$rsid,]
names(fincl_EUR)[8] <- "pval"
names(fincl_EUR)[9] <- "rsid"
rm(F_incl_Eur)
fall_EUR$exposure <- "F-ALL"
fanat_EUR$exposure <- "F-ANAT"
fanov_EUR$exposure <- "F-ANOV"
fincl_EUR$exposure <- "F-INCL"
male_infertility_EUR <- read.table("GCST90483468.tsv.gz",sep = "\t",header = TRUE)
mall_EUR <- male_infertility_EUR[male_infertility_EUR$rs_id %in% male_infertility_5e8$rs_id,]
names(mall_EUR)[8] <- "pval"
names(mall_EUR)[9] <- "rsid"
rm(male_infertility_EUR)
mall_EUR$exposure <- "male infertility"

#Making the exposure list for the MR
exposure_list2 <- list("F-ALL" = fall_EUR,"F-INCL" = fincl_EUR,"F-ANOV" = fanov_EUR,"F-ANAT" = fanat_EUR, "M-ALL" = mall_EUR)

#Exporting the results for use in making supplementary tables
write.csv(fall_EUR,"fall_EUR.csv")
write.csv(fincl_EUR,"fincl_EUR.csv")
write.csv(fanov_EUR,"fanov_EUR.csv")
write.csv(fanat_EUR,"fanat_EUR.csv")
write.csv(mall_EUR,"mall_EUR.csv")
