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

#Loading the required packages for the task
library(LDlinkR)
library(MRutils)

#Logging into LDlink online database
Sys.setenv(LDLINK_TOKEN=#Insert relevant token here)
#Extracting all available proxies for exposure SNPs missing outcome data. Note, we were unable to find any good proxies for rs1228269928, leaving us with only 1 male infertility SNP              
All_proxy_SNPs <- get_proxies(c("rs61768001","rs17378154","rs9902027","rs75045132","rs139862664","rs1228269928"), pop="CEU",results_dir="derived_data")
#Removing duplicated SNPs from proxy SNP list
All_proxy_SNPs <- All_proxy_SNPs[!(All_proxy_SNPs$rsid %in% c("rs61768001","rs17378154","rs9902027","rs75045132","rs139862664","rs1228269928")),]
#Exporting list of proxy SNPs to be used for extracting outcome information from the server
write.csv(All_proxy_SNPs,"All_proxy_SNPs.csv")
