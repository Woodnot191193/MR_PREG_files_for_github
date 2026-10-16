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
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(ggforce)
library(ggforestplot)
library(ggpubr)
library(patchwork)
library(simex)
library(purrr)
library(stringi)

###Preparing Figures 2-6

##Preparing data frame for figures.
#Selecting correct main MR results
results <- all_mr_summary3 %>% filter(method %in% c("Inverse variance weighted","Wald ratio"))
#Changing names of outcomes
results <- results %>% mutate(outcome = case_when(
  outcome == "hdp_subsamp" ~ "HDP",
  outcome == "gh_subsamp" ~ "GH",
  outcome == "gdm_subsamp" ~ "GDM",
  outcome == "pe_subsamp" ~ "PE",
  outcome == "misc_subsamp" ~ "Miscarriage",
  outcome == "sb_subsamp" ~ "Stillbirth",
  outcome == "pretb_all" ~ "PTB",
  outcome == "vpretb_all" ~ "Very PTB",
  outcome == "sga" ~ "SGA",
  outcome == "lga" ~ "LGA",
  outcome == "lbw_all" ~ "Low Birth Weight",
  outcome == "hbw_all" ~ "High Birth Weight",
  TRUE ~ outcome
))
#Numerically ordering outcomes
results <- results %>% mutate(order1 = case_when(
outcome == "Stillbirth" ~ "1",
outcome == "Miscarriage" ~ "2",
outcome == "HDP" ~ "3",
outcome == "GH" ~ "4",
outcome == "PE" ~ "5",
outcome == "GDM" ~ "6",
outcome == "Very PTB" ~ "7",
outcome == "PTB" ~ "8",
outcome == "SGA" ~ "9",
outcome == "LGA" ~ "10",
outcome == "Low Birth Weight" ~ "11",
outcome == "High Birth Weight" ~ "12",
TRUE ~ outcome))

results$order1 <- as.numeric(results$order1)
results <- results[order(results$order1),]

#Making Odds Ratios and 95% confidence intervals for figures
results$OR <- exp(results$b)
results$OR <- round(results$OR,digits = 2)
results$LCI <- results$b - (1.96*results$se)
results$LCI <- exp(results$LCI)
results$LCI <- round(results$LCI,2)
results$UCI <- results$b + (1.96*results$se)
results$UCI <- exp(results$UCI)
results$UCI <- round(results$UCI,2)
#Making Odds Ratios and 95% confidence intervals as text variables (estimate)
results$estimate <- paste0(results$OR," (",results$LCI,",",results$UCI,")")
#Manually changing some of the estimate values to be the correct number of decimal places.
results$estimate <- ifelse(results$estimate == "1.05 (1,1.1)","1.05 (1.00,1.10)",results$estimate)
results$estimate <- ifelse(results$estimate == "1.01 (0.92,1.1)","1.01 (0.92,1.10)",results$estimate)
results$estimate <- ifelse(results$estimate == "1.06 (0.93,1.2)","1.06 (0.93,1.20)",results$estimate)
results$estimate <- ifelse(results$estimate == "0.9 (0.5,1.63)","0.90 (0.50,1.63)",results$estimate)
results$estimate <- ifelse(results$estimate == "0.93 (0.8,1.08)","0.93 (0.80,1.08)",results$estimate)
results$estimate <- ifelse(results$estimate == "1.15 (1.01,1.3)","1.15 (1.01,1.30)",results$estimate)
results$estimate <- ifelse(results$estimate == "1.46 (1.2,1.77)","1.46 (1.20,1.77)",results$estimate)
results$estimate <- ifelse(results$estimate == "1.07 (0.95,1.2)","1.07 (0.95,1.20)",results$estimate)
results$estimate <- ifelse(results$estimate == "0.6 (0.33,1.09)","0.60 (0.33,1.09)",results$estimate)
results$estimate <- ifelse(results$estimate == "1 (0.95,1.04)","1.00 (0.95,1.04)",results$estimate)
results$estimate <- ifelse(results$estimate == "1.01 (0.93,1.1)","1.01 (0.93,1.10)",results$estimate)
results$estimate <- ifelse(results$estimate == "1 (0.91,1.1)","1.00 (0.91,1.10)",results$estimate)
results$estimate <- ifelse(results$estimate == "1 (0.92,1.09)","1.00 (0.92,1.09)",results$estimate)
results$estimate <- ifelse(results$estimate == "1 (0.88,1.13)","1.00 (0.88,1.13)",results$estimate)
results$estimate <- ifelse(results$estimate == "1.02 (0.9,1.16)","1.02 (0.90,1.16)",results$estimate)

##Function to make forest plot figures
myforestplot_2 <- function(df, log_T = TRUE, xlab, limits)
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
  ) + theme(panel.grid.major.x = element_blank(),axis.text.y = element_text(size = 15)) + geom_text(aes(x = 3, label = estimate),hjust = 0) + geom_point(size = 1) + geom_effect(ggplot2::aes(xmin = .data$.xmin, xmax = .data$.xmax), 
                                                                                                                                                                                 position = ggstance::position_dodgev(height = 0.5),
                                                                                                                                                                                 size=0.5)
  colours_BP <- c("black", "gold2", "deepskyblue1", "coral", "springgreen3", "purple", "magenta", "turquoise")
  shapes_BP <- c(23L, 21L, 21L, 21L, 21L, 21L, 21L, 21L)
  x <- x + scale_color_manual(values=colours_BP)
  #x <- x + scale_x_continuous(breaks = c(0.5, 1, 1.5, 2, 2.5, 3, 5))
  x <- x + ggplot2::scale_shape_manual(values = shapes_BP)
  print(x)
}

##Making Figures 2-6. Note, the figures were exported using RStudios "Save as PDF" function.

#Figure 2: Any female infertility
F_ALL_Meta5 <- results %>% filter(exposure == "F-ALL") %>% myforestplot_2(df = ., xlab = "Odds ratio (95% CI) per 1 log-odds increased risk of \n any female infertility", limits = c(0.10,10))
#Figure 3: Anovulatory infertility
F_ANOV_Meta5 <- results %>% filter(exposure == "F-ANOV") %>% myforestplot_2(df = ., xlab = "Odds ratio (95% CI) per 1 log-odds increased risk of \n anovulatory infertility", limits = c(0.10,10))
#Figure 4: Anatomical infertility
F_ANAT_Meta5 <- results %>% filter(exposure == "F-ANAT") %>% myforestplot_2(df = ., xlab = "Odds ratio (95% CI) per 1 log-odds increased risk of \n anatomical infertility", limits = c(0.10,10))
#Figure 5: Unknown cause infertility
F_INCL_Meta5 <- results %>% filter(exposure == "F-INCL") %>% myforestplot_2(df = ., xlab = "Odds ratio (95% CI) per 1 log-odds increased risk of \n unknown cause infertility", limits = c(0.10,10))
#Figure 6: Male infertility
M_ALL_Meta5 <- results %>% filter(exposure == "M-ALL") %>% myforestplot_2(df = ., xlab = "Odds ratio (95% CI) per 1 log-odds increased risk of \n male infertility", limits = c(0.10,10))


###Preparing Figures 7-10
##Preparing data frame for figures
#The core file was prepared in excel from exported data
new_methods_results5 <- read.csv("new_methods_forest_plot_data4.csv")
#Making Odds Ratios and 95% confidence intervals for figures
new_methods_results5$OR <- exp(new_methods_results5$b)
new_methods_results5$OR <- round(new_methods_results5$OR,digits = 2)
new_methods_results5$LCI <- new_methods_results5$b - (1.96*new_methods_results5$se)
new_methods_results5$LCI <- exp(new_methods_results5$LCI)
new_methods_results5$LCI <- round(new_methods_results5$LCI,2)
new_methods_results5$UCI <- new_methods_results5$b + (1.96*new_methods_results5$se)
new_methods_results5$UCI <- exp(new_methods_results5$UCI)
new_methods_results5$UCI <- round(new_methods_results5$UCI,2)
#Making Odds Ratios and 95% confidence intervals as text variables (estimate)
new_methods_results5$estimate <- paste0(new_methods_results5$OR," (",new_methods_results5$LCI,"-",new_methods_results5$UCI,")")

#Numerically ordering methods
new_methods_results5 <- new_methods_results5 %>% mutate(order1 = case_when(
method == "Inverse variance weighted" ~ "1",
method == "MR Egger" ~ "2",
method == "Weighted median" ~ "3",
method == "Fetal and paternal genotype adjusted" ~ "4",
TRUE ~ method))

new_methods_results5$order1 <- as.numeric(new_methods_results5$order1)
new_methods_results5 <- new_methods_results5[order(new_methods_results5$order1),]

#Manually changing some of the estimate values to be the correct number of decimal places.
new_methods_results5$estimate <- ifelse(new_methods_results5$estimate == "0.5 (0.3,0.84)","0.50 (0.30,0.84)",new_methods_results5$estimate)
new_methods_results5$estimate <- ifelse(new_methods_results5$estimate == "5.66 (0.1,306.51)","5.66 (0.10,0.306.51)",new_methods_results5$estimate)
new_methods_results5$estimate <- ifelse(new_methods_results5$estimate == "0.3 (0.11,0.8)","0.30 (0.11,0.80)",new_methods_results5$estimate)
new_methods_results5$estimate <- ifelse(new_methods_results5$estimate == "5.66 (0.10,0.306.51)","5.66 (0.10,306.51)",new_methods_results5$estimate)
new_methods_results5$estimate <- ifelse(new_methods_results5$estimate == "1.23 (0,25222.47)","1.23 (0.00,25222.47)",new_methods_results5$estimate)
new_methods_results5$estimate <- ifelse(new_methods_results5$estimate == "1.15 (1.01,1.3)","1.15 (1.01,1.30)",new_methods_results5$estimate)
new_methods_results5$estimate <- ifelse(new_methods_results5$estimate == "1.07 (0.63,1.8)","1.07 (0.63,1.80)",new_methods_results5$estimate)
new_methods_results5$estimate <- ifelse(new_methods_results5$estimate == "0.98 (0.81,1.2)","0.98 (0.81,1.20)",new_methods_results5$estimate)
new_methods_results5$estimate <- ifelse(new_methods_results5$estimate == "1.05 (1,1.1)","1.05 (1.00,1.10)",new_methods_results5$estimate)
new_methods_results5$estimate <- ifelse(new_methods_results5$estimate == "1.09 (0.9,1.31)","1.09 (0.90,1.31)",new_methods_results5$estimate)


##Function to make forest plots for sensitivity analyses
myforestplot_3 <- function(df, log_T = TRUE, xlab, limits)
{
  x <- forestplot(
    df = df,
    estimate = b,
    se = se,
    pvalue = pval,
    name = method,
    logodds = log_T,
    #title = title,
    xlab = xlab,
    xlim= limits,
    colour = method
  ) + theme(
    panel.grid.major.x = element_blank(),
    legend.position = "bottom",
    legend.justification = c(0.6,0.5),
    legend.text = element_text(size = 7),
    plot.margin = margin(1, 1, 1,1)
  ) + geom_effect(ggplot2::aes(xmin = .data$.xmin, xmax = .data$.xmax,colour = method), 
                  position = ggstance::position_dodgev(height = 0.5),
                  size=0.5) + geom_segment(data = subset(df,UCI > limits[2]),aes(x = limits[2]-0.05,xend = limits[2]+2.5,y = method,yend = method,colour = method),arrow = arrow(length = unit(4,"mm")),show.legend = FALSE) + 
    geom_segment(data = subset(df,LCI < limits[1]),aes(x = limits[1]+0.05,xend = limits[1]-0.03,y = method,yend = method,colour = method),arrow = arrow(length = unit(4,"mm")),show.legend = FALSE)  + 
    ggforce::facet_col(
      facets = ~outcome,
      scales = "free_y",
      space = "free"
    )
  colours_BP <- c("gold2", "deepskyblue1", "coral", "black", "springgreen3", "purple", "magenta", "turquoise")
  shapes_BP <- c(23L, 21L, 21L, 21L, 21L, 21L, 21L, 21L)
  x <- x + scale_color_manual(values=colours_BP)
  #x <- x + scale_x_continuous(breaks = c(0.5, 1, 1.5, 2, 2.5, 3, 5))
  x <- x + ggplot2::scale_shape_manual(values = shapes_BP)
  print(x)
}

#Had to change the function for Anatomical infertility, in order for the colour coding to match
myforestplot_3_anat <- function(df, log_T = TRUE, xlab, limits)
{
  x <- forestplot(
    df = df,
    estimate = b,
    se = se,
    pvalue = pval,
    name = method,
    logodds = log_T,
    #title = title,
    xlab = xlab,
    xlim= limits,
    colour = method
  ) + theme(
    panel.grid.major.x = element_blank(),
    legend.position = "bottom",
    legend.justification = c(0.6,0.5),
    legend.text = element_text(size = 7),
    plot.margin = margin(1, 1, 1,1)
  ) + geom_effect(ggplot2::aes(xmin = .data$.xmin, xmax = .data$.xmax,colour = method), 
                  position = ggstance::position_dodgev(height = 0.5),
                  size=0.5) + geom_segment(data = subset(df,UCI > limits[2]),aes(x = limits[2]-0.05,xend = limits[2]+2.5,y = method,yend = method,colour = method),arrow = arrow(length = unit(4,"mm")),show.legend = FALSE) + 
    geom_segment(data = subset(df,LCI < limits[1]),aes(x = limits[1]+0.05,xend = limits[1]-0.03,y = method,yend = method,colour = method),arrow = arrow(length = unit(4,"mm")),show.legend = FALSE)  + 
    ggforce::facet_col(
      facets = ~outcome,
      scales = "free_y",
      space = "free"
    )
  colours_BP <- c("gold2", "black", "springgreen3", "purple", "magenta", "turquoise")
  shapes_BP <- c(23L, 21L, 21L, 21L, 21L, 21L, 21L, 21L)
  x <- x + scale_color_manual(values=colours_BP)
  #x <- x + scale_x_continuous(breaks = c(0.5, 1, 1.5, 2, 2.5, 3, 5))
  x <- x + ggplot2::scale_shape_manual(values = shapes_BP)
  print(x)
}

##Making Figures 7-10. Note, the figures were exported using RStudios "Save as PDF" function.
#Figure 7: Any female infertility
F_ALL_Methods4 <- new_methods_results5 %>% filter(exposure == "F-ALL") %>% myforestplot_3(df = ., xlab = "Odds ratio (95% CI) per 1 log-odds increased risk of \n any female infertility", limits = c(0.15,10))
p_left <- new_methods_results5 %>% filter(exposure == "F-ALL") %>% ggplot(aes(y = factor(method,levels = rev(unique(method))))) + geom_text(aes(x = 0, label = estimate),hjust = 0.55) + scale_x_continuous(expand = c(0, 0)) + theme_void() + theme( plot.margin = margin(5.5, 0, 5.5, 0),strip.text = element_text(colour = "white")) + ggforce::facet_col(facets = ~outcome, scales = "free_y", space = "free")
F_ALL_Methods4 <- F_ALL_Methods4 + p_left + plot_layout(widths = c(3,2.3))
F_ALL_Methods4
#Figure 8: Anovulatory infertility
F_ANOV_Methods4 <- new_methods_results5 %>% filter(exposure == "F-ANOV") %>% myforestplot_3(df = ., xlab = "Odds ratio (95% CI) per 1 log-odds increased risk of \n anovulatory infertility", limits = c(0.15,10))
p_left <- new_methods_results5 %>% filter(exposure == "F-ANOV") %>% ggplot(aes(y = factor(method,levels = rev(unique(method))))) + geom_text(aes(x = 0, label = estimate),hjust = 0.55) + scale_x_continuous(expand = c(0, 0)) + theme_void() + theme( plot.margin = margin(5.5, 0, 5.5, 0),strip.text = element_text(colour = "white")) + ggforce::facet_col(facets = ~outcome, scales = "free_y", space = "free")
F_ANOV_Methods4 <- F_ANOV_Methods4 + p_left + plot_layout(widths = c(3,2.3))
F_ANOV_Methods4
#Figure 9: Anatomical infertility
F_ANAT_Methods4 <- new_methods_results5 %>% filter(exposure == "F-ANAT") %>% myforestplot_3_anat(df = ., xlab = "Odds ratio (95% CI) per 1 log-odds increased risk of \n anatomical infertility", limits = c(0.15,10))
p_left <- new_methods_results5 %>% filter(exposure == "F-ANAT") %>% ggplot(aes(y = factor(method,levels = rev(unique(method))))) + geom_text(aes(x = 0, label = estimate),hjust = 0.55) + scale_x_continuous(expand = c(0, 0)) + theme_void() + theme( plot.margin = margin(5.5, 0, 5.5, 0),strip.text = element_text(colour = "white")) + ggforce::facet_col(facets = ~outcome, scales = "free_y", space = "free")
F_ANAT_Methods4 <- F_ANAT_Methods4 + p_left + plot_layout(widths = c(3,2.3))
F_ANAT_Methods4
#Figure 10: Unknown cause infertility
F_INCL_Methods4 <- new_methods_results5 %>% filter(exposure == "F-INCL") %>% myforestplot_3(df = ., xlab = "Odds ratio (95% CI) per 1 log-odds increased risk of \n unknown cause infertility", limits = c(0.15,10))
p_left <- new_methods_results5 %>% filter(exposure == "F-INCL") %>% ggplot(aes(y = factor(method,levels = rev(unique(method))))) + geom_text(aes(x = 0, label = estimate),hjust = 0.55) + scale_x_continuous(expand = c(0, 0)) + theme_void() + theme( plot.margin = margin(5.5, 0, 5.5, 0),strip.text = element_text(colour = "white")) + ggforce::facet_col(facets = ~outcome, scales = "free_y", space = "free")
F_INCL_Methods4 <- F_INCL_Methods4 + p_left + plot_layout(widths = c(3,2.3))
F_INCL_Methods4





