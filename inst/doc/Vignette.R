## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  warning = FALSE)

## ----setup, include = FALSE---------------------------------------------------
library(MAGMA.R)
MAGMA_sim_data_tar <- readRDS(file = "data_tar_matched.rds")
MAGMA_sim_data_tar_exact <- readRDS(file ="data_tar_exact_matched.rds")
MAGMA_sim_data_2x2 <- readRDS(file = "data_2_x_2_matched.rds")
MAGMA_sim_data_2x2_exact <- readRDS(file = "data_2_x_2_exact_matched.rds")

Balance_tar <- readRDS(file = "Balance_tar.rds")
Balance_tar_exact <- readRDS(file = "Balance_tar_exact.rds")
Balance_2x2 <- readRDS(file = "Balance_2_x_2.rds")
Balance_2x2_exact <- readRDS(file = "Balance_2_x_2_exact.rds")

## ----data_introduction--------------------------------------------------------
str(MAGMA_sim_data)

## ----covariates_gifted--------------------------------------------------------
covariates_gifted <- c("GPA_school",
                       "IQ_score",
                       "Motivation",
                       "parents_academic",
                       "gender")

## ----unbalance_gifted---------------------------------------------------------
# Estimate overall and group specific descriptive statistics and Cohen’s d 
descs_gifted_pre <- MAGMA_desc(Data = MAGMA_sim_data,
                               group = "gifted_support",
                               covariates = covariates_gifted,
                               filename = "stats_gifted_pre.docx")
descs_gifted_pre %>%
  purrr::set_names(c("Overall N", "Overall Mean", "Overall SD",
                     "No Support N", "No Support Mean", "No Support SD",
                     "Support N", "Support Mean", "Support SD",
                     "d"))

# Estimating the four balance criteria
unbalance_gifted <- initial_unbalance(Data = MAGMA_sim_data,
                                      group = "gifted_support",
                                      covariates = covariates_gifted)
unbalance_gifted

## ----common_support_gifted, fig.height = 5, fig.width = 7.5-------------------
# Density overlap in propensity scores for gifted before matching
Density_overlap(Data = MAGMA_sim_data,
                variable = "ps_gifted", 
                group = "gifted_support",
                variable_name = "Propensity Score",
                group_labels = c("No Support", "Support"),
                group_name = "Gifted Support")

## ----standard_2_group_matching------------------------------------------------
# Conducting matching for gifted support
MAGMA_sim_data_gifted <- MAGMA(Data = MAGMA_sim_data,
                               group = "gifted_support",
                               dist = "ps_gifted",
                               cores = 2)
str(MAGMA_sim_data_gifted)

## ----Balance_standard_2_group_matching----------------------------------------
# Estimating the four balance criteria iteratively over possible sample sizes
Balance_gifted <- Balance_MAGMA(Data = MAGMA_sim_data_gifted,
                                group = "gifted_support",
                                covariates = covariates_gifted,
                                step = "step") 

# Extracting balance criteria for 100 cases per group
Balance_100_criteria <- Balance_extract(Balance = Balance_gifted,
                                        samplesize = 100,
                                        effects = FALSE)
Balance_100_criteria

# Extracting pairwise effects for 100 cases per group
Balance_100_effects <- Balance_extract(Balance = Balance_gifted,
                                       samplesize = 100,
                                       effects = TRUE)
Balance_100_effects

## ----Plots_standard_2_group_matching, fig.height = 5, fig.width = 7.5---------
# Plotting balance trend over sample size
Plot_MAGMA(Balance = Balance_gifted,
           criterion = c("Pillai", "d_ratio", "mean_g", "Adj_d_ratio"))

## ----Table_standard_2_group_matching------------------------------------------
Table_MAGMA(Balance = Balance_gifted,
            filename = "Balance_gifted.docx")

## ----post_matching_standard_2-------------------------------------------------
# Computing descriptive statistics and pairwise effects for 100 cases per group
descs_gifted_post <- MAGMA_desc(Data = MAGMA_sim_data_gifted,
                                group = "gifted_support",
                                covariates = covariates_gifted,
                                covariates_ordinal = "teacher_ability_rating",
                               covariates_nominal = "enrichment",
                                step_num = 100,
                                step_var = "step",
                                filename = "stats_gifted_post.docx")

# Displaying the table with defined colum names
descs_gifted_post %>%
  purrr::set_names(c("Overall N", "Overall Mean", "Overall SD",
                     "No Support N", "No Support Mean", "No Support SD",
                     "Support N", "Support Mean", "Support SD",
                     "d"))

## ----exact_2_group_matching---------------------------------------------------
MAGMA_sim_data_gifted_exact <- MAGMA_exact(Data = MAGMA_sim_data,
                                           group = "gifted_support",
                                           dist = "ps_gifted",
                                           exact = "enrichment",
                                           cores = 2)
str(MAGMA_sim_data_gifted_exact)

## ----Balance_exact_2_group_matching, fig.height = 5, fig.width = 7.5----------
Balance_gifted_exact <- Balance_MAGMA(Data = MAGMA_sim_data_gifted_exact,
                                      group = "gifted_support",
                                      covariates = covariates_gifted,
                                      step = "step") 

# Extracting balance criteria for 100 cases per group
Balance_100_criteria_exact <- Balance_extract(Balance = Balance_gifted_exact,
                                        samplesize = 100,
                                        effects = FALSE)
Balance_100_criteria_exact

# Extracting pairwise effects for 100 cases per group
Balance_100_effects_exact <- Balance_extract(Balance = Balance_gifted_exact,
                                       samplesize = 100,
                                       effects = TRUE)
Balance_100_effects_exact

# Plotting trend over increasing sample size
Plot_MAGMA(Balance = Balance_gifted_exact,
           criterion = c("Pillai", "d_ratio", "mean_g", "Adj_d_ratio")) 

# Creating table
Table_MAGMA(Balance = Balance_gifted_exact,
            filename = "Balance_gifted_exact.docx")

# Computing descriptive statistics and pairwise effects for 100 cases per group
descs_gifted_post_exact <- MAGMA_desc(Data = MAGMA_sim_data_gifted_exact,
                                      group = "gifted_support",
                                      covariates = covariates_gifted,
                                      step_num = 100,
                                      step_var = "step",
                                      filename = "stats_gifted_post_exact.docx")

# Displaying the table with defined colum names
descs_gifted_post_exact %>%
  purrr::set_names(c("Overall N", "Overall Mean", "Overall SD",
                     "No Support N", "No Support Mean", "No Support SD",
                     "Support N", "Support Mean", "Support SD",
                     "d"))

## ----covariates_tar-----------------------------------------------------------
covariates_tar <- c("GPA_school",
                    "IQ_score",
                    "Motivation",
                    "parents_academic",
                    "gifted_support")

## ----unbalance_tar, fig.height = 5, fig.width = 7.5---------------------------
# Computing descriptive statistics and all pairwise effects for three groups
descs_tar_pre <- MAGMA_desc(Data = MAGMA_sim_data,
                            group = "teacher_ability_rating",
                            covariates = covariates_tar,
                            filename = "stats_tar_pre.docx")

descs_tar_pre %>%
  purrr::set_names(c("Overall N", "Overall Mean", "Overall SD",
                     "BA N", "BA Support Mean", "BA Support SD",
                     "A N", "A Mean", "A SD",
                     "AA N", "AA Mean", "AA SD",
                     "d BA-A", "d BA-AA", "d A-AA"))

# Estimating and printing initial unbalance for teacher rated ability
unbalance_tar <- initial_unbalance(Data = MAGMA_sim_data,
                                   group = "teacher_ability_rating",
                                   covariates = covariates_tar)
unbalance_tar

# Estimating and plotting density overlap in teacher rated ability propensity scores
# Returns vector for each pairwise overlap
Density_overlap(Data = MAGMA_sim_data,
                variable = "ps_tar", 
                group = "teacher_ability_rating",
                variable_name = "Propensity Score",
                group_labels = c("Below Average", "Average", "Above Average"),
                group_name = "Teacher Rated Ability")

## ----standard_3_group_matching, eval = FALSE----------------------------------
#  MAGMA_sim_data_tar <- MAGMA(Data = MAGMA_sim_data,
#                              group = "teacher_ability_rating",
#                              dist = "ps_tar",
#                              cores = 2)

## ----standard_3_group_matching_str--------------------------------------------
str(MAGMA_sim_data_tar)

## ----Balance_3_group_matching, eval = FALSE-----------------------------------
#  Balance_tar <- Balance_MAGMA(Data = MAGMA_sim_data_tar,
#                               group = "teacher_ability_rating",
#                               covariates = covariates_tar,
#                               step = "step")

## ----Balance_3_group_matching_results, fig.height = 5, fig.width = 7.5--------
# Balance criteria for 100 cases per group
Balance_100_tar_criteria <- Balance_extract(Balance = Balance_tar,
                                            samplesize = 100,
                                            effects = FALSE)
Balance_100_tar_criteria

# Extracting pairwise effects for 100 cases per group
Balance_100_tar_effects <- Balance_extract(Balance = Balance_tar,
                                           samplesize = 100,
                                           effects = TRUE)
Balance_100_tar_effects

# Plotting trend over increasing sample size
Plot_MAGMA(Balance = Balance_tar,
           criterion = c("Pillai", "d_ratio", "mean_g", "Adj_d_ratio")) 

# Creating table
Table_MAGMA(Balance = Balance_tar,
            filename = "Balance_tar.docx")


## ----post_matching_standard_3-------------------------------------------------
# Computing descriptive statistics and pairwise effects for 100 cases per group
descs_tar_post <- MAGMA_desc(Data = MAGMA_sim_data_tar,
                            group = "teacher_ability_rating",
                            covariates = covariates_tar,
                            step_num = 100,
                            step_var = "step",
                            filename = "stats_tar_post.docx")

# Displaying the table with defined colum names
descs_tar_post %>%
    purrr::set_names(c("Overall N", "Overall Mean", "Overall SD",
                       "BA N", "BA Support Mean", "BA Support SD",
                       "A N", "A Mean", "A SD",
                       "AA N", "AA Mean", "AA SD",
                       "d BA-A", "d BA-AA", "d A-AA"))

## ----exact_3_group_matching, eval = FALSE-------------------------------------
#  MAGMA_sim_data_tar_exact <- MAGMA_exact(Data = MAGMA_sim_data,
#                                          group = "teacher_ability_rating",
#                                          dist = "ps_tar",
#                                          exact = "gender",
#                                          cores = 2)

## ----exact_3_group_matching_str-----------------------------------------------
str(MAGMA_sim_data_tar_exact)

## ----exact_3_group_matching_balance, eval = FALSE-----------------------------
#  Balance_tar_exact <- Balance_MAGMA(Data = MAGMA_sim_data_tar_exact,
#                                     group = "teacher_ability_rating",
#                                     covariates = covariates_tar,
#                                     step = "step")

## ----exact_3_group_matching_results, fig.height = 5, fig.width = 7.5----------
# Balance criteria for 100 cases per group
Balance_100_tar_criteria_exact <- Balance_extract(Balance = Balance_tar_exact,
                                                  samplesize = 100,
                                                  effects = FALSE)
Balance_100_tar_criteria_exact

# Extracting pairwise effects for 100 cases per group
Balance_100_tar_effects_exact <- Balance_extract(Balance = Balance_tar_exact,
                                                 samplesize = 100,
                                                 effects = TRUE)
Balance_100_tar_effects_exact

# Plotting trend over increasing sample size
Plot_MAGMA(Balance = Balance_tar_exact,
           criterion = c("Pillai", "d_ratio", "mean_g", "Adj_d_ratio")) #Could be omitted

# Creating table
Table_MAGMA(Balance = Balance_tar_exact,
            filename = "Balance_tar_exact.docx")

# Computing descriptive statistics and pairwise effects for 100 cases per group
descs_tar_post_exact <- MAGMA_desc(Data = MAGMA_sim_data_tar_exact,
                                   group = "teacher_ability_rating",
                                   covariates = covariates_tar,
                                   step_num = 100,
                                   step_var = "step",
                                   filename = "stats_tar_post_exact.docx")

# Displaying the table with defined column names
descs_tar_post_exact %>%
    purrr::set_names(c("Overall N", "Overall Mean", "Overall SD",
                       "BA N", "BA Support Mean", "BA Support SD",
                       "A N", "A Mean", "A SD",
                       "AA N", "AA Mean", "AA SD",
                       "d BA-A", "d BA-AA", "d A-AA"))

## ----covariates_2x2, fig.height = 5, fig.width = 7.5--------------------------
# Defining the covariates
covariates_2x2 <- c("GPA_school",
                    "IQ_score",
                    "Motivation",
                    "parents_academic",
                    "gender")

# Computing descriptive statistics and all pairwise effects
descs_2x2_pre <- MAGMA_desc(Data = MAGMA_sim_data,
                            group = c("gifted_support", "enrichment"),
                            covariates = covariates_2x2,
                            filename = "stats_2x2_pre.docx")

descs_2x2_pre %>%
  purrr::set_names(c("Overall N", "Overall Mean", "Overall SD",
                     "Sup & No En N", "Sup & No En Mean", 
                     "Sup & No En SD",
                     "Sup & En N", "Sup & En Mean", "Sup & En SD",
                     "No Sup & No En N", "No Sup & No En Mean", "No Sup & No En SD",
                     "No Sup & En N", "No Sup & En Mean", "No Sup & En SD",
                     "d YesNo-YesYes", "d YesNo-NoNo", "d YesNo-NoYes",
                     "d YesYes-NoNo", "d YesYes-YNoYes",
                     "d NoNo-NoYes"))

# Estimating initial unbalance
unbalance_2x2 <- initial_unbalance(Data = MAGMA_sim_data,
                                   group = c("gifted_support", "enrichment"),
                                   covariates = covariates_2x2)
unbalance_2x2

# Estimating and plotting density overlap in gifted support & enrichment propensity score
Density_overlap(Data = MAGMA_sim_data,
                variable = "ps_2x2", 
                group = c("gifted_support", "enrichment"),
                variable_name = "Propensity Score",
                group_labels = c("No Support & No Enrichment",
                                 "No Support & Enrichment",
                                 "Support & No Enrichment",
                                 "Support & Enrichment"),
                group_name = "Gifted Support & Enrichment")

## ----standard_2x2_group_matching, eval = FALSE--------------------------------
#  MAGMA_sim_data_2x2 <- MAGMA(Data = MAGMA_sim_data,
#                              group = c("gifted_support", "enrichment"),
#                              dist = "ps_2x2",
#                              cores = 2)

## ----standard_2x2_group_matching_str------------------------------------------
str(MAGMA_sim_data_2x2)

## ----Balance_2x2_group_matching, eval = FALSE---------------------------------
#  Balance_2x2 <- Balance_MAGMA(Data = MAGMA_sim_data_2x2,
#                               group = c("gifted_support", "enrichment"),
#                               covariates = covariates_2x2,
#                               step = "step")

## ----Balance_2x2_group_matching_results, fig.height = 5, fig.width = 7.5------

# Balance criteria for 100 cases per group
Balance_100_2x2_criteria <- Balance_extract(Balance = Balance_2x2,
                                            samplesize = 100,
                                            effects = FALSE)
Balance_100_2x2_criteria

# Extracting pairwise effects for 100 cases per group
Balance_100_2x2_effects <- Balance_extract(Balance = Balance_2x2,
                                           samplesize = 100,
                                           effects = TRUE)
Balance_100_2x2_effects

# Plotting trend over increasing sample size
Plot_MAGMA(Balance = Balance_2x2,
           criterion = c("Pillai", "d_ratio", "mean_g", "Adj_d_ratio")) 

# Creating table
Table_MAGMA(Balance = Balance_2x2,
            filename = "Balance_2x2.docx")

# Computing descriptive statistics and all pairwise effects after matching
descs_2x2_post <- MAGMA_desc(Data = MAGMA_sim_data_2x2,
                             group = c("gifted_support", "enrichment"),
                             covariates = covariates_2x2,
                             step_num = 100,
                             step_var = "step",
                             filename = "stats_post_pre.docx")

descs_2x2_post %>%
  purrr::set_names(c("Overall N", "Overall Mean", "Overall SD",
                     "Sup & No En N", "Sup & No En Mean", 
                     "Sup & No En SD",
                     "Sup & En N", "Sup & En Mean", "Sup & En SD",
                     "No Sup & No En N", "No Sup & No En Mean", "No Sup & No En SD",
                     "No Sup & En N", "No Sup & En Mean", "No Sup & En SD",
                     "d YesNo-YesYes", "d YesNo-NoNo", "d YesNo-NoYes",
                     "d YesYes-NoNo", "d YesYes-YNoYes",
                     "d NoNo-NoYes"))


## ----exact_2x2_group_matching, eval = FALSE-----------------------------------
#  MAGMA_sim_data_2x2_exact <- MAGMA_exact(Data = MAGMA_sim_data,
#                                          group = c("gifted_support", "enrichment"),
#                                          dist = "ps_2x2",
#                                          exact = "teacher_ability_rating",
#                                          cores = 2)

## ----exact_2x2_group_matching_str---------------------------------------------
str(MAGMA_sim_data_2x2_exact)


## ----exact_2x2_group_matching_balance, eval = FALSE---------------------------
#  # Estimating Balance
#  Balance_2x2_exact <- Balance_MAGMA(Data = MAGMA_sim_data_2x2_exact,
#                                     group = c("gifted_support", "enrichment"),
#                                     covariates = covariates_2x2,
#                                     step = "step") #Not necessary to define here

## ----exact_2x2_group_matching_results, fig.height = 5, fig.width = 7.5--------
# Balance criteria for 100 cases per group
Balance_100_2x2_criteria_exact <- Balance_extract(Balance = Balance_2x2_exact,
                                                  samplesize = 100,
                                                  effects = FALSE)
Balance_100_2x2_criteria_exact

# Extracting pairwise effects for 100 cases per group
Balance_100_2x2_effects_exact <- Balance_extract(Balance = Balance_2x2_exact,
                                                 samplesize = 100,
                                                 effects = TRUE)
Balance_100_2x2_effects_exact

# Plotting trend over increasing sample size
Plot_MAGMA(Balance = Balance_2x2_exact,
           criterion = c("Pillai", "d_ratio", "mean_g", "Adj_d_ratio")) #Could be omitted

# Creating table
Table_MAGMA(Balance = Balance_2x2_exact,
            filename = "Balance_2x2_exact.docx")

# Computing descriptive statistics and all pairwise effects post matching
descs_2x2_post <- MAGMA_desc(Data = MAGMA_sim_data_2x2_exact,
                             group = c("gifted_support", "enrichment"),
                             covariates = covariates_2x2,
                             step_num = 100,
                             step_var = "step",
                             filename = "stats_2x2_post.docx")

descs_2x2_post %>%
  purrr::set_names(c("Overall N", "Overall Mean", "Overall SD",
                     "Sup & No En N", "Sup & No En Mean", 
                     "Sup & No En SD",
                     "Sup & En N", "Sup & En Mean", "Sup & En SD",
                     "No Sup & No En N", "No Sup & No En Mean", "No Sup & No En SD",
                     "No Sup & En N", "No Sup & En Mean", "No Sup & En SD",
                     "d YesNo-YesYes", "d YesNo-NoNo", "d YesNo-NoYes",
                     "d YesYes-NoNo", "d YesYes-YNoYes",
                     "d NoNo-NoYes"))


