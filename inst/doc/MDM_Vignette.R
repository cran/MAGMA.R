## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  warning = FALSE)
options(tibble.width = Inf)

## ----setup, include = FALSE---------------------------------------------------
library(MAGMA.R)
MAGMA_sim_data_gifted <- readRDS(file = "data_gifted_matched_MDM.rds")
MAGMA_sim_data_gifted_exact <- readRDS(file ="data_gifted_matched_exact_MDM.rds")
MAGMA_sim_data_tar <- readRDS(file = "data_tar_matched_MDM.rds")
MAGMA_sim_data_tar_exact <- readRDS(file = "data_tar_matched_exact_MDM.rds")
MAGMA_sim_data_2x2 <- readRDS(file = "data_2x2_matched_MDM.rds")
MAGMA_sim_data_2x2_exact <- readRDS(file = "data_2x2_matched_exact_MDM.rds")

Balance_gifted <- readRDS(file = "Balance_gifted_MDM.rds")
Balance_gifted_exact <- readRDS(file = "Balance_gifted_exact_MDM.rds")
Balance_tar <- readRDS(file = "Balance_tar_MDM.rds")
Balance_tar_exact<- readRDS(file = "Balance_tar_exact_MDM.rds")
Balance_2x2 <- readRDS(file = "Balance_2x2_MDM.rds")
Balance_2x2_exact <- readRDS(file = "Balance_2x2_exact_MDM.rds")


## ----covariates---------------------------------------------------------------
covariates <- c("GPA_school",
                       "IQ_score",
                       "Motivation")

## ----unbalance_gifted---------------------------------------------------------
# Estimate overall and group specific descriptive statistics and Cohen’s d 
descs_gifted_pre <- MAGMA_desc(Data = MAGMA_sim_data,
                               group = "gifted_support",
                               covariates = covariates,
                               filename = "stats_gifted_pre.docx")
descs_gifted_pre %>%
  purrr::set_names(c("Overall N", "Overall Mean", "Overall SD",
                     "No Support N", "No Support Mean", "No Support SD",
                     "Support N", "Support Mean", "Support SD",
                     "d"))

# Estimating the four balance criteria
unbalance_gifted <- initial_unbalance(Data = MAGMA_sim_data,
                                      group = "gifted_support",
                                      covariates = covariates,
                                      round = 3)
unbalance_gifted

## ----standard_2_group_matching, eval = FALSE----------------------------------
# # Conducting matching for gifted support
# MAGMA_sim_data_gifted <- MAGMA(Data = MAGMA_sim_data,
#                                group = "gifted_support",
#                                covs = covariates,
#                                weights = c(1, 1.5, 0.5),
#                                cores = 2,
#                                verbose = TRUE)
# 

## ----Balance_standard_2_group_matching, eval = FALSE--------------------------
# # Estimating the four balance criteria iteratively over possible sample sizes
# Balance_gifted <- Balance_MAGMA(Data = MAGMA_sim_data_gifted,
#                                 group = "gifted_support",
#                                 covariates = covariates,
#                                 step = "step")

## ----Balance_standard_2_group_matching_out,fig.width=6, fig.height=4----------
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

# Plotting balance trend over sample size
Plot_MAGMA(Balance = Balance_gifted,
           criterion = c("Pillai", "d_ratio", "mean_g", "Adj_d_ratio"),
           print = FALSE)

# Sample Sizes with optimized balance criteria
Table_MAGMA(Balance = Balance_gifted,
            round = 3)

# Computing descriptive statistics and pairwise effects for 100 cases per group
descs_gifted_post <- MAGMA_desc(Data = MAGMA_sim_data_gifted,
                                group = "gifted_support",
                                covariates = covariates,
                                step_num = 255,
                                step_var = "step")

# Displaying the table with defined colum names
descs_gifted_post %>%
  purrr::set_names(c("Overall N", "Overall Mean", "Overall SD",
                     "No Support N", "No Support Mean", "No Support SD",
                     "Support N", "Support Mean", "Support SD",
                     "d"))

## ----exact_2_group_matching, eval = FALSE-------------------------------------
# # Conducting matching for gifted support
# MAGMA_sim_data_gifted_exact <- MAGMA_exact(Data = MAGMA_sim_data,
#                                            group = "gifted_support",
#                                            covs = covariates,
#                                            weights = c(1, 1.5, 0.5),
#                                            cores = 2,
#                                            exact = "gender")
# 

## ----Balance_exact_2_group_matching, eval = FALSE-----------------------------
# # Estimating the four balance criteria iteratively over possible sample sizes
# Balance_gifted_exact <- Balance_MAGMA(Data = MAGMA_sim_data_gifted_exact,
#                                       group = "gifted_support",
#                                       covariates = covariates,
#                                       step = "step")
# 

## ----Balance_sexact2_group_matching_out,fig.width=6, fig.height=4-------------

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

# Plotting balance trend over sample size
Plot_MAGMA(Balance = Balance_gifted_exact,
           criterion = c("Pillai", "d_ratio", "mean_g", "Adj_d_ratio"),
           print = FALSE)

# Sample Sizes with optimized balance criteria
Table_MAGMA(Balance = Balance_gifted_exact,
            round = 3)

# Computing descriptive statistics and pairwise effects for 100 cases per group
descs_gifted_post_exact <- MAGMA_desc(Data = MAGMA_sim_data_gifted_exact,
                                      group = "gifted_support",
                                      covariates = covariates,
                                      step_num = 179,
                                      step_var = "step")

# Displaying the table with defined colum names
descs_gifted_post_exact %>%
  purrr::set_names(c("Overall N", "Overall Mean", "Overall SD",
                     "No Support N", "No Support Mean", "No Support SD",
                     "Support N", "Support Mean", "Support SD",
                     "d"))

## ----unbalance_tar------------------------------------------------------------
# Computing descriptive statistics and all pairwise effects for three groups
descs_tar_pre <- MAGMA_desc(Data = MAGMA_sim_data,
                            group = "teacher_ability_rating",
                            covariates = covariates)

descs_tar_pre %>%
  purrr::set_names(c("Overall N", "Overall Mean", "Overall SD",
                     "BA N", "BA Support Mean", "BA Support SD",
                     "A N", "A Mean", "A SD",
                     "AA N", "AA Mean", "AA SD",
                     "d BA-A", "d BA-AA", "d A-AA"))

# Estimating and printing initial unbalance for teacher rated ability
unbalance_tar <- initial_unbalance(Data = MAGMA_sim_data,
                                   group = "teacher_ability_rating",
                                   covariates = covariates,
                                   round = 3)
unbalance_tar

## ----standard_3_group_matching, eval = FALSE----------------------------------
# # Conducting matching for gifted support
# MAGMA_sim_data_tar <- MAGMA(Data = MAGMA_sim_data,
#                             group = "teacher_ability_rating",
#                             covs = covariates,
#                             weights = c(1, 1.5, 0.5),
#                             cores = 2)
# 

## ----Balance_standard_3_group_matching, eval = FALSE--------------------------
# # Estimating the four balance criteria iteratively over possible sample sizes
# Balance_tar <- Balance_MAGMA(Data = MAGMA_sim_data_tar,
#                                 group = "teacher_ability_rating",
#                                 covariates = covariates,
#                                 step = "step")

## ----Balance_standard_3_group_matching_out,fig.width=6, fig.height=4----------
# Extracting balance criteria for 100 cases per group
Balance_100_criteria_tar <- Balance_extract(Balance = Balance_tar,
                                            samplesize = 100,
                                            effects = FALSE)
Balance_100_criteria_tar

# Extracting pairwise effects for 100 cases per group
Balance_100_effects_tar <- Balance_extract(Balance = Balance_tar,
                                           samplesize = 100,
                                           effects = TRUE)
Balance_100_effects_tar

# Plotting balance trend over sample size
Plot_MAGMA(Balance = Balance_tar,
           criterion = c("Pillai", "d_ratio", "mean_g", "Adj_d_ratio"),
           print = FALSE)

# Sample Sizes with optimized balance criteria
Table_MAGMA(Balance = Balance_tar,
            round = 3)

# Computing descriptive statistics and pairwise effects for 100 cases per group
descs_tar_post <- MAGMA_desc(Data = MAGMA_sim_data_tar,
                             group = "teacher_ability_rating",
                             covariates = covariates,
                             step_num = 120,
                             step_var = "step")

# Displaying the table with defined colum names
descs_tar_post %>%
  purrr::set_names(c("Overall N", "Overall Mean", "Overall SD",
                     "BA N", "BA Support Mean", "BA Support SD",
                     "A N", "A Mean", "A SD",
                     "AA N", "AA Mean", "AA SD",
                     "d BA-A", "d BA-AA", "d A-AA"))

## ----exact_3_group_matching, eval = FALSE-------------------------------------
# # Conducting matching for gifted support
# MAGMA_sim_data_tar_exact <- MAGMA_exact(Data = MAGMA_sim_data,
#                                         group = "teacher_ability_rating",
#                                         covs = covariates,
#                                         weights = c(1, 1.5, 0.5),
#                                         cores = 2,
#                                         exact = "gender")
# 

## ----Balance_exact_3_group_matching, eval = FALSE-----------------------------
# # Estimating the four balance criteria iteratively over possible sample sizes
# Balance_tar_exact <- Balance_MAGMA(Data = MAGMA_sim_data_tar_exact,
#                                    group = "teacher_ability_rating",
#                                    covariates = covariates,
#                                    step = "step")

## ----Balance_exact_3_group_matching_out,fig.width=6, fig.height=4-------------
# Extracting balance criteria for 100 cases per group
Balance_100_criteria_tar_exact <- Balance_extract(Balance = Balance_tar_exact,
                                                  samplesize = 100,
                                                  effects = FALSE)
Balance_100_criteria_tar_exact

# Extracting pairwise effects for 100 cases per group
Balance_100_effects_tar_exact <- Balance_extract(Balance = Balance_tar_exact,
                                                 samplesize = 100,
                                                 effects = TRUE)
Balance_100_effects_tar_exact

# Plotting balance trend over sample size
Plot_MAGMA(Balance = Balance_tar_exact,
           criterion = c("Pillai", "d_ratio", "mean_g", "Adj_d_ratio"),
           print = FALSE)

# Sample Sizes with optimized balance criteria
Table_MAGMA(Balance = Balance_tar_exact,
            round = 3)

# Computing descriptive statistics and pairwise effects for 100 cases per group
descs_tar_post_exact <- MAGMA_desc(Data = MAGMA_sim_data_tar_exact,
                                   group = "teacher_ability_rating",
                                   covariates = covariates,
                                   step_num = 116,
                                   step_var = "step")

# Displaying the table with defined colum names
descs_tar_post_exact %>%
  purrr::set_names(c("Overall N", "Overall Mean", "Overall SD",
                     "BA N", "BA Support Mean", "BA Support SD",
                     "A N", "A Mean", "A SD",
                     "AA N", "AA Mean", "AA SD",
                     "d BA-A", "d BA-AA", "d A-AA"))

## ----unbalance_2x2------------------------------------------------------------
# Computing descriptive statistics and all pairwise effects for three groups
descs_2x2 <- MAGMA_desc(Data = MAGMA_sim_data,
                        group = c("gifted_support", "enrichment"),
                        covariates = covariates)

descs_2x2  %>%
  purrr::set_names(c("Overall N", "Overall Mean", "Overall SD",
                     "Sup & No En N", "Sup & No En Mean", 
                     "Sup & No En SD",
                     "Sup & En N", "Sup & En Mean", "Sup & En SD",
                     "No Sup & No En N", "No Sup & No En Mean", "No Sup & No En SD",
                     "No Sup & En N", "No Sup & En Mean", "No Sup & En SD",
                     "d YesNo-YesYes", "d YesNo-NoNo", "d YesNo-NoYes",
                     "d YesYes-NoNo", "d YesYes-YNoYes",
                     "d NoNo-NoYes"))

# Estimating and printing initial unbalance for teacher rated ability
unbalance_2x2 <- initial_unbalance(Data = MAGMA_sim_data,
                                   group = c("gifted_support", "enrichment"),
                                   covariates = covariates,
                                   round = 3)
unbalance_2x2

## ----standard_2x2_group_matching, eval = FALSE--------------------------------
# # Conducting matching for gifted support
# MAGMA_sim_data_2x2 <- MAGMA(Data = MAGMA_sim_data,
#                             group = c("gifted_support", "enrichment"),
#                             covs = covariates,
#                             weights = c(1, 1.5, 0.5),
#                             cores = 2)
# 

## ----Balance_standard_2x2_group_matching, eval = FALSE------------------------
# # Estimating the four balance criteria iteratively over possible sample sizes
# Balance_2x2 <- Balance_MAGMA(Data = MAGMA_sim_data_2x2,
#                                 group = c("gifted_support", "enrichment"),
#                                 covariates = covariates,
#                                 step = "step")

## ----Balance_standard_2x2_group_matching_out,fig.width=6, fig.height=4--------
# Extracting balance criteria for 100 cases per group
Balance_100_criteria_2x2 <- Balance_extract(Balance = Balance_2x2,
                                            samplesize = 100,
                                            effects = FALSE)
Balance_100_criteria_2x2

# Extracting pairwise effects for 100 cases per group
Balance_100_effects_2x2 <- Balance_extract(Balance = Balance_2x2,
                                           samplesize = 100,
                                           effects = TRUE)
Balance_100_effects_2x2

# Plotting balance trend over sample size
Plot_MAGMA(Balance = Balance_2x2,
           criterion = c("Pillai", "d_ratio", "mean_g", "Adj_d_ratio"),
           print = FALSE)

# Sample Sizes with optimized balance criteria
Table_MAGMA(Balance = Balance_2x2,
            round = 3)

# Computing descriptive statistics and pairwise effects for 100 cases per group
descs_2x2_post <- MAGMA_desc(Data = MAGMA_sim_data_2x2,
                             group = c("gifted_support", "enrichment"),
                             covariates = covariates,
                             step_num = 112,
                             step_var = "step")

# Displaying the table with defined colum names
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
# # Conducting matching for gifted support
# MAGMA_sim_data_2x2_exact <- MAGMA_exact(Data = MAGMA_sim_data,
#                                         group = c("gifted_support", "enrichment"),
#                                         covs = covariates,
#                                         weights = c(1, 1.5, 0.5),
#                                         cores = 2,
#                                         exact = "gender")
# 

## ----Balance_exact_2x2_group_matching, eval = FALSE---------------------------
# # Estimating the four balance criteria iteratively over possible sample sizes
# Balance_2x2_exact <- Balance_MAGMA(Data = MAGMA_sim_data_2x2_exact,
#                                    group = c("gifted_support", "enrichment"),
#                                    covariates = covariates,
#                                    step = "step")

## ----Balance_exact_2x2_group_matching_out,fig.width=6, fig.height=4-----------
# Extracting balance criteria for 100 cases per group
Balance_100_criteria_2x2_exact <- Balance_extract(Balance = Balance_2x2_exact,
                                                  samplesize = 100,
                                                  effects = FALSE)
Balance_100_criteria_2x2_exact

# Extracting pairwise effects for 100 cases per group
Balance_100_effects_2x2_exact <- Balance_extract(Balance = Balance_2x2_exact,
                                                 samplesize = 100,
                                                 effects = TRUE)
Balance_100_effects_2x2_exact

# Plotting balance trend over sample size
Plot_MAGMA(Balance = Balance_2x2_exact,
           criterion = c("Pillai", "d_ratio", "mean_g", "Adj_d_ratio"),
           print = FALSE)

# Sample Sizes with optimized balance criteria
Table_MAGMA(Balance = Balance_2x2_exact,
            round = 3)

# Computing descriptive statistics and pairwise effects for 100 cases per group
descs_2x2_post_exact <- MAGMA_desc(Data = MAGMA_sim_data_2x2_exact,
                                   group = c("gifted_support", "enrichment"),
                                   covariates = covariates,
                                   step_num = 116,
                                   step_var = "step")

# Displaying the table with defined colum names
descs_2x2_post_exact %>%
  purrr::set_names(c("Overall N", "Overall Mean", "Overall SD",
                     "Sup & No En N", "Sup & No En Mean", 
                     "Sup & No En SD",
                     "Sup & En N", "Sup & En Mean", "Sup & En SD",
                     "No Sup & No En N", "No Sup & No En Mean", "No Sup & No En SD",
                     "No Sup & En N", "No Sup & En Mean", "No Sup & En SD",
                     "d YesNo-YesYes", "d YesNo-NoNo", "d YesNo-NoYes",
                     "d YesYes-NoNo", "d YesYes-YNoYes",
                     "d NoNo-NoYes"))

