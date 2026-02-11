library(scales)
library(tidyverse) # include ggplot2, dplyr, tidyr, purrr, readr

# 1. import data
setwd("/Users/livz/Desktop/CtTrajectories_shedding")
source("code/data_import.R")

# raw_data <- read.csv("data/ct_dat_refined.csv", stringsAsFactors = FALSE)
# refined_data <- raw_data |>
#   dplyr::select(-all_of(c('PersonID','NovelPersistent','B117Status'))) |>
#   select(c('PersonIDClean'),  
#          everything())


# 2. set global parameter
source("code/set_global_pars.R")

# 3. load self-defined functions
source('code/utils.R')

# 4. data analysis and data fitting
source("code/fit_posteriors_preamble.R")
source("code/fit_posteriors.R")

# 5. generate outcomes of figures and distribution summary
source("code/make_figures.R")
source("code/summarise_dists.R")

# 6. save the results
savedir <- paste0("figures/")
source("code/save_figures.R")
save(ct_fit, file="output/ct_fit.RData")