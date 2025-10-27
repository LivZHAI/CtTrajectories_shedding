library(tidyverse) 
library(scales)

# 1. load self-defined functions and global para
source('code/utils.R')
source("code/set_global_pars.R")

# 2. read data (already read in create_function.R)
refined_data

# 3. data analysis and data fitting
source("code/fit_posteriors_preamble.R")
source("code/fit_posteriors.R")

# 4. generate outcomes of figures and distribution summary
source("code/make_figures.R")
source("code/summarise_dists.R")

# 5. save the results
savedir <- paste0("figures/")
source("code/save_figures.R")
save(ct_fit, file="output/ct_fit.RData")