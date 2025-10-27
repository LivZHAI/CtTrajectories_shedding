# jackknife是一种敏感性分析
# 它的核心目的是：系统地评估移除任何一个B.1.1.7感染个体对模型整体结果（如病毒动力学参数估计）的影响。
# 这种分析旨在回答一个非常重要的敏感性问题：“我们的结论是否过分依赖于某一个特殊的个体？”
# 具体来说，研究者可以比较这7份 dist_summary_no....csv 文件。如果排除任何一个B.1.1.7个体后，B.1.1.7组的参数估计（如平均峰值病毒载量 peak.geml.B117_mean）没有发生剧烈变化，
# 那么就说明模型的发现是稳健的（robust），不依赖于单个数据点。

library(tidyverse) 
library(scales)

source('code/utils.R')
source("code/set_global_pars.R")

B117inds <- c(1368,1371,1374,1375,3229,4399,4447)

for(indexZ in 1:length(B117inds)){

	ct_dat_refined <- read_csv("data/ct_dat_refined.csv")

	ct_dat_refined <- ct_dat_refined %>% 
		filter(PersonID!=B117inds[indexZ]) %>%
		select(-PersonIDClean) %>%
		clean_person_id

	source("code/fit_posteriors_preamble.R")
	source("code/fit_posteriors.R")

	source("code/make_figures.R")
	source("code/summarise_dists.R")
	savedir <- paste0("figures/no",B117inds[indexZ],"/")
	source("code/save_figures.R")
	write.csv(dist_summary, file=paste0("output/dist_summary_no",B117inds[indexZ],".csv"), row.names=FALSE)
	save(ct_fit, file=paste0("output/ct_fit_no",B117inds[indexZ],".RData"))

}

