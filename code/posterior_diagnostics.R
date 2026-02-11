library(bayesplot)

# 提取后验样本
posterior_ct_fit <- as.array(ct_fit)

# all and core parameter names
all_params <- dimnames(posterior_ct_fit)$parameters
core_params <- all_params[!grepl("mu\\[|process_sd\\[", all_params) & all_params != "lp__"]
#print(core_params)


# 模型参数在当前值下计算出的未归一化的对数化的posterior probability
lp_ct_fit <- log_posterior(ct_fit)
head(lp_ct_fit)

# sampler internal diagnostic criterion
nuts_ct_fit <- nuts_params(ct_fit)
head(nuts_ct_fit)

######### 
# acceptance rate line chart: 0.7 - 0.9 good
accept_stats <- nuts_ct_fit[nuts_ct_fit$Parameter == "accept_stat__", ]

ggplot(accept_stats, aes(x = Iteration, y = Value)) +
  geom_line(color = "steelblue", alpha = 0.8) +
  geom_hline(yintercept = c(0.7, 0.9), linetype = "dashed", color = "red") +
  facet_wrap(~ Chain, ncol = 1) + 
  labs(title = "Acceptance Statistics by Chain",
       subtitle = "Dashed lines show ideal range (0.7-0.9)",
       x = "Iteration",
       y = "accept_stat__") +
  theme_bw()


######### 
# mcmc_parcoord: each line represents one iteration
# population parameter有outlier，导致有尖峰
color_scheme_set("darkgray")
#p <- mcmc_parcoord(posterior_ct_fit[100:150, , ], pars = core_params, np = nuts_ct_fit)
p <- mcmc_parcoord(posterior_ct_fit, pars = core_params, np = nuts_ct_fit)
p + theme(
  axis.text.x = element_text(
    size = 8,           
    angle = 45,    
    hjust = 1,     
    vjust = 1           
  )
)


######### 
# mcmc_pairs
# removing dpsd,wrsd,wpsd outlier
pop_vars <- c("dpmeanW","wpmeanW","wrmeanW","dpsd","wpsd","wrsd")
sd_vars <- c("dpsd", "wpsd", "wrsd")
n <- n_distinct(refined_data$PersonIDClean)
dp_vars <- c(paste0("dp[", 1:n, "]"))
wp_vars <- c(paste0("wp[", 1:n, "]"))
wr_vars <- c(paste0("wr[", 1:n, "]"))

# param names are in dimension 3
param_names <- dimnames(posterior_ct_fit)[[3]]
sd_idx <- match(sd_vars, param_names)
n_iter <- dim(posterior_ct_fit)[1]
n_chain <- dim(posterior_ct_fit)[2]

remove_outliers <- function(x) {
  q1  <- quantile(x, 0.25)
  q3  <- quantile(x, 0.75)
  iqr <- q3 - q1
  lower <- q1 - 5 * iqr    
  upper <- q3 + 5 * iqr
  x >= lower & x <= upper
}


keep_list <- vector("list", n_chain)
for (ch in 1:n_chain) {
  keep_idx <- rep(TRUE, n_iter)
  for (sdv in sd_idx) {
    values <- posterior_ct_fit[, ch, sdv]
    keep_idx <- keep_idx & remove_outliers(values)
  }
  keep_list[[ch]] <- keep_idx
}

# mutual iterations within four chains after screening out the outliers
len_per_chain <- sapply(keep_list, sum)
min_len <- min(len_per_chain)
min_len


set.seed(123)
posterior_clean <- array(NA,
                         dim = c(min_len, n_chain, length(param_names)),
                         dimnames = dimnames(posterior_ct_fit))
# 采用随机抽样使得所有链保留相同迭代数
# 用于存储每个链的抽样位置
sampled_indices_list <- vector("list", n_chain)

for (ch in 1:n_chain) {
  # 当前链保留的 iteration 位置
  keep_positions <- which(keep_list[[ch]])  
  
  # 随机抽取 min_len 个位置
  if (length(keep_positions) > min_len) {
    idx <- sample(keep_positions, min_len, replace = FALSE)
  } else {
    idx <- keep_positions
  }
  
  idx <- sort(idx)
  sampled_indices_list[[ch]] <- idx
  
  # 填充 posterior_clean
  posterior_clean[, ch, ] <- posterior_ct_fit[idx, ch, ]
}

# transform nuts_ct_fit into wide format
nuts_filtered_list <- list()

for (ch in 1:n_chain) {
  tmp <- nuts_ct_fit %>%
    filter(Chain == ch)
  
  # 使用与上面相同的抽样位置
  keep_idx <- sampled_indices_list[[ch]]
  
  # 使用 slice 按原始位置索引选择行
  nuts_filtered_list[[ch]] <- tmp %>%
    slice(keep_idx)
}

nuts_clean <- bind_rows(nuts_filtered_list)

# defined function of mcmc_pair for better comparing different groups of variables
group_mcmc_pairs <- function (variables) {
  mcmc_pairs(
  posterior_clean,
  np = nuts_clean,
  pars = variables,
  off_diag_args = list(size = 0.75)
)
}

group_mcmc_pairs(pop_vars)
# pair 2: dpmeanW, dpsd, dp[1:6]
group_mcmc_pairs(c(pop_vars[1],sd_vars[1], dp_vars[1:6]))
# dpmeanW, dpsd, dp[7:12]
group_mcmc_pairs(c(pop_vars[1],sd_vars[1], dp_vars[7:12]))
# pair 3: wpmeanW, wpsd, wp[1:6]
group_mcmc_pairs(c(pop_vars[2],sd_vars[2], wp_vars[1:6]))
# wpmeanW, wpsd, wp[7:12]
group_mcmc_pairs(c(pop_vars[2],sd_vars[2], wp_vars[7:12]))
# pair 4: wrmeanW, wrsd, wr[1:6]
group_mcmc_pairs(c(pop_vars[3],sd_vars[3], wr_vars[1:6]))
# wrmeanW, wrsd, wr[7:12]
group_mcmc_pairs(c(pop_vars[3],sd_vars[3], wr_vars[7:12]))



######### 全部用posterior_clean和nuts_clean?
# mcmc_trace
color_scheme_set("mix-brightblue-gray")
mcmc_trace(posterior_ct_fit, pars = "wpsd", np = nuts_ct_fit) + 
  xlab("Post-warmup iteration")

color_scheme_set("mix-brightblue-gray")
mcmc_trace(posterior_ct_fit, pars = "wpmeanW", np = nuts_ct_fit) + 
  xlab("Post-warmup iteration")


######### 代码不成功
# mcmc_nuts_divergence
color_scheme_set("red")
mcmc_nuts_divergence(nuts_ct_fit,posterior_ct_fit)


######### 
# mcmc_nuts_energy
color_scheme_set("red")
mcmc_nuts_energy(nuts_ct_fit)


######### 
# extract R-hat
summ <- summary(ct_fit)$summary
summ[, c("n_eff", "Rhat")]
summary(summ)

# mcmc_rhat
rhats <- rhat(ct_fit)
core_rhats <- rhats[grepl("^(dp|wp|wr|tp)", names(rhats))]
print(core_rhats)

color_scheme_set("brightblue") # see help("color_scheme_set")
mcmc_rhat(core_rhats) + yaxis_text(hjust = 1)

# mcmc_neff
ratios <- neff_ratio(ct_fit)   # 生成有效样本比例
core_ratios <- ratios[grepl("^(dp|wp|wr|tp)", names(ratios))] 
print(core_ratios)
mcmc_neff(core_ratios, size = 2) + yaxis_text(hjust = 1)

######### 
# mcmc_areas
# population-level parameter
mcmc_areas(posterior,
           pars = c("dpmeanW", "wpmeanW", "wrmeanW", 
                    "dpsd", "wpsd", "wrsd"),
           prob = 0.8) +
  ggtitle("Population-level Posterior Distributions (80% Credible Intervals)") +
  scale_x_continuous(limits = c(0, 90)) +   
  theme_minimal(base_size = 14)


# individual-level
mcmc_areas(posterior, regex_pars = "dp\\[", prob = 0.8) +
  ggtitle("Individual-level dp Posterior Distributions")

mcmc_areas(posterior, regex_pars = "wp\\[", prob = 0.8) +
  ggtitle("Individual-level wp Posterior Distributions")

mcmc_areas(posterior, regex_pars = "wr\\[", prob = 0.8) +
  ggtitle("Individual-level wr Posterior Distributions")

mcmc_areas(posterior, regex_pars = "tp\\[", prob = 0.8) +
  ggtitle("Individual-level tp Posterior Distributions")