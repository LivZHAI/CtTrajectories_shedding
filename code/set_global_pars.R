# global_pars <- c(
# 	lod = 40							# Limit of detection
# 	)

global_pars <- c(
lod = ifelse((names(refined_data)[3] == 'CtValue'), 40, min(log(refined_data[3])))
)