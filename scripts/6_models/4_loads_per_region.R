##### Bayesian models for LMS: total load #####

pacman::p_load(bayesplot, brms, dplyr, data.table)

args <- commandArgs(trailingOnly = TRUE)
method <- args[[1]]
region <- args[[2]]
out <- args[[3]]
iter <- args[[4]]
warm <- args[[5]]
thin <- args[[6]]

# load loads
load(file = "output/4_load/loads_per_region.RData")

# subset only the relevant method/loadtype
type == paste0(method, "_", region)

#### model ####
fit <- brm(scale(total_load) ~ age + (1|site), data = subset(loads, loadtype == type),
           family = "gaussian",
           prior = prior(normal(0,1), class = b),
           cores =8, control = list(adapt_delta = 0.99, max_treedepth = 15),
           iter = iter, thin = thin, warmup = burn, seed = 1908)

save(fit, file = out)


