### Here we run bayesian models to compare inbreeding and mutation loads between chicks and adults ####
pacman::p_load(brms, tidyverse, bayesplot, ggridges)
source("scripts/theme_ggplot.R")

## load data
load(file = "output/4_load/loads.RData")

## relevel
loads$age <- factor(loads$age, levels = c("chick", "adult"))

## parameters
iter = 1000000
burn = 500000
thin = 1000

## total
brm_gerp <- brm(scale(total_load) ~ age + (1|site), data = subset(loads, loadtype == "gerp"),
                  family = "gaussian",
                  prior = prior(normal(0,1), class = b),
                  cores =8, control = list(adapt_delta = 0.99, max_treedepth = 15),
                  iter = iter, thin = thin, warmup = burn, seed = 1908)

save(brm_gerp, file = "output/5_models/brm_gerp_total_chicks.RData")

brm_high <- brm(scale(total_load) ~ age + (1|site), data = subset(loads, loadtype == "high"),
                family = "gaussian",
                prior = prior(normal(0,1), class = b),
                cores =8, control = list(adapt_delta = 0.99, max_treedepth = 15),
                iter = iter, thin = thin, warmup = burn, seed = 1908)

save(brm_high, file = "output/5_models/brm_high_total_chicks.RData")

## hom
brm_gerp_hom <- brm(scale(hom_load) ~ age + (1|site), data = subset(loads, loadtype == "gerp"),
                family = "gaussian",
                prior = prior(normal(0,1), class = b),
                cores =8, control = list(adapt_delta = 0.99, max_treedepth = 15),
                iter = iter, thin = thin, warmup = burn, seed = 1908)

save(brm_gerp_hom, file = "output/5_models/brm_gerp_hom_chicks.RData")

brm_high_hom <- brm(scale(hom_load) ~ age + (1|site), data = subset(loads, loadtype == "high"),
                family = "gaussian",
                prior = prior(normal(0,1), class = b),
                cores =8, control = list(adapt_delta = 0.99, max_treedepth = 15),
                iter = iter, thin = thin, warmup = burn, seed = 1908)

save(brm_high_hom, file = "output/5_models/brm_high_hom_chicks.RData")

## het
brm_gerp_het <- brm(scale(het_load) ~ age + (1|site), data = subset(loads, loadtype == "gerp"),
                family = "gaussian",
                prior = prior(normal(0,1), class = b),
                cores =8, control = list(adapt_delta = 0.99, max_treedepth = 15),
                iter = iter, thin = thin, warmup = burn, seed = 1908)

save(brm_gerp_het, file = "output/5_models/brm_gerp_het_chicks.RData")

brm_high_het <- brm(scale(het_load) ~ age + (1|site), data = subset(loads, loadtype == "high"),
                family = "gaussian",
                prior = prior(normal(0,1), class = b),
                cores =8, control = list(adapt_delta = 0.99, max_treedepth = 15),
                iter = iter, thin = thin, warmup = burn, seed = 1908)

save(brm_high_het, file = "output/5_models/brm_high_het_chicks.RData")


