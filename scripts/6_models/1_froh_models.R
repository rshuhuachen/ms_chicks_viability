#### Load packages ####
pacman::p_load(brms, bayesplot, tidyverse, data.table)

#### Load Froh data ####
load(file = "output/2_inbreeding/froh_chick_adult.RData")

### Load metadata ###
load(file = "data/metadata/metadata_adult_chick.RData")

## merge
froh <- left_join(froh, meta, by = "id")
froh$age <- factor(froh$age, levels = c("chick", "adult"))

### brms model chick vs adult ###
## parameters
iter = 1000000
burn = 500000
thin = 1000

## total
summary(lmerTest::lmer(scale(froh) ~ age + (1|site), data = subset(froh, id != "C09")))

brm_froh_chick <- brm(scale(froh) ~ age + (1|site), data = subset(froh, id != "C09"),
                      family = "gaussian",
                      prior = prior(normal(0,1), class = b),
                      cores =8, control = list(adapt_delta = 0.99, max_treedepth = 15),
                      iter = iter, thin = thin, warmup = burn, seed = 1908)

save(brm_froh_chick, file = "output/5_models/brm_froh_chicks.RData")

