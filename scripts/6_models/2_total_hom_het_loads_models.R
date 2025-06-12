### Here we run bayesian models to compare inbreeding and mutation loads between chicks and adults ####
pacman::p_load(brms, tidyverse, bayesplot, ggridges)
source("scripts/theme_ggplot.R")

## load data
load(file = "output/loads.RData")

#### early life #####
## relevel
loads$age <- factor(loads$age, levels = c("adult", "chick"))

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


#### yearling vs adult ####
load("data/phenotypic/phenotypes_lifetime.RData")
loads <- left_join(loads, pheno_wide[,c("id", "lifespan")], by = "id")
loads <- loads %>% mutate(lifespan_cat = as.factor(case_when(
  lifespan == 1 ~ "yearling",
  lifespan > 1 ~ "adult"
)))
loads$lifespan_cat <- factor(loads$lifespan_cat, levels = c("adult", "yearling"))

## total
brm_gerp_yearling <- brm(scale(total_load) ~ lifespan_cat + (1|site), data = subset(loads, loadtype == "gerp"),
                family = "gaussian",
                prior = prior(normal(0,1), class = b),
                cores =8, control = list(adapt_delta = 0.99, max_treedepth = 15),
                iter = iter, thin = thin, warmup = burn, seed = 1908)

save(brm_gerp_yearling, file = "output/5_models/brm_gerp_total_yearlings.RData")

brm_high_yearling <- brm(scale(total_load) ~ lifespan_cat + (1|site), data = subset(loads, loadtype == "high"),
                family = "gaussian",
                prior = prior(normal(0,1), class = b),
                cores =8, control = list(adapt_delta = 0.99, max_treedepth = 15),
                iter = iter, thin = thin, warmup = burn, seed = 1908)

save(brm_high_yearling, file = "output/5_models/brm_high_total_yearlings.RData")

## hom
brm_gerp_hom_yearling <- brm(scale(hom_load) ~ lifespan_cat + (1|site), data = subset(loads, loadtype == "gerp"),
                    family = "gaussian",
                    prior = prior(normal(0,1), class = b),
                    cores =8, control = list(adapt_delta = 0.99, max_treedepth = 15),
                    iter = iter, thin = thin, warmup = burn, seed = 1908)

save(brm_gerp_hom_yearling, file = "output/5_models/brm_gerp_hom_yearlings.RData")

brm_high_hom_yearling <- brm(scale(hom_load) ~ lifespan_cat + (1|site), data = subset(loads, loadtype == "high"),
                    family = "gaussian",
                    prior = prior(normal(0,1), class = b),
                    cores =8, control = list(adapt_delta = 0.99, max_treedepth = 15),
                    iter = iter, thin = thin, warmup = burn, seed = 1908)

save(brm_high_hom_yearling, file = "output/5_models/brm_high_hom_yearlings.RData")

## het
brm_gerp_het_yearling <- brm(scale(het_load) ~ lifespan_cat + (1|site), data = subset(loads, loadtype == "gerp"),
                    family = "gaussian",
                    prior = prior(normal(0,1), class = b),
                    cores =8, control = list(adapt_delta = 0.99, max_treedepth = 15),
                    iter = iter, thin = thin, warmup = burn, seed = 1908)

save(brm_gerp_het_yearling, file = "output/5_models/brm_gerp_het_yearlings.RData")

brm_high_het_yearling <- brm(scale(het_load) ~ lifespan_cat + (1|site), data = subset(loads, loadtype == "high"),
                    family = "gaussian",
                    prior = prior(normal(0,1), class = b),
                    cores =8, control = list(adapt_delta = 0.99, max_treedepth = 15),
                    iter = iter, thin = thin, warmup = burn, seed = 1908)

save(brm_high_het_yearling, file = "output/5_models/brm_high_het_yearlings.RData")


