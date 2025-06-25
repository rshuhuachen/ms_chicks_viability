#### Load packages ####
pacman::p_load(brms, bayesplot, tidyverse, data.table, lmerTest)

#### Load Froh data ####
load(file = "output/froh_chick_adult.RData")

### Load metadata ###
load(file = "data/metadata/metadata_adult_chick.RData")

## merge
froh <- left_join(froh, meta, by = "id")
froh$age <- factor(froh$age, levels = c("adult", "chick"))

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

### brms model yearling vs adult ####
load("data/phenotypic/phenotypes_lifetime.RData")
froh <- left_join(froh, pheno_wide[,c("id", "lifespan")], by = "id")
froh <- froh %>% mutate(lifespan_cat = as.factor(case_when(
  lifespan == 1 ~ "yearling",
  lifespan > 1 ~ "adult"
)))
froh$lifespan_cat <- factor(froh$lifespan_cat, levels = c("adult", "yearling"))

summary(lmerTest::lmer(scale(froh) ~ lifespan_cat + (1|site),
                       data = froh))


brm_froh_yearling <- brm(scale(froh) ~ lifespan_cat + (1|site), data = subset(froh, id != "C09"),
                      family = "gaussian",
                      prior = prior(normal(0,1), class = b),
                      cores =8, control = list(adapt_delta = 0.99, max_treedepth = 15),
                      iter = iter, thin = thin, warmup = burn, seed = 1908)

save(brm_froh_yearling, file = "output/5_models/brm_froh_yearling.RData")



# rerun model early life without single most inbred chick


brm_froh_chick_minus1 <- brm(scale(froh) ~ age + (1|site), data = subset(froh, id != "C09" & id != "C37"),
                      family = "gaussian",
                      prior = prior(normal(0,1), class = b),
                      cores =8, control = list(adapt_delta = 0.99, max_treedepth = 15),
                      iter = iter, thin = thin, warmup = burn, seed = 1908)

save(brm_froh_chick_minus1, file = "output/5_models/brm_froh_chicks_minus1.RData")

summary(brm_froh_chick_minus1)
