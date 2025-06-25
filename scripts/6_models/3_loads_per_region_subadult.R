##### Bayesian models for LMS: total load #####

pacman::p_load(bayesplot, brms, dplyr, data.table)

args <- commandArgs(trailingOnly = TRUE)
method <- args[[1]]
region <- args[[2]]
out <- args[[3]]
iter <- args[[4]]
burn <- args[[5]]
thin <- args[[6]]

# load loads
load(file = "output/loads_per_region.RData")

# load meta
load(file = "data/metadata/metadata_adult_chick.RData")

load_per_region <- left_join(load_per_region, meta, by = "id")
load_per_region <- load_per_region %>% mutate(age = as.factor(case_when(
  grepl("C", load_per_region$id) ~ "chick",
  grepl("D", load_per_region$id) ~ "adult"
)))

load_per_region$age <- factor(load_per_region$age, levels = c("adult", "chick"))

# load survival data (sub)-adults
load("data/phenotypic/phenotypes_lifetime.RData")
load_per_region <- left_join(load_per_region, pheno_wide[,c("id", "lifespan")])

load_per_region <- load_per_region %>% mutate(lifespan_cat = as.factor(case_when(
  lifespan == 1 ~ "Yearling",
  lifespan > 1 ~"Adult"
)))
load_per_region$lifespan_cat <- factor(load_per_region$lifespan_cat, levels = c("Yearling", "Adult"))

# subset only the relevant method/loadtype
type = paste0(method, "_", region)

#### model ####
fit <- brm(scale(total_load) ~ lifespan_cat + (1|site), data = subset(load_per_region, loadtype == type),
           family = "gaussian",
           prior = prior(normal(0,1), class = b),
           cores =8, control = list(adapt_delta = 0.99, max_treedepth = 15),
           iter = iter, thin = thin, warmup = burn, seed = 1908)

save(fit, file = out)


