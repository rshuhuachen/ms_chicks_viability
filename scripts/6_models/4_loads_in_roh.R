##### here we isolate the load within ROHs and model their fitness effects #####

#### load packages ####
extrafont::loadfonts(device="all")
pacman::p_load(tidyverse, data.table, lme4)

### load theme
source("scripts/theme_ggplot.R")

###  load ROHs ####
roh <- read.csv("output/2_inbreeding/bcftools_roh_chick_adult.txt")
# exclude sex scaf and poor quality ones 
scaf <- fread("data/metadata/30_largest.scafs.tsv")
scaf$scaf <- gsub(":", ";", scaf$scaf)
scaf$scaf <- gsub("\\.", "=", scaf$scaf)

rohs <- subset(roh, qual > 30 & length >= 100000 & nsnp >= 100 & chr != scaf$scaf[which(scaf$scaf_no==4)]) 

rohs <- rohs %>% mutate(roh_nr = row_number())
roh_list <- rohs %>% group_split(roh_nr)
test <- roh_list[1:10]

#### GERP #### 
 
load(file = "output/4_load/gerp/gerps_all.RData")

## recode mutations
gt <- c(11:ncol(gerp))
select_n3 <- function(x){x = substr(x,1,3)}
gerp[gt] <- lapply(gerp[gt], select_n3)

## make it long
gerps_long <- gather(gerp, id, genotype, C01:D229777, factor_key = T)
gerps_long <- gerps_long %>% mutate(zygosity = case_when(
  genotype == "0/0" ~ "hom_anc",
  genotype == "0/1" | genotype == "1/0" ~ "het",
  genotype == "1/1" ~ "hom_der",
  TRUE ~ NA
))

gerps_long_hom_der <- subset(gerps_long, zygosity == "hom_der") #7636432

### count how many gerp mutations per ROH

# GERP table (assuming columns id, chr, pos)
setDT(gerps_long_hom_der)
gerps_long_hom_der[, start := pos]
gerps_long_hom_der[, end := pos]

# ROHs need start/end columns too
rohs_all <- rohs
setDT(rohs_all)

# set keys for overlap join
setkey(rohs_all, id, chr, start, end)
setkey(gerps_long_hom_der, id, chr, start, end)

# do overlap join
overlaps <- foverlaps(
  gerps_long_hom_der, rohs_all,
  by.x = c("id", "chr", "start", "end"),
  by.y = c("id", "chr", "start", "end"),
  nomatch = 0L
)

# count mutations per ROH
counts <- overlaps[, .N, by = .(id, chr, start, end, length)]

# merge back
gerp_per_roh <- merge(rohs_all, counts,
                      by = c("id", "chr", "start", "end", "length"),
                      all.x = TRUE)

gerp_per_roh[is.na(N), N := 0]

# sum number of hom mutations in ROH to get a homozygous FROH load

froh_load_gerp <- gerp_per_roh %>% group_by(id) %>% summarise(froh_load_gerp = sum(N))

#### SnpEff #####

high <- read.table("output/4_load/snpeff/ltet_filtered_ann_aa_chick_high.vcf.gz")
#rename
names <- fread("output/2_inbreeding/chicks_adults_samples.txt", header = F)
names(names) <- "file_id"
ids <- read.csv("data/metadata/file_list_all_bgi_clean.csv")
ids$file_id <- gsub(".sorted.bam", "", ids$file)
names <- left_join(names, ids[,c("file_id", "id")], by = c("file_id"))
names <- names %>% mutate(id = case_when(
  is.na(id) ~ file_id,
  TRUE ~ id
))

names(high)[1:9] <- c("chr", "pos", "ID", "ref", "alt", "QUAL", "FILTER", "INFO","FORMAT")
names(high)[c(10:ncol(high))] <- names$id

## recode mutations
gt <- c(10:ncol(high))
select_n3 <- function(x){x = substr(x,1,3)}
high[gt] <- lapply(high[gt], select_n3)

## make it long
high_long <- gather(high, id, genotype, C01:D229777, factor_key = T)
high_long <- high_long %>% mutate(zygosity = case_when(
  genotype == "0/0" ~ "hom_anc",
  genotype == "0/1" | genotype == "1/0" ~ "het",
  genotype == "1/1" ~ "hom_der",
  TRUE ~ NA
))

high_long_hom_der <- subset(high_long, zygosity == "hom_der")

### count how many high mutations per ROH

# high table (assuming columns id, chr, pos)
setDT(high_long_hom_der)
high_long_hom_der[, start := pos]
high_long_hom_der[, end := pos]

# set keys for overlap join
setkey(rohs_all, id, chr, start, end)
setkey(high_long_hom_der, id, chr, start, end)

# do overlap join
overlaps <- foverlaps(
  high_long_hom_der, rohs_all,
  by.x = c("id", "chr", "start", "end"),
  by.y = c("id", "chr", "start", "end"),
  nomatch = 0L
)

# count mutations per ROH
counts <- overlaps[, .N, by = .(id, chr, start, end, length)]

# merge back
high_per_roh <- merge(rohs_all, counts,
                      by = c("id", "chr", "start", "end", "length"),
                      all.x = TRUE)

high_per_roh[is.na(N), N := 0]

# sum number of hom mutations in ROH to get a homozygous FROH load

froh_load_high <- high_per_roh %>% group_by(id) %>% summarise(froh_load_high = sum(N))

#### combine ####
# link it to froh
load("output/froh_chick_adult.RData")

froh_load <- left_join(froh[,c("id", "froh", "age")], froh_load_gerp, by = "id")
froh_load <- left_join(froh_load, froh_load_high, by = "id")

ggplot(froh_load, aes(x = froh, y = froh_load_gerp)) + geom_point(size=3, col = "black", fill = clr_high, shape = 21) + 
  geom_smooth(method='lm', col = clr_high) + 
  labs(x = expression(F[ROH]), y = "# GERP mutations in ROH") -> load_in_roh_gerp

ggplot(froh_load, aes(x = froh, y = froh_load_high)) + geom_point(size=3, col = "black", fill = clr_high, shape = 21) + 
  geom_smooth(method='lm', col = clr_high) + 
  labs(x = expression(F[ROH]), y = "# SnpEff mutations in ROH")-> load_in_roh_high

plot_grid(load_in_roh_gerp, load_in_roh_high, ncol = 2, align = "hv", axis = "lb",
          labels = c("A", "B"), 
          label_fontface = "plain", label_size = 22) -> load_in_roh_plot

ggsave(load_in_roh_plot, file = "plots/load_in_ROH.png", width=12, height=10)

cor.test(froh_load$froh, froh_load$froh_load_gerp)
cor.test(froh_load$froh, froh_load$froh_load_high)
