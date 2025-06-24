##### Supplementary descriptive figure ####
pacman::p_load(tidyverse, data.table, cowplot)
source("scripts/theme_ggplot.R")
extrafont::loadfonts(device="all")

### Sup fig 1: inbreeding ####
# load rohs
rohs <- fread("output/2_inbreeding/bcftools_roh_chick_adult.txt")

# load froh
load(file = "output/froh_chick_adult.RData")
froh$age <- gsub("chick", "Chick", froh$age)

#add longevity data 
load("data/phenotypic/phenotypes_lifetime.RData")
froh <- left_join(froh, pheno_wide[,c("id", "lifespan")], by = "id")
froh <- froh %>% mutate(age = as.factor(case_when(
  age == "Chick" ~ "Chick",
  lifespan == 1 ~ "Yearling",
  lifespan > 1 ~ "Adult"
)))

froh$age <- factor(froh$age, levels = c("Chick", "Yearling", "Adult"))
rohs <- left_join(rohs, froh[,c("id", "age")], by = c("id"))

### Sup fig 1a: froh distribution ####

ggplot(froh, aes(x = froh, fill = age)) +  geom_density() +
  labs(x = expression(italic(F)[ROH]), y = "Density", fill = "Age") +
  scale_fill_manual(values = c(alpha(clrs_hunting[2], 0.7),
                               alpha(clrs_hunting[3], 0.7), 
                               alpha(clrs_hunting[4], 0.7))) -> froh_dist

froh_dist

### Sup fig 1b: roh length distribution ####
num_roh_per_ind <- rohs %>% group_by(id, age.y) %>% tally() 

ggplot(num_roh_per_ind, aes(n)) +
  geom_histogram(aes(fill = age.y), col = "black", bins = 20,  position = "dodge") + 
  scale_fill_manual(values = c(alpha(clrs_hunting[2], 0.7),
                               alpha(clrs_hunting[3], 0.7), 
                               alpha(clrs_hunting[4], 0.7)))+
  ylab("Number of individuals") +labs(fill = "Age")+
  xlab("Number of ROHs")  -> fig_hist_nroh

fig_hist_nroh

#### Figure 1b: ROH classes #### 

rohs_class <- rohs %>% mutate(class = as.factor(case_when(
  length < 1*10^6 ~ "< 1 Mb",
  length >= 1*10^6 & length < 2*10^6 ~ "1 - 2Mb",
  length >= 2*10^6 ~ "> 2Mb"
)))

rohs_class$class <- factor(rohs_class$class, levels = c("", "< 1 Mb", "1 - 2Mb", "> 2Mb"))
rohs_class_tally <- rohs_class %>% group_by(id, age.y, class) %>% tally()

rohs_class_tally %>% filter(id != "C06") %>% group_by(class, age.y) %>%
  ggplot(aes(x = age.y, y = n)) + 
  geom_boxplot(outlier.shape=NA, size = 1)+ 
  geom_point(aes(col=age.y), position = "jitter", size = 2)+
  scale_color_manual(values = c(alpha(clrs_hunting[2], 0.7),
                               alpha(clrs_hunting[3], 0.7), 
                               alpha(clrs_hunting[4], 0.7)))+
  facet_wrap(~ class, scales="free")+
  theme(legend.position="none",
        panel.border = element_blank(),
        panel.grid = element_blank(),
        strip.background = element_blank())+
  labs(x = "ROH length", y = "Number of ROHs per individual", col = "Age class") -> fig_roh_lengths

fig_roh_lengths

summary(rohs_class_tally)
rohs_class_tally %>% filter(id != "C06") %>% group_by(class, age.y) %>% summarise(mean_n = mean(n))
rohs_class_tally %>% filter(id != "C06") %>% group_by(class) %>% summarise(n = sum(n))

#### Combine in one figure ####

cowplot::plot_grid(froh_dist, fig_hist_nroh,
                   ncol = 2, labels = c("a", "b"), label_fontface = "plain", label_size = 22,
                   align = "hv", axis = "lb") -> fig_froh_top

cowplot::plot_grid(fig_froh_top, fig_roh_lengths,
                   ncol = 1, labels = c("", "c"), label_fontface = "plain", label_size = 22) -> fig_froh

ggsave(fig_froh, file = "plots/sup_fig_inbreeding.png", width=16,height=12)


#### Fig 2: mutations ####
### b: snpeff distribution ####
# scaf info
# scaf <- fread("data/metadata/30_largest.scafs.tsv")
# scaf$scaf <- gsub(":", ";", scaf$scaf)
# scaf$scaf <- gsub("\\.", "=", scaf$scaf)
# scaf <- subset(scaf, scaf_no != 4)
# 
# vcf <- read.table("output/3_annotated_genome/ltet_filtered_ann_aa_chick_correct.vcf.gz")
# vcf_29 <- subset(vcf, V1 %in% scaf$scaf)
# vcf_29_nowarn <- subset(vcf_29, !grepl("WARNING", V8))
# 
# #subset regions
# high <- subset(vcf_29_nowarn, grepl("HIGH", vcf_29_nowarn$V8))
# mod <- subset(vcf_29_nowarn, grepl("MODERATE", vcf_29_nowarn$V8))
# low <- subset(vcf_29_nowarn, grepl("LOW", vcf_29_nowarn$V8))
# modifier <- subset(vcf_29_nowarn, grepl("MODIFIER", vcf_29_nowarn$V8))
# 
# sum <- data.frame(type = c("total_29scaf", "total_high", "total_mod", "total_low", "total_modifier",
#                            "downstream_gene_variant","five_prime_UTR_premature_start_codon_gain_variant",
#                            "five_prime_UTR_variant", "initiator_codon_variant",  
#                            "intron_variant", "LOF", "missense_variant" ,"NMD","splice_acceptor_variant"   , 
#                            "splice_donor_variant" ,  "splice_region_variant" , "start_lost" ,   
#                            "stop_gained" , "stop_lost",
#                            "stop_retained_variant" ,"synonymous_variant",
#                            "three_prime_UTR_variant", "upstream_gene_variant", 
#                            "intergenic_region"),
#                   n_mutations = c(nrow(vcf_29_nowarn), nrow(high), nrow(mod), nrow(low), nrow(modifier),
#                                   nrow(subset(vcf_29_nowarn, grepl("downstream_gene_variant", vcf_29_nowarn$V8))),
#                                   nrow(subset(vcf_29_nowarn, grepl("5_prime_UTR_premature_start_codon_gain_variant", vcf_29_nowarn$V8))),
#                                   nrow(subset(vcf_29_nowarn, grepl("5_prime_UTR_variant", vcf_29_nowarn$V8))),
#                                   nrow(subset(vcf_29_nowarn, grepl("initiator_codon_variant", vcf_29_nowarn$V8))),
#                                   nrow(subset(vcf_29_nowarn, grepl("intron_variant", vcf_29_nowarn$V8))),
#                                   nrow(subset(vcf_29_nowarn, grepl("LOF", vcf_29_nowarn$V8))),
#                                   nrow(subset(vcf_29_nowarn, grepl("missense_variant", vcf_29_nowarn$V8))),
#                                   nrow(subset(vcf_29_nowarn, grepl("NMD", vcf_29_nowarn$V8))),
#                                   nrow(subset(vcf_29_nowarn, grepl("splice_acceptor_variant", vcf_29_nowarn$V8))),
#                                   nrow(subset(vcf_29_nowarn, grepl("splice_donor_variant", vcf_29_nowarn$V8))),
#                                   nrow(subset(vcf_29_nowarn, grepl("splice_region_variant", vcf_29_nowarn$V8))),
#                                   nrow(subset(vcf_29_nowarn, grepl("start_lost", vcf_29_nowarn$V8))),
#                                   nrow(subset(vcf_29_nowarn, grepl("stop_gained", vcf_29_nowarn$V8))),
#                                   nrow(subset(vcf_29_nowarn, grepl("stop_lost", vcf_29_nowarn$V8))),
#                                   nrow(subset(vcf_29_nowarn, grepl("stop_retained_variant", vcf_29_nowarn$V8))),
#                                   nrow(subset(vcf_29_nowarn, grepl("synonymous_variant", vcf_29_nowarn$V8))),
#                                   nrow(subset(vcf_29_nowarn, grepl("3_prime_UTR_variant", vcf_29_nowarn$V8))),
#                                   nrow(subset(vcf_29_nowarn, grepl("upstream_gene_variant", vcf_29_nowarn$V8))),
#                                   nrow(subset(vcf_29_nowarn, grepl("intergenic_region", vcf_29_nowarn$V8)))),
#                   impact = c("NA", "NA", "NA", "NA", "NA", "modify","low", "modify", "low", "modify", "high", "moderate", "high", "high","high",    
#                              "low","high", "high", "high","low", "low", "modify","modify", "low"))
# 
# write.csv(sum, file = "output/4_load/snpeff/n_mutations_per_type", quote=F, row.names = F)

sum <- read.csv(file = "output/4_load/snpeff/n_mutations_per_type")

n_mutations_pertype <- subset(as.data.frame(sum), impact != "NA")
n_mutations_pertype$impact <- gsub("low", "Low", n_mutations_pertype$impact)
n_mutations_pertype$impact <- gsub("moderate", "Moderate", n_mutations_pertype$impact)
n_mutations_pertype$impact <- gsub("high", "High", n_mutations_pertype$impact)
n_mutations_pertype$impact <- gsub("modify", "Modifier", n_mutations_pertype$impact)
n_mutations_pertype$impact <- factor(n_mutations_pertype$impact, 
                                     levels = c("High", "Moderate", "Low", "Modifier"))

n_mutations_pertype$abb <- c("Downstream gene variant",
                             "5' UTR premature start codon",
                             "5' UTR variant",
                             "Initiator codon variant",
                             "Intron variant",
                             "Loss of Function",
                             "Missense variant",
                             "Nonsense mediated decay",
                             "Splice acceptor variant",
                             "Splice donor variant",
                             "Splice region variant",
                             "Start codon lost",
                             "Stop codon gained",
                             "Stop codon lost",
                             "Stop codon retained",
                             "Synonymous variant",
                             "3' UTR varaint",
                             "Upstream gene variant",
                             "Intergenic region")


ggplot(n_mutations_pertype, aes(x = reorder(abb, desc(n_mutations)), 
                                y = n_mutations)) + 
  geom_col(aes(fill = impact), col = "black") + 
  scale_y_log10(labels = c(expression(paste(~10^1)), expression(paste(~10^3)), expression(paste(~10^5))), 
                breaks = c(10, 1000, 100000))+
  labs(y = expression('Number of SNPs (log'[10]*')'), fill = "Impact Class", x = "Mutation type")+
  scale_fill_manual(values = alpha(c(clr_high, "#8EA4CC","#703D57",  "#FFCD70"), 0.7))+
  scale_color_manual(values = c(clr_high, "#8EA4CC","#703D57",  "#FFCD70"))+
  coord_flip() +geom_text(aes(label = prettyNum(n_mutations, big.mark = ",", scientific=F)), 
                          hjust = 1.5, size = 6)+
  guides(col="none") +
  theme(legend.position = c(0.7,0.8)) -> fig_countsnpef_cat

fig_countsnpef_cat

## simple
snpeff_cat <- sum[c(2:5),]
snpeff_cat$type_clean <- c("High", "Moderate", "Low", "Modifier")
snpeff_cat$clr <- c("s", "ns", "ns", "ns")

snpeff_cat$type_clean <- factor(snpeff_cat$type_clean, levels = c("Modifier", "Low", "Moderate", "High"))

ggplot(snpeff_cat, aes(x = type_clean, y = n_mutations)) + geom_col(aes(fill = clr), col = "black") + 
  labs(x = "SnpEff impact category", y= expression('Number of SNPs (log'[10]*')')) +
  scale_y_log10(labels = c(expression(paste(~10^1)), expression(paste(~10^4)), 
                           expression(paste(~10^7))),
                breaks=c(10, 10000, 10000000))+
  geom_text(aes(label = prettyNum(n_mutations, big.mark=","), y = n_mutations), 
            hjust=1.5, size = 6) +
  scale_fill_manual(values = c(alpha(clr_grey, 0.7), alpha(clr_high, 0.7)))+
  coord_flip() + theme(legend.position = "none") -> simple_snpeff_cat

#### c: gerp score distribution #####
# gerp <- data.frame()
# for (i in c(1:3,5:15,17:28,30)){
#   scaf <- read.table(paste0("output/4_load/gerp/overlap/gerp_overlapSNP_scaf_",i, ".tsv.gz"))
#   scaf <- scaf[,c(1:10)]
#   gerp <- rbind(gerp, scaf)
# }
#  
# save(gerp, file = "output/3_annotated_genome/gerp_scores_all.RData")

load(file = "output/3_annotated_genome/gerp_scores_all.RData")
ggplot(gerp, aes(x = V5)) + geom_histogram(fill = clr_grey, col = "black") + 
  scale_y_continuous(breaks = c(0, 200000, 400000, 600000), labels = c("0","2", "4", "6"))+
  labs(x = "GERP score", y = expression("Count (x"~10^5~")")) -> gerp_score_dist

gerp_score_dist

gerp_per_cat <- data.frame(cat = c("< 0", "0-1", "1-2", "2-3", "3-4", "≥ 4"),
                           n = c(nrow(subset(gerp, V5 < 0)),
                                   nrow(subset(gerp, V5 >= 0 & V5 < 1)),
                                   nrow(subset(gerp, V5 >= 1 & V5 < 2)),
                                   nrow(subset(gerp, V5 >= 2 & V5 < 3)),
                                   nrow(subset(gerp, V5 >= 3 & V5 < 4)),
                                 nrow(subset(gerp, V5 >= 4))),
                           clr = c("ns", "ns", "ns", "ns", "ns", "s"))

gerp_per_cat$cat <- factor(gerp_per_cat$cat, levels = c("< 0", "0-1", "1-2", "2-3", "3-4", "≥ 4"))

ggplot(gerp_per_cat, aes(x = cat, y = n)) + geom_col(aes(fill = clr), col = "black") + 
  labs(x = "GERP score category", y= expression('Number of SNPs (log'[10]*')')) +
  scale_y_log10(limits=c(1,10000000), labels = c(expression(10^1), expression(10^3), expression(10^6)),
                breaks=c(1, 1000,1000000)) +
  geom_text(aes(label = prettyNum(n, big.mark=","), y = n), 
            hjust=1.5, size = 6) +
  scale_fill_manual(values = c(alpha(clr_grey, 0.7), alpha(clr_high, 0.7)))+
  coord_flip() + theme(legend.position = "none") -> simple_gerp_cat

#### d: gerp load distribution #####
load(file = "output/loads.RData")

loads <- left_join(loads, pheno_wide[,c("id", "lifespan")], by = "id")
loads <- loads %>% mutate(age = as.factor(case_when(
  is.na(lifespan) ~ "Chick",
  lifespan == 1 ~ "Yearling",
  lifespan > 1 ~ "Adult"
)))

loads$age <- factor(loads$age, levels = c("Chick", "Yearling", "Adult"))

ggplot(subset(loads, loadtype == "gerp"), aes(x = total_load, fill = age)) + geom_density() +
  labs(x = "Total GERP load", y = "Density", fill = "Age") +
  scale_fill_manual(values = c(alpha(clrs_hunting[2], 0.7),
                               alpha(clrs_hunting[3], 0.7), 
                               alpha(clrs_hunting[4], 0.7))) -> gerp_dist

gerp_dist

#### e: snpeff load distribution #####

ggplot(subset(loads, loadtype == "high"), aes(x = total_load, fill = age)) + geom_density() +
  labs(x = "Total SnpEff load", y = "Density", fill = "Age") +
  scale_fill_manual(values = c(alpha(clrs_hunting[2], 0.7),
                               alpha(clrs_hunting[3], 0.7), 
                               alpha(clrs_hunting[4], 0.7))) -> high_dist

high_dist

#### combine ####
plot_grid(gerp_dist, high_dist,
          ncol = 1, align = "hv", axis = "lb",
          labels = c("D", "E"), label_fontface = "plain", label_size = 22) -> right

plot_grid(simple_gerp_cat, simple_snpeff_cat,  
          ncol = 2, align = "hv", axis = "lb",
          labels = c("A", "B"), label_fontface = "plain", label_size = 22) -> left_1
          
plot_grid(left_1, fig_countsnpef_cat, rel_heights = c(0.7, 1),
          ncol = 1, 
          labels = c("", "C"), label_fontface = "plain", label_size = 22) -> left_2

plot_grid(left_2, right, rel_widths = c(1, 0.8),
          ncol = 2, #align = "hv", axis = "lb",
          labels = c("", ""), label_fontface = "plain", label_size = 22) -> sup

sup
ggsave(sup, file = "plots/sup_fig_dist_froh_mutations.png", width=18,height=14)

