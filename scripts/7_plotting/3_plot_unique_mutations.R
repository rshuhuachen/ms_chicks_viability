## load packages ####
pacman::p_load(dplyr, data.table, ggplot2)

source("scripts/theme_ggplot.R")

## load summary data ####
load(file= "output/3_annotated_genome/summary_unique_mutations_40ids.RData")

ggplot(summary, aes(x = rxy)) + geom_histogram(fill = clr_grey, col = "black") + 
  labs(x = expression(R[XY]), y = "Count") + 
  facet_grid(~method, scales="free")+
  geom_vline(xintercept = 1, col = "darkred", linetype = "dotted") -> plot_rxy

plot_rxy

ggplot(summary, aes(x = n_unique_in_chick/n_unique_in_ad)) + geom_histogram(fill = clr_grey, col = "black") + 
  labs(x = "Ratio unique mutations in chicks:(sub)adults", y = "Count") + 
  facet_grid(~method, scales="free")+
  geom_vline(xintercept = 1, col = "darkred", linetype = "dotted") -> plot_ratio

plot_ratio

ggsave(plot_rxy, file = "plots/unique_mutations_rxy.png", width = 8, height = 6)
ggsave(plot_ratio, file = "plots/unique_mutations_ratio.png", width = 8, height = 6)
