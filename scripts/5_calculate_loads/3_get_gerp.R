library(dplyr); library(data.table)

gerp <- data.frame()
for (i in c(1:3,5:15,17:28,30)){
  scaf <- read.table(paste0("output/4_load/gerp/overlap/gerp_overlapSNP_scaf_",i, ".tsv.gz"))
  scaf <- subset(scaf, V5 >= 4) #don't forget the =!
  gerp <- rbind(gerp, scaf)
}

#### rename columns #####
### load headers
names <- fread("output/2_inbreeding/chicks_adults_samples.txt", header = F)
names(names) <- "file_id"

### load id info
#add real id
ids <- read.csv("data/metadata/file_list_all_bgi_clean.csv")
ids$file_id <- gsub(".sorted.bam", "", ids$file)
names <- left_join(names, ids[,c("file_id", "id")], by = c("file_id"))
names <- names %>% mutate(id = case_when(
  is.na(id) ~ file_id,
  TRUE ~ id
))

names(gerp) <- c(c("chr", "start", "pos", "neutral_rate_n", "rs_score", "ancestral", "derived", "qual", "info","format"), names$id)

save(gerp, file = "output/4_load/gerp/gerps_all.RData")
