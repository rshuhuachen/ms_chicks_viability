pacman::p_load(tidyverse, data.table)

list_files <- fread("data/genomic/genomes/list_bam_all.txt", header = F)
list_files <- c(list_files$V1)

n_reads <- data.frame()
for (i in 1:length(list_files)){
  n_reads_file <- system(paste0("samtools view -c ", list_files[i]), intern = TRUE)
  n_reads_i <- data.frame(file = list_files[i], n_reads = n_reads_file)
  n_reads <- rbind(n_reads, n_reads_i)
} 

