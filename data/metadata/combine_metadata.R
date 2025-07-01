#combining data from Carl to chick IDs
pacman::p_load(readxl, dplyr, lubridate)

carl <- read_excel("data/metadata/Chick samples to send.xlsx")
id <- read_excel("data/metadata/Chick_sample_list_extractions.xlsx")

#combine
chick <- left_join(id[,c(1,2,4,5,7)], carl[,c(1,4:9)], by = c("Location in box"))
names(chick) <- c("box_location", "year", "day_collect", "extract", "id", 
                  "site", "mother", "father", "no_egg", "date_hatch_incorrect", 
                  "no_chicks")

chick$date_collect <- as.Date(paste0(substr(chick$day_collect, 2,3), "-0", 
                                     substr(chick$day_collect, 1,1), "-",
                                     substr(chick$year, 1, 4)), format="%d-%m-%y")

chick$date_hatch <- as.Date(paste0(day(chick$date_hatch_incorrect), "-", 
                                   month(chick$date_hatch_incorrect), "-",
                                   substr(chick$year, 1, 4)), format="%d-%m-%y")

chick$age <- chick$date_collect - chick$date_hatch

#check repeated moms
length(unique(chick$mother)) == nrow(chick)

#replace 'UNK' fathers with NA
chick$father <- gsub("UNK", NA, chick$father)
summary(as.factor(chick$father))
#some fathers sired multiple chicks
summary(as.factor(chick$site))

#combine with adult data
load("metadata/phenotypes_wide_extra.RData")
#combine
meta_chick <- chick %>% select(c(id, site, father, year)) 
meta_adult <- pheno %>% filter(!is.na(LMS_min)) %>% select(c(id, site, born))
meta_adult <- meta_adult %>% mutate(father = NA, .before = born) %>%
  rename(year = born)

meta <- rbind(meta_chick, meta_adult)

#remove C47 and C48 because they didn't reach quality threshold
meta <- subset(meta, id != "C47" & id != "C48")
save(meta, file = "metadata/metadata_adult_chick.RData")
