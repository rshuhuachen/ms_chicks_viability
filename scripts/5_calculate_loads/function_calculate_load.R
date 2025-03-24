calculate_load_snpeff <- function(vcf, output_vcf, loadtype){
  # only select 29 largest autosomal scaffolds
  scaf <- fread("data/metadata/30_largest.scafs.tsv")
  scaf$scaf <- gsub(":", ";", scaf$scaf)
  scaf$scaf <- gsub("\\.", "=", scaf$scaf)
  scaf <- subset(scaf, scaf_no != 4)
  
  vcf <- subset(vcf, CHROM %in% scaf$scaf)
  
  # exclude warning messages
  
  vcf <- subset(vcf, !grepl("WARNING", INFO))
  
  # only get GT info, PL and DP are filtered by already anyway 
  gt <- c(10:ncol(vcf))
  select_n3 <- function(x){x = substr(x,1,3)}
  vcf[gt] <- lapply(vcf[gt], select_n3)
  
  # calculate load
  load <- list()
  # loop over ids
  for( id in 10:(ncol(vcf))){
    # subset per id
    subset_id <- vcf[,c(1:9, id)]
    
    # filter for snps in het and hom
    het_data <- subset(subset_id, subset_id[[10]] == "1/0" | subset_id[[10]] == "0/1")
    hom_data <- subset(subset_id, subset_id[[10]] == "1/1")
    
    # count amount of snps in het and hom
    het_load_sum <- nrow(het_data)
    hom_load_sum <- nrow(hom_data)
    
    # count no of snps successfully genotyped
    n_genotyped <- nrow(subset_id) - nrow(subset(subset_id, subset_id[[10]] == "./."))
    n_total <- nrow(subset_id)
    
    # collect data in df
    df <- data.frame(id = colnames(vcf[id]),
                     n_total = n_total,
                     n_genotyped = n_genotyped,
                     het_load = het_load_sum / n_genotyped,
                     hom_load = hom_load_sum / n_genotyped,
                     total_load = (het_load_sum*0.5 + hom_load_sum) / n_genotyped,
                     loadtype = loadtype)
    load[[id]] <- df
  }
  # convert list to df
  load <- do.call(rbind.data.frame, load)
  load <- subset(load, id != "C09") # remove the individual that is probably sampled as an adult
  if(output_vcf == TRUE){
    out <- list(load = load, vcf = vcf)
    return(out)
  }
  
  if(output_vcf==FALSE){
    return(load)}
}

calculate_load_gerp <- function(vcf, output_vcf, loadtype){
  
  # only get GT info, PL and DP are filtered by already anyway 
  gt <- c(11:ncol(vcf))
  select_n3 <- function(x){x = substr(x,1,3)}
  vcf[gt] <- lapply(vcf[gt], select_n3)
  
  # calculate load
  load <- list()
  # loop over ids
  for( id in 11:(ncol(vcf))){
    # subset per id
    subset_id <- vcf[,c(1:10, id)]
    
    # filter for snps in het and hom
    het_data <- subset(subset_id, subset_id[[11]] == "1/0" | subset_id[[11]] == "0/1")
    hom_data <- subset(subset_id, subset_id[[11]] == "1/1")
    
    # count amount of snps in het and hom
    het_load_sum <- nrow(het_data)
    hom_load_sum <- nrow(hom_data)
    
    # count no of snps successfully genotyped
    n_genotyped <- nrow(subset_id) - nrow(subset(subset_id, subset_id[[11]] == "./."))
    n_total <- nrow(subset_id)
    
    # collect data in df
    df <- data.frame(id = colnames(vcf[id]),
                     n_total = n_total,
                     n_genotyped = n_genotyped,
                     het_load = het_load_sum / n_genotyped,
                     hom_load = hom_load_sum / n_genotyped,
                     total_load = (het_load_sum*0.5 + hom_load_sum) / n_genotyped,
                     loadtype = loadtype)
    load[[id]] <- df
  }
  # convert list to df
  load <- do.call(rbind.data.frame, load)
  load <- subset(load, id != "C09")
  if(output_vcf == TRUE){
    out <- list(load = load, vcf = vcf)
    return(out)
  }
  
  if(output_vcf==FALSE){
  return(load)}
}


