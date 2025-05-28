#!/usr/bin/env Rscript

packages.list <- c("tidyverse","vcfR")
for(pkg in packages.list){ eval(bquote(library(.(pkg)))) }

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)

# Function to display usage
print_usage <- function() {
  cat("Usage: Rscript Process_multi-VCF_to_CSV.R <input_file>\n")
  cat("  <input_file>   Path to input file (required)\n")
}

# Error checking
if (length(args) < 1) {
  cat("Error: Not enough arguments provided.\n\n")
  print_usage()
  quit(status = 1)
}

# assign argument to input file
multivcf.file <- args[1]

# Validate input file exists
if (!file.exists(multivcf.file)) {
  cat("Error: Input file does not exist.\n")
  quit(status = 1)
}

#multivcf.file <- "/Users/mb29/Papers/AmpliSeq_Treponema_paper_2024/Revision_04-2025/Gono_work/clean.SNPs.vcf"

#######################################
# Function to process vcf file

# read_and_process_vcf <- function(myvcf.path){
#   sim.dataset.vcf <- read.vcfR(myvcf.path, verbose = FALSE)
#   sim.dataset.vcf.fix <- getFIX(sim.dataset.vcf)
#   sim.dataset.vcf.fix <- data.frame(sim.dataset.vcf.fix[,c(2,4,5)], stringsAsFactors = F)
#   sim.dataset.vcf.fix$Key <- 1:nrow(sim.dataset.vcf.fix)
#   
#   sim.dataset.vcf.gt <- extract_gt_tidy(sim.dataset.vcf)
#   sim.dataset.vcf.gt.f <- plyr::join(sim.dataset.vcf.gt, sim.dataset.vcf.fix[,c("Key","POS","REF")], by="Key", type="left")
#   sim.dataset.vcf.gt.f$POS <- as.numeric(sim.dataset.vcf.gt.f$POS)
#   sim.dataset.vcf.gt.f$gt_GT <- as.numeric(sim.dataset.vcf.gt.f$gt_GT)
#   return(sim.dataset.vcf.gt.f)
# }

# Streamlined version of code
read_and_process_vcf <- function(myvcf.path){
  sim.dataset.vcf.file <- read.vcfR(myvcf.path, verbose = FALSE, limit=1e10)
  sim.dataset.vcf <- vcfR2tidy(sim.dataset.vcf.file, single_frame=T, toss_INFO_column = TRUE)$dat
  sim.dataset.vcf <- sim.dataset.vcf %>% select(Indiv, POS, REF, ALT, gt_GT_alleles, gt_GT)
  sim.dataset.vcf$POS <- as.numeric(sim.dataset.vcf$POS)
  sim.dataset.vcf$gt_GT <- as.numeric(sim.dataset.vcf$gt_GT)
  return(sim.dataset.vcf)
}

# code to spread data
spread_processed_data <- function(processed_data){
  processed.spread.temp <- tidyr::spread(processed_data[,c("POS","Indiv","gt_GT")], POS, gt_GT)
  return(processed.spread.temp)
}


  
#######################################

multivcf.f <- read_and_process_vcf(multivcf.file)

multivcf.f.spread <- spread_processed_data(multivcf.f)

write.csv(multivcf.f, file=paste0(multivcf.file,"_processed.tsv"), row.names = F, quote=F)
write.csv(multivcf.f.spread, file=paste0(multivcf.file,"_processed_spread.tsv"), row.names = F, quote=F)


