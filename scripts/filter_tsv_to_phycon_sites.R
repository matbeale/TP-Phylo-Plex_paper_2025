#!/usr/bin/env Rscript

packages.list <- c("dplyr","tidyr")
for(pkg in packages.list){ eval(bquote(library(.(pkg)))) }

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)

# Function to display usage
print_usage <- function() {
  cat("Usage: Rscript filter_tsv_to_phycon_sites.R <input_file> <position_list>\n")
  cat("  <SNP_table>   Path to input SNP_table file (required)\n")
  cat("  <position_list>   Path to input position list (required)\n")
}

# Error checking
if (length(args) < 2) {
  cat("Error: Not enough arguments provided.\n\n")
  print_usage()
  quit(status = 1)
}

# assign argument to input file
SNP_table.file <- args[1]
position_list.file <- args[2]

# Validate SNP_table file exists
if (!file.exists(SNP_table.file)) {
  cat("Error: SNP_table.file does not exist.\n")
  quit(status = 1)
}

# Validate position_list file exists
if (!file.exists(position_list.file)) {
  cat("Error: position_list.file does not exist.\n")
  quit(status = 1)
}

#Taouk_gono.processed.vcf.file <- "/Users/mb29/Papers/AmpliSeq_Treponema_paper_2024/Github/TP-Phylo-Plex_paper_2025/Gono_data/clean.SNPs.vcf.gz_processed_spread.tsv.gz"


######################################################################
# Read in SNPs table file
SNP_table.processed.vcf <- read.csv(SNP_table.file, check.names = F, row.names = NULL, header=T, col.names=c("Indiv","POS","REF","ALT","ALT2","gt_GT_alleles","gt_GT")) %>% mutate(POS.n = as.character(POS))

print("check SNP table to evaluate")
head(SNP_table.processed.vcf)
class(SNP_table.processed.vcf$POS)
class(SNP_table.processed.vcf$POS.n)
names(SNP_table.processed.vcf)
head(SNP_table.processed.vcf$POS)


# pivot longer
#SNP_table.processed.vcf.melt <- SNP_table.processed.vcf %>%
#  pivot_longer(-Indiv, names_to="POS", values_to="Allele")
  
# Read in positions list file
position_list <- read.csv(position_list.file, check.names = F, header=F, col.names = c("POS")) %>% mutate(POS.n = as.character(POS))

print("check positions to extract")
head(position_list)
class(position_list$POS)
class(position_list$POS.n)
names(position_list)

cat("Rows in SNP table:", nrow(SNP_table.processed.vcf), "\n")
cat("Rows in position list:", nrow(position_list), "\n")


# filter to sites
#SNP_table.processed.vcf.melt_in.phycons <- SNP_table.processed.vcf.melt %>% 
#  filter(POS %in% position_list$V1)

SNP_table.processed.vcf.melt_in.phycons <- SNP_table.processed.vcf %>%
  filter(POS.n %in% (position_list$POS.n))


cat("Rows after filtering:", nrow(SNP_table.processed.vcf.melt_in.phycons), "\n")

######################################################################
# Write out to file
write.csv(SNP_table.processed.vcf.melt_in.phycons, file=paste0(SNP_table.file,"_filtered_phycons.tsv"), row.names = F, quote=F)


