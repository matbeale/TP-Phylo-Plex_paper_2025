#!/usr/bin/env Rscript

packages.list <- c("tidyverse","hierfstat", "foreach", "doParallel")
for(pkg in packages.list){ eval(bquote(library(.(pkg)))) }

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)

# Function to display usage
print_usage <- function() {
  cat("Usage: Rscript Calculate_FST_from_populations.R  <variants.csv> <populations.csv> <FST.threshold>\n")
  cat("  <variants.csv>  Path to variants file (already generated from vcf - spread version) (required)\n")
  cat("  <populations.csv>   Path to populations file (2 columns: sample,population) (required)\n")
  cat("  <fst.threshold>   Float (required)\n")
}

# Error checking
if (length(args) < 3) {
  cat("Error: Not enough arguments provided.\n\n")
  print_usage()
  quit(status = 1)
}

# assign argument to input file
processed.vcf.file <- args[1]
populations.file <- args[2]
fst.threshold <- as.numeric(args[3])



################################################################################################ 
# Functions
################################################################################################

#### Function to process into population matrix
convert_populations_to_matrix <- function(input.table){
  # relabel csv headers
  colnames(input.table) <- c("sampleid", "population")
  # Get list of unique populations
  input.pops.list <- unique(input.table$population)
  # loop to convert into binary population matrix
  pops.matrix <- input.table
  for (current in input.pops.list){
    pops.matrix$sublin.class <- ifelse(pops.matrix$population==current,0,1)
    colnames(pops.matrix)[ncol(pops.matrix)] <- current
  }
  return(pops.matrix)
}
#convert_populations_to_matrix(input.pops)
################################################################################################
#
#
# Function to calculate site-by-site fst for one set of populations
calculate_fst_for_pop <- function(poplabel, my.locus.df, mypop.data, FST.cuttoff){
  # Extract relevant columns into separate dataframes
  # Dataframe containing multiallelic codes as integers
  #my.locus.df.gt.f.spread.temp <- tidyr::spread(my.locus.df[,c("POS","Indiv","gt_GT")], POS, gt_GT) 
  my.locus.df.gt.f.spread.temp <- my.locus.df
  # Fix discrepancy between names (vcfR trims sample names)
  mypop.data.mod <- mypop.data %>% mutate(sampleid = gsub("\\_.+$","", sampleid))
  # combine allele information with population definitions
  fst.join.meta <- left_join(data.frame(sampleid=my.locus.df.gt.f.spread.temp$Indiv), mypop.data.mod, by="sampleid")
  fst.join.pops <- cbind(fst.join.meta[,2], my.locus.df.gt.f.spread.temp)
  # calculate site by site fst for pops using Weir and Cockrham estimates of Fstatistics (hierfstat package)
  fst.run <- hierfstat::wc(fst.join.pops[,c(1,3:ncol(fst.join.pops))], diploid=F)
  # fst.out contains the fst values for each site evaluated (and specifies whether they meet the defined cuttoff)
  fst.out <- data.frame(POS=as.numeric(gsub("X","",row.names(fst.run$per.loc))), FST=fst.run$per.loc$FST, stringsAsFactors = F)
  fst.out$population <- poplabel
  fst.out$sig <- ifelse(fst.out$FST>FST.cuttoff, "sig","other")
  fst.out.list <- fst.out
  return(fst.out.list)
}
###
# Usage Example: 
# FST.species.TEN.1.1 <- calculate_fst_for_pop("TEN",variants.df, pop.matrix, 0.99)
################################################################################################
#
#
## function to loop through all populations and assess variants
assess_variants_loop <- function(variants.df, pop.matrix, fst.threshold){
  FST.combined.output <- NULL
  for (current in sort(unique(pop.matrix$population))){
    print(paste0("Testing variants against population:",current))
    FST.sublin.current <- calculate_fst_for_pop(current,variants.df, pop.matrix[,c("sampleid",current)],fst.threshold)
    FST.combined.output <- rbind(FST.combined.output, FST.sublin.current)
  }
  FST.combined.output <- FST.combined.output[!is.na(FST.combined.output$sig),]
  return(FST.combined.output)
}
###
# Usage Example: 
# assess_variants_loop(variants.df, pop.matrix, 0.99)
################################################################################################



#
#
## function to loop through all populations and assess variants

# Parallel version of assess_variants_loop
assess_variants_parallel <- function(variants.df, pop.matrix, fst.threshold, n.cores = 4){
  # Setup parallel backend
  cl <- makeCluster(n.cores)
  registerDoParallel(cl)
  clusterExport(cl, varlist = c("calculate_fst_for_pop"))
  
  populations <- sort(unique(pop.matrix$population))
  
  # Run in parallel
  FST.list <- foreach(current = populations, .combine = rbind, .packages = c("tidyverse","hierfstat")) %dopar% {
    print(paste0("Testing variants against population: ", current))
    FST.sublin.current <- calculate_fst_for_pop(
      current,
      variants.df,
      pop.matrix[, c("sampleid", current)],
      fst.threshold
    )
    FST.sublin.current
  }
  
  # Stop the cluster
  stopCluster(cl)
  
  # Filter NAs
  FST.list <- FST.list[!is.na(FST.list$sig), ]
  return(FST.list)
}


###
# Usage Example: 
# assess_variants_loop(variants.df, pop.matrix, 0.99)
################################################################################################





###################




# Main:
print("Reading in populations file")
populations.input <- read.csv(populations.file)

print("Generating population matrix")
populations.matrix <- convert_populations_to_matrix(populations.input)

print("Reading in processed variant data")
processed.vcf <- read.csv(processed.vcf.file, check.names = F)


print("Assessing variant sites for FST against each population")
# Example:
# FST.sites.by.population.test <- calculate_fst_for_pop("TPA-1", processed.vcf, populations.matrix[,c("sampleid","TPA-1")], 0.99)

#FST.sites.by.population <- assess_variants_loop(processed.vcf, populations.matrix, fst.threshold)
FST.sites.by.population <- assess_variants_parallel(processed.vcf, populations.matrix, fst.threshold, n.cores = 12)


output.filename <- paste0(processed.vcf.file,".FST.scoring.csv")

print(paste0("Writing analysis out to file: ", output.filename))
write.csv(FST.sites.by.population, file=output.filename, row.names = F)
