#!/usr/bin/env python
import sys
import io
import subprocess
import os
import pandas as pd
import numpy as np
import re
import argparse
import csv

# Mat Beale, Wellcome Sanger Institute, July 2022

description="Takes a bam file and a list of regions, and uses samtools depth to determine the mean depth per region, the normalised depth (useful for CNV detection), and the number of missing (zero coverage) sites per gene. Note that this script currently assumes a single reference chromosome \n"

parser = argparse.ArgumentParser(description=description)
parser.add_argument("-b", type=str, dest="bam_file",help="Specify input bam file")
parser.add_argument("-w", dest="wgs", action="store_true", help="Calculate stats for Whole Genome" )
parser.add_argument("-r", type=str, dest="region", help="Comma separated region positions (<start>,<end>,<name>)" )
parser.add_argument("-c", type=str, dest="csv", help="CSV file containing comma separated region positions per line (<start>,<end>,<name>)" )
parser.set_defaults(feature=False)

args = parser.parse_args()

bam_file = args.bam_file



# Get sequence name 
mysamplename = re.sub("^.+/","",re.sub(".bam","", bam_file))


# Specify location or path of samtools version
samtools_exe = 'samtools'

# Use samtools to generate depth data and read from stdin:
run_samtools = subprocess.Popen([samtools_exe,"depth","-aa",bam_file],stdout=subprocess.PIPE)
csv = io.StringIO()
for line in run_samtools.stdout:
    csv.write(line.decode().strip('"\n') + '\n')
csv.seek(0)
mydepth = pd.read_csv(csv, sep='\t',index_col=0, header=None,names=['Ref','Pos','Cov'])
csv.close()

# Anything that isn't covered (is na) make 0
mydepth[pd.isnull(mydepth['Cov'])] = 0


# Calculate genomewide median and mean coverage (median is probably more reliable here)
full_median = np.median(mydepth['Cov'])
full_mean = np.mean(mydepth['Cov'])
wgs_length = len(list(mydepth['Pos']))

# Now define a function to extract a region and calculate some summary coverage stats
def get_covwindow(all_cov, mystart, myend, region_name):
    genedepth = mydepth.loc[(mydepth['Pos']>=int(mystart)-1) & (mydepth['Pos']<=int(myend)-1)]
    gene_length = len(list(genedepth['Cov']))
    missing_sites = len(list(genedepth['Cov'][genedepth['Cov']==0]))
    gene_median = round(np.median(genedepth['Cov']),1)
    gene_mean = np.mean(genedepth['Cov'])
    mincov = np.min(genedepth['Cov'])
    maxcov = np.max(genedepth['Cov'])
    meancov = round(gene_mean,1)
    normcov = round(gene_mean/full_median,1)
    perc1x = round((len(list(genedepth['Cov'][genedepth['Cov']>=1]))/gene_length)*100,1)
    perc5x = round((len(list(genedepth['Cov'][genedepth['Cov']>=5]))/gene_length)*100,1)
    perc8x = round((len(list(genedepth['Cov'][genedepth['Cov']>=8]))/gene_length)*100,1)
    perc20x = round((len(list(genedepth['Cov'][genedepth['Cov']>=20]))/gene_length)*100,1)
    perc100x = round((len(list(genedepth['Cov'][genedepth['Cov']>=100]))/gene_length)*100,1)
    perc250x = round((len(list(genedepth['Cov'][genedepth['Cov']>=250]))/gene_length)*100,1)

    return[mysamplename, region_name, str(mystart), str(myend), str(gene_length), str(missing_sites), str(gene_median), str(meancov), str(mincov), str(maxcov), str(normcov), str(perc1x), str(perc5x), str(perc8x), str(perc20x), str(perc100x), str(perc250x)]



header = ["Sample", "Region", "StartPos", "EndPos", "Length", "Missing_Sites", "Median_Cov", "Mean_Cov", "Min_Cov", "Max_Cov", "Genome_Norm_Cov", "cov1x.perc", "cov5x.perc", "cov8x.perc", "cov20x.perc", "cov100x.perc", "cov250x.perc"]
dfout = pd.DataFrame(columns = header)


# Run function on a whole genome
if args.wgs :
    region_cov_summary = get_covwindow(mydepth, "1", wgs_length, "WGS")
    dfout = dfout.append(pd.DataFrame([region_cov_summary], columns=header), ignore_index=True)


# Run function on a single region
if args.region :
    myregion = [s.strip() for s in args.region.split(",")]
    region_cov_summary = get_covwindow(mydepth, myregion[0], myregion[1], myregion[2])
    dfout = dfout.append(pd.DataFrame([region_cov_summary], columns=header), ignore_index=True)



# Run function on a csv list of regions
if args.csv : 
    csv_file = pd.read_csv(args.csv, header=None)
    for myline in range(0, len(csv_file.index)):
        region_cov_summary = get_covwindow(mydepth, csv_file.iloc[myline,0], csv_file.iloc[myline,1], csv_file.iloc[myline,2])
        dfout = dfout.append(pd.DataFrame([region_cov_summary], columns=header), ignore_index=True)


#print(dfout)
dfout.to_csv(mysamplename + ".depth.tsv", sep="\t", index=False, header=True)

