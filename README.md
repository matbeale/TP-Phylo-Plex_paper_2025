# TP-Phylo-Plex_paper_2025

Raw & intermediate data, code, and analysis underlying design of a phylogenetically informed multiplex amplicon sequencing assay "TP-Phylo-Plex" for Treponema pallidum. Includes a complete RNotebook including primary analysis and figure generation. 
 
## Code testing
To use this code, first clone the repo
```
git clone https://github.com/matbeale/TP-Phylo-Plex_paper_2025
```

Navigate to the repository directory and unzip the largest sim files:
``` 
gunzip data/simulated_mutations_250-sims/*
```

Then open the RNotebook (`Treponema_Phylo-Seq_Discrimintory_amplicon_analysis_02-2025b.Rmd`) using Rstudio.
\
<br />

The Rnotebook contains internal paths for the repo, and should find all data. In some cases, lengthy steps are set to ```eval=F``` or commented out, with the intermediate results files available to continue running the code (although those steps can still be run).


<br />
Ensure all R package dependencies are installed and available:
\

| R package  | Version | URL                                    |
|------------|---------|----------------------------------------|
| ggplot2    | 3.5.1   | CRAN                                   |
| dplyr      | 1.1.4   | CRAN                                   |
| tidyverse  | 2.0.0   | CRAN                                   |
| treeio     | 1.28.0  | https://github.com/YuLab-SMU/treeio    |
| ggtree     | 3.12.0  | https://github.com/YuLab-SMU/ggtree    |
| phytools   | 2.1-1   | CRAN                                   |
| readxl     | 1.4.3   | CRAN                                   |
| ggnewscale | 0.4.10  | CRAN                                   |
| ape        | 5.8     | CRAN                                   |
| adegenet   | 2.1.10  | CRAN                                   |
| pegas      | 1.3     | CRAN                                   |
| mmod       | 1.3.3   | CRAN                                   |
| vcfR       | 1.15.0  | CRAN                                   |
| hierfstat  | 0.5-11  | CRAN                                   |
| network    | 1.18.2  | CRAN                                   |
| ggnetwork  | 0.5.13  | CRAN                                   |
| intergraph | 2.0-4   | CRAN                                   |
| igraph     | 2.0.3   | CRAN                                   |
| seqinr     | 4.2-36  | CRAN                                   |
| IRanges    | 2.38.0  | CRAN                                   |
| ggbeeswarm | 0.7.2   | CRAN                                   |
| primerTree | 1.0.6   | CRAN                                   |
| taxize     | 0.9.100 | CRAN                                   |
| cowplot    | 1.1.3   | CRAN                                   |
| ggstance   | 0.3.7   | CRAN                                   |
| ggforce    | 0.4.2   | CRAN                                   |
| vegan      | 2.6-6.1 | CRAN                                   |
| Polychrome | 1.5.1   | CRAN                                   |
| ggrastr    | 1.0.2   | CRAN                                   |
| treemapify | 2.5.6   | CRAN                                   |
| treespace  | 1.1.4.3 | CRAN                                   |
| reshape2   | 1.4.4   | CRAN                                   |
| doMC       | 1.3.8   | CRAN                                   |
| pairsnp    | 0.1.0   | https://github.com/gtonkinhill/pairsnp |
| rpinecone  | 0.1.0   | https://github.com/alexwailan/rpinecone|


<br />


