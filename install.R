pkgs = c("knitr","rmarkdown","vegan","ggplot2","ggpubr","reshape2","cowplot","superheat","plyr","dplyr")
ncores = parallel::detectCores()
install.packages(pkgs, Ncpus = ncores)

source("http://www.bioconductor.org/biocLite.R")
biocLite("DESeq2")
biocLite("phyloseq")

install.packages("devtools")
devtools::install_github("benjjneb/decontam")
