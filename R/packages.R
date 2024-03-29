install.packages("tidyverse")

if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install("edgeR")
BiocManager::install("DESeq2")
