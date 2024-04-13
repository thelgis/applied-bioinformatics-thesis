library(edgeR)
library(ggplot2)

setwd("/Users/thelgis/Documents/source/applied-bioinformatics-thesis")

source("R/differential_expression_libs.R")

run_edgeR('RA_RNA_SEQ_SYNOVIAL_MEMBRANE.csv')
run_edgeR('SLE_RNA_SEQ_WHOLE_BLOOD.csv', logFC_threshold=0.5)  # logFC_threshold=1 returns only one gene
