library(edgeR)
library(ggplot2)

setwd("/Users/thelgis/Documents/source/applied-bioinformatics-thesis")

source("R/differential_expression_libs.R")

run_edgeR('RA_RNA_SEQ_SYNOVIAL_MEMBRANE.csv')
# Positive Results:
# - Good separation of Healthy against disease before and after DE
# - Explained variation before DE was 28%/12% and moved to 51%/5%
# - 148 differentially expressed genes
# - Depends on two datasets GSE89408(good initial separation) and GSE90081(bad initial separation - consider excluding)

run_edgeR('SLE_RNA_SEQ_WHOLE_BLOOD.csv', logFC_threshold=0.5)  # logFC_threshold=1 returns only one gene
# Negative Results:
# - Bad separation of Healthy against disease before and after DE
# - Only 7 differentially expressed genes and even that requires to lower the logFC threshold

run_edgeR('SSc_RNA_SEQ_WHOLE_BLOOD.csv')
# Negative Results:
# - No differentially expressed genes (also hinted by initial PCA)

run_edgeR('SSc_RNA_SEQ_PERIPHERAL_BLOOD.csv')
# Negative Results:
# - No differentially expressed genes (also hinted by initial PCA)
