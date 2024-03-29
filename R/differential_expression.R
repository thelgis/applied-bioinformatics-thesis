library(edgeR)

setwd("/Users/thelgis/Documents/source/applied-bioinformatics-thesis")

RA_RNA_SEQ_SYNOVIAL_MEMBRANE_with_metadata <- read.table(
  "data/r_pre_processed_datasets/RA_RNA_SEQ_SYNOVIAL_MEMBRANE.csv",
  header = TRUE,
  sep = ",", 
  row.names = 1
)

# Determine groups vector 
groups_df <- RA_RNA_SEQ_SYNOVIAL_MEMBRANE_with_metadata[c("Condition"),]
groups <- as.character(groups_df[1,])
rm(groups_df)

# Remove metadata columns 
RA_RNA_SEQ_SYNOVIAL_MEMBRANE <- RA_RNA_SEQ_SYNOVIAL_MEMBRANE_with_metadata[1:(nrow(RA_RNA_SEQ_SYNOVIAL_MEMBRANE_with_metadata) - 11), ]
rm(RA_RNA_SEQ_SYNOVIAL_MEMBRANE_with_metadata)

# Convert to numeric matrix 
RA_RNA_SEQ_SYNOVIAL_MEMBRANE_matrix <- as.matrix(sapply(RA_RNA_SEQ_SYNOVIAL_MEMBRANE, as.numeric))  
row.names(RA_RNA_SEQ_SYNOVIAL_MEMBRANE_matrix) <- row.names(RA_RNA_SEQ_SYNOVIAL_MEMBRANE)
