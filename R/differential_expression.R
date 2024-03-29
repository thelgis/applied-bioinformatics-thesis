library(edgeR)
library(ggplot2)

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

# Try DE
y <- DGEList(counts = RA_RNA_SEQ_SYNOVIAL_MEMBRANE_matrix, group = groups)

y <- calcNormFactors(y) # Should we do that?

y <- estimateDisp(y)
sqrt(y$common.dispersion) # biological coefficient of variation ??? Should the normal be ~0.4?

et <- exactTest(y)
results_edgeR <- topTags(et, n = nrow(RA_RNA_SEQ_SYNOVIAL_MEMBRANE_matrix))
results_of_DE <- results_edgeR$table

FDR_threshold = 0.05

sum(
  results_edgeR$table$FDR < FDR_threshold & 
  (results_edgeR$table$logFC >= 1 | results_edgeR$table$logFC <= -1)
) # why does he use .1 for FDR?

DE_genes = rownames(results_edgeR)[results_edgeR$table$FDR < FDR_threshold  & (results_edgeR$table$logFC >= 1 | results_edgeR$table$logFC <= -1)]

plotSmear(et, de.tags = rownames(results_edgeR)[results_edgeR$table$FDR < FDR_threshold  & (results_edgeR$table$logFC >= 1 | results_edgeR$table$logFC <= -1)])
abline(h = c(-1, 1), col = "blue")

# Experimenting with volcano plot
results_of_DE$topDE <- "NA"
results_of_DE$topDE[results_of_DE$logFC > 1 & results_of_DE$FDR < 0.05] <- "Up"
results_of_DE$topDE[results_of_DE$logFC < -1 & results_of_DE$FDR < 0.05] <- "Down"

ggplot(data=results_of_DE, aes(x=logFC, y=-log10(FDR), color = topDE)) +
  geom_point() +
  theme_minimal() +
  scale_colour_discrete(breaks = c("Up", "Down"))

