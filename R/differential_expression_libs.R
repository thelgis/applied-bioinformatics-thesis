run_edgeR <- function (input_csv, fdr_threshold = 0.05, logFC_threshold = 1) {

  dataset_with_metadata <- read.table(
    paste0("./data/r_pre_processed_datasets/", input_csv),
    header = TRUE,
    sep = ",",
    row.names = 1
  )

  # Determine groups vector
  groups_df <- dataset_with_metadata[c("Condition"),]
  groups <- as.character(groups_df[1,])
  print("Groups:")
  print(groups)

  # Remove metadata columns
  dataset <- dataset_with_metadata[1:(nrow(dataset_with_metadata) - 1), ] # There is only one metadata column i.e. 'Condition'

  # Convert to numeric matrix
  dataset_as_matrix <- as.matrix(sapply(dataset, as.numeric))
  row.names(dataset_as_matrix) <- row.names(dataset)

  # Run DE
  y <- DGEList(counts = dataset_as_matrix, group = groups)
  y <- calcNormFactors(y) # TODO Should we do that in our case?
  y <- estimateDisp(y)

  dispersion = sqrt(y$common.dispersion) # TODO biological coefficient of variation? Should the normal be ~0.4?
  print(paste("Dispersion for", input_csv, "is:", dispersion))

  et <- exactTest(y)
  results_edgeR <- topTags(et, n = nrow(dataset_as_matrix))

  number_of_DE_genes = sum(
    results_edgeR$table$FDR < fdr_threshold &
      (results_edgeR$table$logFC >= logFC_threshold | results_edgeR$table$logFC <= -logFC_threshold)
  ) # TODO why does he use .1 for FDR?

  print(paste("Number of DE genes in", input_csv, "is:", number_of_DE_genes))
  de_genes <- rownames(results_edgeR)[results_edgeR$table$FDR < fdr_threshold  & (results_edgeR$table$logFC >= logFC_threshold | results_edgeR$table$logFC <= -logFC_threshold)]

  results_table <- results_edgeR$table
  results_table$topDE <- "NA"
  results_table$topDE[results_table$logFC >= logFC_threshold & results_table$FDR < fdr_threshold] <- "Up"
  results_table$topDE[results_table$logFC <= -logFC_threshold & results_table$FDR < fdr_threshold] <- "Down"

  # Save results
  subpath <- substr(input_csv, 1, nchar(input_csv) - 4)
  path <- file.path("results", subpath)
  dir.create(path)

  # Save DE data
  write.csv(results_table, paste0(path, "/edgeR_results.csv"), row.names=FALSE)
  write.csv(de_genes, paste0(path, "/edgeR_de_genes.csv"), row.names=FALSE)

  # Save smear plot
  png(paste0(path,"/edgeR_smear_plot.png"))
  plotSmear(et, de.tags = rownames(results_edgeR)[results_edgeR$table$FDR < fdr_threshold  & (results_edgeR$table$logFC >= logFC_threshold | results_edgeR$table$logFC <= -logFC_threshold)])
  abline(h = c(-logFC_threshold, logFC_threshold), col = "blue")
  dev.off()

  # Save volcano plot
  # TODO there seems to be an issue here -> debug
  png(paste0(path,"/edgeR_volcano_plot.png"))
  ggplot(data=results_table, aes(x=logFC, y=-log10(FDR), color = topDE)) +
    geom_point() +
    theme_minimal() +
    scale_colour_discrete(breaks = c("Up", "Down"))
  dev.off()

}