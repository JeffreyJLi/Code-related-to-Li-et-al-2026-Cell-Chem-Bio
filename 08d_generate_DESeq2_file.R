generate_DESeq2_file <- function(processed_featureCounts, featureCounts_colData, output_DESeq2_filename) {
  if (endsWith(output_DESeq2_filename, ".Rds") != TRUE) {
    output_DESeq2_filename <- paste0(output_DESeq2_filename,".Rds")
  }
  # Check if DESeq2 output .R file exists already. If so, it loads it in. If not, it runs DESeq2.
  # I wrote this to save time since DESeq2 takes ~30 sec-1 min to run otherwise.
  if (file.exists(output_DESeq2_filename)) {
    print(paste0(output_DESeq2_filename, " file found. Loading .Rds object as 'dds'."))
    dds <- readRDS(output_DESeq2_filename)
  } else {
    print(paste0(output_DESeq2_filename," file NOT found. Running DESeq2."))
    # Function to get dataframe of human-readable gene names to relabel rows (which use the ensembl gene ids) later
    dds <- DESeqDataSetFromMatrix(countData = processed_featureCounts, colData=featureCounts_colData, design=~drug + inhibitor + drug:inhibitor)
    dds <- DESeq(dds)
    
    print("Saving .Rds file of DESeq2 output data")
    saveRDS(dds, file = output_DESeq2_filename)
  }
  return(dds)
}