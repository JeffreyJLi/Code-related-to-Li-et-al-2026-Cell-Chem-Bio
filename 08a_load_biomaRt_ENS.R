load_biomaRt_ENS <- function(output_bm_filename) {
  if (endsWith(output_bm_filename, ".Rds") != TRUE) {
    output_bm_filename <- paste0(output_bm_filename,".Rds")
  }
  if (file.exists(output_bm_filename)) {
    print(paste0(output_bm_filename, " file found. Loading .Rds object as 'bm'."))
    bm <- readRDS(output_bm_filename)
  } else {
    print(paste0(output_bm_filename, " file NOT found. Downloading ensembl gene names for human GRCh38.p14 genome using biomaRt"))
    # Function to get dataframe of human-readable gene names to relabel rows (which use the ensembl gene ids) later
    get_ENS <- function(dset="hsapiens_gene_ensembl") {
      mart <- useMart("ENSEMBL_MART_ENSEMBL", dataset = dset)
      attributes <- c("ensembl_gene_id","ensembl_gene_id_version", "external_gene_name",
                      "chromosome_name", "start_position", "end_position")
      bm <- getBM(attributes, mart = mart)
      return(bm)
    }
    bm <- get_ENS()
    
    # I saved the bm list object containing all the ensembl gene names for the human genes (GRCh38.p14).
    # This is because I had recent issues downloading the gene names using biomaRt (I suspect Ensembl had server issues)
    print("Saving .Rds file of human ensembl gene names pulled by biomaRt")
    saveRDS(bm, file = output_bm_filename)
  }
  print("Dropping rows in biomaRt data with no external_gene_names")
  # Drop any rows in bm dataframe with no external_gene_names
  bm <- dplyr::filter(bm, ! external_gene_name == '')
  return(bm)
}