#### 20260826 Written. Reviewer asked for GSEA analysis.
#### I had done some cursory analysis, but honestly I didn't know what I was doing back then. Let's do it right.
#### Followed this tutorial https://alexslemonade.github.io/refinebio-examples/03-rnaseq/pathway-analysis_rnaseq_02_gsea.html

library("tidyverse")
library("dplyr")
library("here")
library("clusterProfiler")
library("msigdbr")


# Import DESeq2 data
ANS_genes <- read_tsv(file = "script_outputs/deseq2_results/ANS_deseq2_results.tsv")
EME_genes <- read_tsv(file = "script_outputs/deseq2_results/EME_deseq2_results.tsv")
DDB_genes <- read_tsv(file = "script_outputs/deseq2_results/DDB_deseq2_results.tsv")

ANSZAKi_genes <- read_tsv(file = "script_outputs/deseq2_results/ANS+ZAKi_deseq2_results.tsv")
EMEZAKi_genes <- read_tsv(file = "script_outputs/deseq2_results/EME+ZAKi_deseq2_results.tsv")
DDBZAKi_genes <- read_tsv(file = "script_outputs/deseq2_results/DDB+ZAKi_deseq2_results.tsv")

#### GSEA analysis using clusterProfiler ---------
# Followed this tutorial https://alexslemonade.github.io/refinebio-examples/03-rnaseq/pathway-analysis_rnaseq_02_gsea.html

# Get list of all genes in each condition. Do NOT pre-filter.
# Convert to a named vector of Log2FC
ANS_DEG <-ANS_genes %>% drop_na(external_gene_name) %>% select(external_gene_name, log2FoldChange) %>% arrange(desc(log2FoldChange)) %>% deframe()
EME_DEG <- EME_genes %>% drop_na(external_gene_name) %>% select(external_gene_name, log2FoldChange) %>% arrange(desc(log2FoldChange)) %>% deframe()
DDB_DEG <- DDB_genes %>% drop_na(external_gene_name) %>% select(external_gene_name, log2FoldChange) %>% arrange(desc(log2FoldChange)) %>% deframe()

ANSZAKi_DEG <- ANSZAKi_genes %>% drop_na(external_gene_name) %>% select(external_gene_name, log2FoldChange) %>% arrange(desc(log2FoldChange)) %>% deframe()
EMEZAKi_DEG <- EMEZAKi_genes %>% drop_na(external_gene_name) %>% select(external_gene_name, log2FoldChange) %>% arrange(desc(log2FoldChange)) %>% deframe()
DDBZAKi_DEG <- DDBZAKi_genes  %>% drop_na(external_gene_name) %>% select(external_gene_name, log2FoldChange) %>% arrange(desc(log2FoldChange)) %>% deframe()

# GSEA for only the curated "hallmark" sets
human_hallmark_sets <- msigdbr(
  species = "Homo sapiens", 
  collection = "H"
)

# Write function to abstract GSEA call for readability
run_gsea <- function(gene_list) {
  # Set seed for reproducibility
  set.seed(8888)
  
  GSEA(
    geneList = gene_list, # Ordered ranked gene list
    minGSSize = 25, # Minimum gene set size
    maxGSSize = 500, # Maximum gene set set
    pvalueCutoff = 0.05, # p-value cutoff
    eps = 0, # Boundary for calculating the p value
    seed = TRUE, # Set seed to make results reproducible
    pAdjustMethod = "BH", # Benjamini-Hochberg correction
    TERM2GENE = dplyr::select(
      human_hallmark_sets,
      gs_name,
      gene_symbol))
}

ANS_gsea_results <- run_gsea(ANS_DEG)
# We can access the results from our `gsea_results` object using `@result`
ANS_gsea_results <- ANS_gsea_results@result

EME_gsea_results <- run_gsea(EME_DEG)
# We can access the results from our `gsea_results` object using `@result`
EME_gsea_results <- EME_gsea_results@result

DDB_gsea_results <- run_gsea(DDB_DEG)
# We can access the results from our `gsea_results` object using `@result`
DDB_gsea_results <- DDB_gsea_results@result

ANSZAKi_gsea_results <- run_gsea(ANSZAKi_DEG)
# We can access the results from our `gsea_results` object using `@result`
ANSZAKi_gsea_results <- ANSZAKi_gsea_results@result

EMEZAKi_gsea_results <- run_gsea(EMEZAKi_DEG)
# We can access the results from our `gsea_results` object using `@result`
EMEZAKi_gsea_results <- EMEZAKi_gsea_results@result

DDBZAKi_gsea_results <- run_gsea(DDBZAKi_DEG)
# We can access the results from our `gsea_results` object using `@result`
DDBZAKi_gsea_results <- DDBZAKi_gsea_results@result


#### Write GSEA results to .tsv
output_dir <- "script_outputs/reviewer_analysis/gsea/Hallmark_sets/"

write_tsv(ANS_gsea_results, paste0(output_dir, "ANS_gsea_results_hallmark.tsv"))
write_tsv(EME_gsea_results,  paste0(output_dir, "EME_gsea_results_hallmark.tsv"))
write_tsv(DDB_gsea_results,  paste0(output_dir, "DDB_gsea_result_hallmarks.tsv"))
write_tsv(ANSZAKi_gsea_results, paste0(output_dir,"ANSZAKi_gsea_results_hallmark.tsv"))
write_tsv(EMEZAKi_gsea_results, paste0(output_dir,"EMEZAKi_gsea_results_hallmark.tsv"))
write_tsv(DDBZAKi_gsea_results,  paste0(output_dir, "DDBZAKi_gsea_results_hallmark.tsv"))

###### ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
##### # GSEA for only the curated GO + HPO sets
human_hallmark_sets <- msigdbr(
  species = "Homo sapiens", 
  collection = "C5"
)

# Write function to abstract GSEA call for readability
run_gsea <- function(gene_list) {
  # Set seed for reproducibility
  set.seed(8888)
  
  GSEA(
    geneList = gene_list, # Ordered ranked gene list
    minGSSize = 25, # Minimum gene set size
    maxGSSize = 500, # Maximum gene set set
    pvalueCutoff = 0.05, # p-value cutoff
    eps = 0, # Boundary for calculating the p value
    seed = TRUE, # Set seed to make results reproducible
    pAdjustMethod = "BH", # Benjamini-Hochberg correction
    TERM2GENE = dplyr::select(
      human_hallmark_sets,
      gs_name,
      gene_symbol))
}

ANS_gsea_results <- run_gsea(ANS_DEG)
# We can access the results from our `gsea_results` object using `@result`
ANS_gsea_results <- ANS_gsea_results@result

EME_gsea_results <- run_gsea(EME_DEG)
# We can access the results from our `gsea_results` object using `@result`
EME_gsea_results <- EME_gsea_results@result

DDB_gsea_results <- run_gsea(DDB_DEG)
# We can access the results from our `gsea_results` object using `@result`
DDB_gsea_results <- DDB_gsea_results@result

ANSZAKi_gsea_results <- run_gsea(ANSZAKi_DEG)
# We can access the results from our `gsea_results` object using `@result`
ANSZAKi_gsea_results <- ANSZAKi_gsea_results@result

EMEZAKi_gsea_results <- run_gsea(EMEZAKi_DEG)
# We can access the results from our `gsea_results` object using `@result`
EMEZAKi_gsea_results <- EMEZAKi_gsea_results@result

DDBZAKi_gsea_results <- run_gsea(DDBZAKi_DEG)
# We can access the results from our `gsea_results` object using `@result`
DDBZAKi_gsea_results <- DDBZAKi_gsea_results@result


#### Write GSEA results to .tsv
output_dir <- "script_outputs/reviewer_analysis/gsea/C5_sets/"

write_tsv(ANS_gsea_results, paste0(output_dir, "ANS_gsea_results_C5collection.tsv"))
write_tsv(EME_gsea_results,  paste0(output_dir, "EME_gsea_results_C5collection.tsv"))
write_tsv(DDB_gsea_results,  paste0(output_dir, "DDB_gsea_results_C5collection.tsv"))
write_tsv(ANSZAKi_gsea_results, paste0(output_dir,"ANSZAKi_gsea_results_C5collection.tsv"))
write_tsv(EMEZAKi_gsea_results, paste0(output_dir,"EMEZAKi_gsea_results_C5collection.tsv"))
write_tsv(DDBZAKi_gsea_results,  paste0(output_dir, "DDBZAKi_gsea_results_C5collection.tsv"))