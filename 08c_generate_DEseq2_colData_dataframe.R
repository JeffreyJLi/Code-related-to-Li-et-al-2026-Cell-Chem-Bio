generate_DESeq2_colData_dataframe <- function(featureCounts_dataframe) {
  # Create new "samples" dataframe of the conditions based on the column names in the count Matrix
  # This will serve as the "coldata" argument for using the DESeqDataSetFromMatrix function
  samples <- colnames(featureCounts_dataframe)
  samples <- as.data.frame(samples)
  
  samples <- samples %>% 
    # Add ZAKi column to samples dataframe based on if sample was +ZAKi or not.
    mutate(inhibitor = case_when(grepl("\\.[A-D]$", samples) ~ "noZAKi", 
                                 grepl("\\.[E-H]$", samples) ~ "ZAKi")) %>%
    # Add Replicate column to samples dataframe
    mutate(Rep = case_when(grepl("^Rep1\\.", samples) ~ 1, 
                           grepl("^Rep2\\.", samples) ~ 2, 
                           grepl("^Rep3\\.", samples) ~ 3)) %>% 
    # Add drug column to samples dataframe based on what drug was applied to each sample
    mutate(drug = case_when(grepl("\\.[A,E]", samples) ~ "untreated", 
                            grepl("\\.[B,F]", samples) ~ "ANS", 
                            grepl("\\.[C,G]", samples) ~ "EME", 
                            grepl("\\.[D,H]", samples) ~ "DDB"))
  
  # Assign rownames of samples as the samples column (i.e. the column names of the featureCounts_dataframe matrix)
  rownames(samples) <- samples$samples
  
  # Reassign samples$drug_treatment variable (which are character type) as factor type in specific order
  samples$drug <- factor(samples$drug, levels = c("untreated","ANS","EME","DDB"))
  samples$inhibitor <- factor(samples$inhibitor, levels = c("noZAKi","ZAKi"))
  
  return(samples)
}