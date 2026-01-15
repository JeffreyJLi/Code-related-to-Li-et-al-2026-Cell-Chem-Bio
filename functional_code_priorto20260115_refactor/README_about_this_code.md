# What the heck happened with this code on 20260115?

The code that is functional and used to generate the publication figures for my JL318 RNA-seq analysis is found in the "08c_featureCounts_multiVariate_DESeq.R" file. However, it was a massive script containing so many functions that made it difficult to actually understand what was going on. As a record, this code is kept in the /functional_code_priorto20260115_refactor directory.

On 20260115, I decided to refactor this "08c_featureCounts_multiVariate_DESeq.R" file into separate .R scripts, with each script containing self-contained functions that would be called by the main "08_main_multiVariate_DESeq2.R" file. This greatly reduced the number of lines in the main analysis file, which has made things way easier to follow.