library(Rsubread)

sample_file_list <- grep('.sortedByCoord.out.bam$',list.files("../05_STAR_alignments"),value=TRUE)
sample_file_list <- paste("../05_STAR_alignments/",sample_file_list, sep="")
sample_file_list

featureCounts <- Rsubread::featureCounts(files = c(sample_file_list), nthreads=36, isPairedEnd = TRUE, countReadPairs = TRUE, strandSpecific=1, minMQS=10, isGTFAnnotationFile=TRUE, GTF.attrType = "gene_id", annot.ext="../gencode_human_genome/gencode.v46.primary_assembly.annotation.gtf")
str(Rsubread::featureCounts)
featureCounts$counts

saveRDS(featureCounts, "featureCounts_R.rds")

featureCounts <- readRDS("06d_featureCounts_R.rds")
featureCounts$count
