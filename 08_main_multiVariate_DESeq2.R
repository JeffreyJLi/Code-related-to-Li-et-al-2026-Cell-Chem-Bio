library("tidyverse")
library("DESeq2")
library("dplyr")
library("ggrepel")
library("biomaRt")
library("conflicted")
library("ggvenn")
library("here")
library("plotly")

# Prefer dplyr::select() over biomaRt::select()
conflicted::conflicts_prefer(dplyr::select)

# Prefer dplyr::filter() over stats::filter()
conflicted::conflicts_prefer(dplyr::filter)


#### Prepare featureCounts data and run DESeq2 -----------------------------------

source("08b_process_featureCounts.R")

featureCounts <- process_featureCounts("../06c_featureCounts_countMatrix/JL318_featureCounts_countMatrix.txt")

source("08c_generate_DEseq2_colData_dataframe.R")

samples <- generate_DESeq2_colData_dataframe(featureCounts)

# Verifies if the rownames of the samples dataframe (i.e. the ColData input for DESeq2) are in the same order as the featureCounts columns.
# Per the DESeq2 vignette, this is REQUIRED for proper DESeq2 function. Otherwise, things break!
stopifnot(all(rownames(samples) == colnames(featureCounts)))

source("08d_generate_DESeq2_file.R")

dds <- generate_DESeq2_file(featureCounts, samples, output_DESeq2_filename = "20250611_JL318_raw_DESeq2_output.Rds")

# Removes rows with less than 100 reads, but only for 3 samples (the smallest group size - each condition has 3 replicates)
# I previously dropped rows with less than 100 reads in ALL 24 samples, which led to me losing major ZAK-responsive genes such as FOS and JUNB.
idx <- rowSums(counts(dds, normalized=TRUE) >= 100 ) >= 3
dds <- dds[idx,]

res <- results(dds)

####  DESeq2 diagnostic plots ------------------------------------------------------

# Plots dispersion plot for diagnostics. Nicolle mentioned that the fit should be as flat as possible (it tends to tail up for <100 mean of normalized counts).
# Removing rows with less than 100 reads will help.
plotDispEsts(dds)
plotMA(dds)
plotMA(dds, ylim = c(-5,5))
resLFC <- lfcShrink(dds, coef="drug_ANS_vs_untreated", type="apeglm")
plotMA(resLFC)


####  Plot PCA plot -----------------------------------------------------------

# Generate PCA plot of the DESeq results
source("08e_plot_PCAplot.R")

# Return DESeq2 PCAplot data as a separate dataframe for custom ggplot
vsd <- vst(dds,blind=FALSE)
PCAplot_data <- plotPCA(vsd,intgroup=c("drug", "inhibitor"), returnData=TRUE)
percentVar <- round(100 * attr(PCAplot_data, "percentVar"))

plot_PCAplot(PCAplot_data, percentVar)

# ggsave("script_outputs/20260303_JL318_PCA_plot.svg", plot = plot_PCAplot(PCAplot_data, percentVar), width = 7, height = 5, unit = "in")


####  Plot DESeq2 data --------------------------------------------------------

source("08f_plotting_functions.R")

# SEE AFFINITY FILE IN JL318/scripts/ignore for experimental design explanation

ANS_resultsnames <- c("drug_ANS_vs_untreated", "inhibitor_ZAKi_vs_noZAKi", "drugANS.inhibitorZAKi")
EME_resultsnames <- c("drug_EME_vs_untreated", "inhibitor_ZAKi_vs_noZAKi", "drugEME.inhibitorZAKi")
DDB_resultsnames <- c("drug_DDB_vs_untreated", "inhibitor_ZAKi_vs_noZAKi", "drugDDB.inhibitorZAKi")

ANS_results <- wrapper_generate_results_dataframes(dds, ANS_resultsnames)
EME_results <- wrapper_generate_results_dataframes(dds, EME_resultsnames)
DDB_results <- wrapper_generate_results_dataframes(dds, DDB_resultsnames)

# This is equivalent to -ZAKi[ANS/UT]
plot_DESeq2_volcano(ANS_results$drug_vs_untreated, plot_title = "ANS treatment")
plot_DESeq2_volcano(EME_results$drug_vs_untreated, plot_title = "EME treatment")
plot_DESeq2_volcano(DDB_results$drug_vs_untreated, plot_title = "DDB treatment")

# This is equivalent to +ZAKi[ANS/UT]
plot_DESeq2_volcano(ANS_results$drug_inhibitor_vs_untreated_inhibitor, plot_title = "ANS+ZAKi (ANS with ZAKi pre-treatment)")
plot_DESeq2_volcano(EME_results$drug_inhibitor_vs_untreated_inhibitor, plot_title = "EME+ZAKi (EME with ZAKi pre-treatment)")
plot_DESeq2_volcano(DDB_results$drug_inhibitor_vs_untreated_inhibitor, plot_title = "DDB+ZAKi (DDB with ZAKi pre-treatment)")

# This is equivalent to UT[+ZAKi/-ZAKi]
plot_DESeq2_volcano(ANS_results$untreated_inhibitor_vs_untreated, plot_title = "UT+ZAKi / UT")
plot_DESeq2_volcano(EME_results$untreated_inhibitor_vs_untreated, plot_title = "UT[+ZAKi/-ZAKi]")
plot_DESeq2_volcano(DDB_results$untreated_inhibitor_vs_untreated)

# This is equivalent to ANS[+ZAKi/-ZAKi]
# This is the ZAKi effect in ANS-treated cells
plot_DESeq2_volcano(ANS_results$drug_inhibitor_vs_drug, xaxis = c(-10,0))
plot_DESeq2_volcano(EME_results$drug_inhibitor_vs_drug, xaxis = c(-2,2), plot_title = "EME+ZAKi / EME")
plot_DESeq2_volcano(DDB_results$drug_inhibitor_vs_drug, xaxis = c(-10,0))

# This is equivalent to ANS[+ZAKi/-ZAKi]/UT[+ZAKi/-ZAKi]
# Do the ZAKi effect different in EME-treated cells compared to UT cells?
# This ignores the contribution of ZAKi treatment in UT cells
plot_DESeq2_volcano(ANS_results$drug_inhibitor_interaction, xaxis = c(-10,5))
plot_DESeq2_volcano(EME_results$drug_inhibitor_interaction, xaxis = c(-2,2), plot_title = "[EME+ZAKi / EME)] / [UT+ZAKi / UT)]")
plot_DESeq2_volcano(DDB_results$drug_inhibitor_interaction, xaxis = c(-10,5))


####  Saving .html interactive volcano plots ----------------------------------

# This is equivalent to -ZAKi[ANS/UT]
save_interactive_plot_as_html(plot_DESeq2_volcano_interactive(ANS_results$drug_vs_untreated, plot_title = "ANS treatment"), output_dir = "interactive_volcanos")
save_interactive_plot_as_html(plot_DESeq2_volcano_interactive(EME_results$drug_vs_untreated, plot_title = "EME treatment"), output_dir = "interactive_volcanos")
save_interactive_plot_as_html(plot_DESeq2_volcano_interactive(DDB_results$drug_vs_untreated, plot_title = "DDB treatment"), output_dir = "interactive_volcanos")

# This is equivalent to +ZAKi[ANS/UT]
save_interactive_plot_as_html(plot_DESeq2_volcano_interactive(ANS_results$drug_inhibitor_vs_untreated_inhibitor, plot_title = "ANS+ZAKi treatment"), output_dir = "interactive_volcanos")
save_interactive_plot_as_html(plot_DESeq2_volcano_interactive(EME_results$drug_inhibitor_vs_untreated_inhibitor, plot_title = "EME+ZAKi treatment"), output_dir = "interactive_volcanos")
save_interactive_plot_as_html(plot_DESeq2_volcano_interactive(DDB_results$drug_inhibitor_vs_untreated_inhibitor, plot_title = "DDB+ZAKi treatment"), output_dir = "interactive_volcanos")


####  What # of genes are significantly up-regulated? -------------------------

# Nicolle used p < 0.01 for her analysis in her paper.
# This is equivalent to -ZAKi[ANS/UT]
ANS_results$drug_vs_untreated %>% filter((padj < 0.01) & (log2FoldChange > 1)) %>% nrow()
EME_results$drug_vs_untreated %>% filter((padj < 0.01) & (log2FoldChange > 1)) %>% nrow()
DDB_results$drug_vs_untreated %>% filter((padj < 0.01) & (log2FoldChange > 1)) %>% nrow()

# This is equivalent to +ZAKi[ANS/UT]
ANS_results$drug_inhibitor_vs_untreated_inhibitor %>% filter((padj < 0.01) & (log2FoldChange > 1)) %>% nrow()
EME_results$drug_inhibitor_vs_untreated_inhibitor %>% filter((padj < 0.01) & (log2FoldChange > 1)) %>% nrow()
DDB_results$drug_inhibitor_vs_untreated_inhibitor %>% filter((padj < 0.01) & (log2FoldChange > 1)) %>% nrow()


####  Generating .tsv files of DESeq2 data (fold changes, padj, etc.) ---------

write_tsv(ANS_results$drug_vs_untreated, file = "script_outputs/deseq2_results/ANS_deseq2_results.tsv")
write_tsv(EME_results$drug_vs_untreated, file = "script_outputs/deseq2_results/EME_deseq2_results.tsv")
write_tsv(DDB_results$drug_vs_untreated, file = "script_outputs/deseq2_results/DDB_deseq2_results.tsv")

write_tsv(ANS_results$drug_inhibitor_vs_untreated_inhibitor, file = "script_outputs/deseq2_results/ANS+ZAKi_deseq2_results.tsv")
write_tsv(EME_results$drug_inhibitor_vs_untreated_inhibitor, file = "script_outputs/deseq2_results/EME+ZAKi_deseq2_results.tsv")
write_tsv(DDB_results$drug_inhibitor_vs_untreated_inhibitor, file = "script_outputs/deseq2_results/DDB+ZAKi_deseq2_results.tsv")


####  Generating cross-correlation plots for each condition -------------------

# Cross-correlations
plot_crossCorrelation(ANS_results$drug_vs_untreated, DDB_results$drug_vs_untreated, "ANS/UT","DDB/UT")
plot_crossCorrelation(ANS_results$drug_vs_untreated, EME_results$drug_vs_untreated, "ANS/UT","EME/UT")
plot_crossCorrelation(ANS_results$drug_vs_untreated, ANS_results$drug_inhibitor_vs_untreated_inhibitor,"ANS/UT", "ANS+ZAKi / UT+ZAKi")
plot_crossCorrelation(EME_results$drug_vs_untreated, ANS_results$drug_inhibitor_vs_untreated_inhibitor,"EME/UT", "ANS+ZAKi / UT+ZAKi")

# Cross-correlation plots for RNA salon talk (no labels)
# Save using RStudio's export as .pdf file feature to JL318/FIGURES folder.
RNA_salon_plot_crossCorrelation(ANS_results$drug_vs_untreated, EME_results$drug_vs_untreated, "ANS/UT","EME/UT", plot_title = "ANS- vs EME-treated samples")
RNA_salon_plot_crossCorrelation_nolabel(ANS_results$drug_vs_untreated, EME_results$drug_vs_untreated, "ANS/UT","EME/UT", plot_title = "ANS- vs EME-treated samples")

RNA_salon_plot_crossCorrelation(ANS_results$drug_vs_untreated, ANS_results$drug_inhibitor_vs_untreated_inhibitor,"ANS/UT", "ANS+ZAKi / UT+ZAKi", plot_title = "ANS- vs ANS+ZAKi-treated samples")
RNA_salon_plot_crossCorrelation_nolabel(ANS_results$drug_vs_untreated, ANS_results$drug_inhibitor_vs_untreated_inhibitor,"ANS/UT", "ANS+ZAKi / UT+ZAKi", plot_title = "ANS- vs ANS+ZAKi-treated samples")

RNA_salon_plot_crossCorrelation(EME_results$drug_vs_untreated, ANS_results$drug_inhibitor_vs_untreated_inhibitor, "EME/UT", "ANS+ZAKi / UT+ZAKi", plot_title = "EME- vs ANS+ZAKi-treated samples")

# Save using RStudio's export as .pdf file feature to JL318/FIGURES folder.

# Cross-correlation plots for paper
# Save using RStudio's export as .pdf file feature to JL318/FIGURES folder. SET AS SQUARE dimensions prior to saving!!


plot_crossCorrelation(ANS_results$drug_vs_untreated, DDB_results$drug_vs_untreated, "ANS/UT","DDB/UT", plot_title = "ANS- vs DDB-treated samples")
plot_crossCorrelation(ANS_results$drug_vs_untreated, EME_results$drug_vs_untreated, "ANS/UT","EME/UT", plot_title = "ANS- vs EME-treated samples")
plot_crossCorrelation(DDB_results$drug_vs_untreated, EME_results$drug_vs_untreated, "DDB/UT","EME/UT", plot_title = "DDB- vs EME-treated samples")

plot_crossCorrelation(ANS_results$drug_vs_untreated, ANS_results$drug_inhibitor_vs_untreated_inhibitor,"ANS/UT", "ANS+ZAKi / UT+ZAKi", plot_title = "ANS- vs ANS+ZAKi-treated samples")
plot_crossCorrelation(DDB_results$drug_vs_untreated, DDB_results$drug_inhibitor_vs_untreated_inhibitor,"DDB/UT", "DDB+ZAKi / UT+ZAKi", plot_title = "DDB- vs DDB+ZAKi-treated samples")
plot_crossCorrelation(EME_results$drug_vs_untreated, EME_results$drug_inhibitor_vs_untreated_inhibitor,"EME/UT", "EME+ZAKi / UT+ZAKi", plot_title = "EME- vs EME+ZAKi-treated samples")

plot_crossCorrelation(EME_results$drug_vs_untreated, ANS_results$drug_inhibitor_vs_untreated_inhibitor, "EME/UT", "ANS+ZAKi / UT+ZAKi", plot_title = "EME- vs ANS+ZAKi-treated samples")
plot_crossCorrelation(EME_results$drug_vs_untreated, DDB_results$drug_inhibitor_vs_untreated_inhibitor, "EME/UT", "DDB+ZAKi / UT+ZAKi", plot_title = "EME- vs DDB+ZAKi-treated samples")
plot_crossCorrelation(EME_results$drug_inhibitor_vs_untreated_inhibitor, ANS_results$drug_inhibitor_vs_untreated_inhibitor, "EME+ZAKi/UT", "ANS+ZAKi / UT+ZAKi", plot_title = "EME+ZAKi- vs ANS+ZAKi-treated samples")

# Save using RStudio's export as .pdf file feature to JL318/FIGURES folder.


####  Saving .html interactive cross-correlation plots ------------------------

plot_crossCorrelation_interactive(ANS_results$drug_vs_untreated, DDB_results$drug_vs_untreated, "ANS/UT","DDB/UT", plot_title = "ANS vs DDB cross-correlation")
plot_crossCorrelation_interactive(ANS_results$drug_vs_untreated, EME_results$drug_vs_untreated, "ANS/UT","EME/UT", plot_title = "ANS vs EME cross-correlation")
plot_crossCorrelation_interactive(ANS_results$drug_vs_untreated, ANS_results$drug_inhibitor_vs_untreated_inhibitor,"ANS/UT", "ANS+ZAKi / UT+ZAKi", plot_title = "ANS vs ANS+ZAKi cross-correlation")
plot_crossCorrelation_interactive(EME_results$drug_vs_untreated, ANS_results$drug_inhibitor_vs_untreated_inhibitor,"EME/UT", "ANS+ZAKi / UT+ZAKi", plot_title = "EME vs ANS+ZAKi cross-correlation")

# Save as .html file
save_interactive_plot_as_html(plot_crossCorrelation_interactive(ANS_results$drug_vs_untreated, DDB_results$drug_vs_untreated, "ANS/UT","DDB/UT", plot_title = "ANS vs DDB cross-correlation"), output_dir = "interactive_cross-correlations")
save_interactive_plot_as_html(plot_crossCorrelation_interactive(ANS_results$drug_vs_untreated, EME_results$drug_vs_untreated, "ANS/UT","EME/UT", plot_title = "ANS vs EME cross-correlation"), output_dir = "interactive_cross-correlations")
save_interactive_plot_as_html(plot_crossCorrelation_interactive(ANS_results$drug_vs_untreated, ANS_results$drug_inhibitor_vs_untreated_inhibitor,"ANS/UT", "ANS+ZAKi / UT+ZAKi", plot_title = "ANS vs ANS+ZAKi cross-correlation"), output_dir = "interactive_cross-correlations")
save_interactive_plot_as_html(plot_crossCorrelation_interactive(EME_results$drug_vs_untreated, ANS_results$drug_inhibitor_vs_untreated_inhibitor,"EME/UT", "ANS+ZAKi / UT+ZAKi", plot_title = "EME vs ANS+ZAKi cross-correlation"), output_dir = "interactive_cross-correlations")



#### Identifying what genes are ZAK-dependent and EDF1-dependent -------------

# Identifying ZAK-dependent genes (ZAKi-responsive)
# I played with the Log2FC and padj cutoffs to get the best overlap between ZAK-dependent genes and EDF1-dependent genes
ANS_ZAKi_dep_genes <- ANS_results$drug_inhibitor_vs_drug %>% filter(log2FoldChange < -2, padj_padded<1e-4) %>% select(ensembl_gene_id, external_gene_name, log2FoldChange, padj_padded)
EME_ZAKi_dep_genes <- EME_results$drug_inhibitor_vs_drug %>% filter(log2FoldChange < 0, padj_padded<1e-4) %>% select(ensembl_gene_id, external_gene_name, log2FoldChange, padj_padded)
DDB_ZAKi_dep_genes <- DDB_results$drug_inhibitor_vs_drug %>% filter(log2FoldChange < -2, padj_padded<1e-4) %>% select(ensembl_gene_id, external_gene_name, log2FoldChange, padj_padded)

# Find which ZAK-dependent genes are present in BOTH ANS + DDB
# This would be the most stringent
common_ZAKi_genes <- inner_join(ANS_ZAKi_dep_genes, DDB_ZAKi_dep_genes, by=c("ensembl_gene_id", "external_gene_name"), suffix = c(".ANS",".DDB"))

# Get EDF1-dependent genes (generated in a different RStudio project)
EDF1_dep_genes <- readRDS(file = "/Volumes/home/jli333/sinha2020_rnaseq/EDF1_dep_genes.rds")
EDF1_vs_WT_EL120 <- readRDS(file = "/Volumes/home/jli333/sinha2020_rnaseq/EDF1_vs_WT_EL120.rds")

# cleaning up dataset for ease of viewing
EDF1_genes <- EDF1_dep_genes %>% select(external_gene_name, log2FoldChange, padj_padded)
DDB_ZAKi_genes <- DDB_ZAKi_dep_genes %>% select(external_gene_name, log2FoldChange, padj_padded)

# What genes are both ZAK-dependent and EDF1-dependent
View(inner_join(EDF1_genes, DDB_ZAKi_genes, by="external_gene_name", suffix = c(".EDF1",".ZAKi")))

# Cross-correlation and Venn diagrams with EDF1-dependent and ZAK-dependent genes defined above

EDF1_ZAK_dep_genes_join <- full_join(EDF1_genes, DDB_ZAKi_genes, by = "external_gene_name", suffix = c(".EDF1", ".ZAK"))
EDF1_ZAK_dep_genes_join <- EDF1_ZAK_dep_genes_join %>% mutate(up_rows = ifelse(abs(log2FoldChange.EDF1) > 0 & abs(log2FoldChange.ZAK) >0, external_gene_name, NA))

# Cross-correlation plot for plotting ZAK-dependent genes vs EDF1-dependent genes
ggplot(EDF1_ZAK_dep_genes_join, aes(x=log2FoldChange.EDF1,
                                    y=log2FoldChange.ZAK)) +
  ylim(-5,5) +
  xlim(-5,5) +
  geom_point(alpha=0.5, pch=16, color="red") +
  # geom_text_repel(data=up_rows, aes(label=external_gene_name), point.padding=0.25, max.overlaps=15, size = 5) +
  geom_text_repel(data=EDF1_ZAK_dep_genes_join, aes(label=up_rows), point.padding=0.25, max.overlaps=25, size = 5) +
  geom_hline(yintercept=0, linetype="dashed", color = "black") +
  geom_vline(xintercept=0, linetype="dashed", color = "black") +
  geom_abline(intercept=0, slope=1, linetype="dashed") +
  theme_bw() +
  theme(axis.text.y = element_text(size = 18), 
        axis.text.x = element_text(size = 18), 
        axis.title = element_text(size = 20, face="bold"), 
        plot.title = element_text(size = 24, face = "bold"),
        plot.subtitle = element_text(size = 18))
#+ annotate("label",x = -4, y = 4, label = "EDF1-dependent genes = 65")

# Generate Venn Diagram of ZAK-dependent genes vs EDF1-dependent genes
ggvenn(list("EDF1 dep genes \n (Log2FC < 1, padj < 1e-4)" = EDF1_genes$external_gene_name,
            "ZAK_dep \n (Log2FC < -2, padj < 1e-4) " = common_ZAKi_genes$external_gene_name), auto_scale =TRUE, show_percentage = FALSE, set_name_size = 4) + coord_cartesian(xlim = c(-2, 2), ylim = c(-1.5, 1.5))


#### Plotting the log2FC of top ZAK-dep genes in ANS/DDB and EDF1 datasets --------

# Selecting the top (lowest log2FC in ANS_results$drug_inhibitor_vs_drug dataset) ZAK-dependent genes
ZAK_genes <- c("EGR1","DUSP1","FOSB","JUN","ATF3","IER2", "SOX4")

ANS_ZAK_top <- ANS_results$drug_inhibitor_vs_drug %>% filter(external_gene_name %in% ZAK_genes) %>% select(external_gene_name,log2FoldChange) %>% mutate(condition = "ANS+ZAKi / ANS") %>% mutate(external_gene_name = as.factor(external_gene_name))
EDF1_ZAK_top <- EDF1_genes %>% filter(external_gene_name %in% ZAK_genes) %>% select(external_gene_name,log2FoldChange) %>% mutate(condition = "∆EDF1/WT") %>% mutate(external_gene_name = as.factor(external_gene_name))
DDB_ZAK_top <- DDB_results$drug_inhibitor_vs_drug %>% filter(external_gene_name %in% ZAK_genes) %>% select(external_gene_name,log2FoldChange) %>% mutate(condition = "DDB+ZAKi / DDB") %>% mutate(external_gene_name = as.factor(external_gene_name))
EME_ZAK_top <- EME_results$drug_inhibitor_vs_drug %>% filter(external_gene_name %in% ZAK_genes) %>% select(external_gene_name,log2FoldChange) %>% mutate(condition = "EME+ZAKi / EME") %>% mutate(external_gene_name = as.factor(external_gene_name))


ZAK_top_genes_only <- full_join(ANS_ZAK_top, EDF1_ZAK_top, by = c("external_gene_name", "log2FoldChange", "condition"), suffix = c(".ANS",".EDF1"))
ZAK_top_genes_only <- full_join(ZAK_top_genes_only, DDB_ZAK_top, by = c("external_gene_name", "log2FoldChange", "condition"))
ZAK_top_genes_only <- full_join(ZAK_top_genes_only, EME_ZAK_top, by = c("external_gene_name", "log2FoldChange", "condition"))

ZAK_top_genes_only

ggplot(ZAK_top_genes_only, aes(x = fct_reorder(external_gene_name, log2FoldChange), y = log2FoldChange, fill = condition)) +
  geom_col(position = position_dodge(preserve = "single"), width = 0.6) +
  theme_bw() +
  theme(axis.text.y = element_text(size = 18), 
        axis.text.x = element_text(size = 18), 
        axis.title = element_text(size = 20, face="bold"), 
        plot.title = element_text(size = 20, face = "bold"),
        plot.subtitle = element_text(size = 18),
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 16)) +
  labs(x = "Gene name",
       y = expression(bold("log"[2]*"FoldChange")),
       fill = "Log2FC condition",
       title = "Log2FC of top ZAK-dependent genes across +ZAKi and ∆EDF1 datasets"
  )



ANS_ZAK_top <- ANS_results$drug_inhibitor_interaction %>% filter(external_gene_name %in% ZAK_genes) %>% select(external_gene_name,log2FoldChange) %>% mutate(condition = "ANS [+ZAKi / -ZAKi] / UT [+ZAKi / -ZAKi]") %>% mutate(external_gene_name = as.factor(external_gene_name))
EDF1_ZAK_top <- EDF1_genes %>% filter(external_gene_name %in% ZAK_genes) %>% select(external_gene_name,log2FoldChange) %>% mutate(condition = "∆EDF1/WT") %>% mutate(external_gene_name = as.factor(external_gene_name))
DDB_ZAK_top <- DDB_results$drug_inhibitor_interaction %>% filter(external_gene_name %in% ZAK_genes) %>% select(external_gene_name,log2FoldChange) %>% mutate(condition = "DDB [+ZAKi / -ZAKi] / UT [+ZAKi / -ZAKi]") %>% mutate(external_gene_name = as.factor(external_gene_name))
EME_ZAK_top <- EME_results$drug_inhibitor_interaction %>% filter(external_gene_name %in% ZAK_genes) %>% select(external_gene_name,log2FoldChange) %>% mutate(condition = "EME [+ZAKi / -ZAKi] / UT [+ZAKi / -ZAKi]") %>% mutate(external_gene_name = as.factor(external_gene_name))

ZAK_top_genes_only <- full_join(ANS_ZAK_top, EDF1_ZAK_top, by = c("external_gene_name", "log2FoldChange", "condition"), suffix = c(".ANS",".EDF1"))
ZAK_top_genes_only <- full_join(ZAK_top_genes_only, DDB_ZAK_top, by = c("external_gene_name", "log2FoldChange", "condition"))
ZAK_top_genes_only <- full_join(ZAK_top_genes_only, EME_ZAK_top, by = c("external_gene_name", "log2FoldChange", "condition"))

ZAK_top_genes_only

ggplot(ZAK_top_genes_only, aes(x = fct_reorder(external_gene_name, log2FoldChange), y = log2FoldChange, fill = condition)) +
  geom_col(position = position_dodge(preserve = "single"), width = 0.6) +
  theme_bw() +
  theme(axis.text.y = element_text(size = 18), 
        axis.text.x = element_text(size = 18), 
        axis.title = element_text(size = 20, face="bold"), 
        plot.title = element_text(size = 20, face = "bold"),
        plot.subtitle = element_text(size = 18),
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 16)) +
  labs(x = "Gene name",
       y = expression(bold("log"[2]*"FoldChange")),
       fill = "Log2FC condition",
       title = "Log2FC of top ZAK-dependent genes across +ZAKi and ∆EDF1 datasets"
  )

View(EDF1_ZAK_dep_genes_join %>% filter(is.na(log2FoldChange.EDF1) == FALSE & is.na(log2FoldChange.ZAK) == FALSE))

#### Generate ZAK-independent genes in ANS+ZAKi vs EME+ZAKi vs DDB+ZAKi for Venn Diagrams ---------

ANS_genes <- ANS_results$drug_inhibitor_vs_untreated_inhibitor %>% filter((padj < 0.01) & (log2FoldChange > 1)) %>% select(ensembl_gene_id)
EME_genes <- EME_results$drug_inhibitor_vs_untreated_inhibitor %>% filter((padj < 0.01) & (log2FoldChange > 1)) %>% select(ensembl_gene_id)
DDB_genes <- DDB_results$drug_inhibitor_vs_untreated_inhibitor %>% filter((padj < 0.01) & (log2FoldChange > 1)) %>% select(ensembl_gene_id)

write_tsv(ANS_genes, file = "script_outputs/ANS+ZAKi_genes_lfc1.tsv")
write_tsv(EME_genes, file = "script_outputs/EME+ZAKi_genes_lfc1.tsv")
write_tsv(DDB_genes, file = "script_outputs/DDB+ZAKi_genes_lfc1.tsv")

# Used DeepVenn webtool to visualize data

#### Generate gene lists in ANS vs EME vs DDB datasets for Venn Diagram ---------

ANS_genes <- ANS_results$drug_vs_untreated %>% filter((padj < 0.01) & (log2FoldChange > 1)) %>% select(ensembl_gene_id)
EME_genes <- EME_results$drug_vs_untreated %>% filter((padj < 0.01) & (log2FoldChange > 1)) %>% select(ensembl_gene_id)
DDB_genes <- DDB_results$drug_vs_untreated %>% filter((padj < 0.01) & (log2FoldChange > 1)) %>% select(ensembl_gene_id)

write_tsv(ANS_genes, file = "script_outputs/ANS_genes_lfc1.tsv")
write_tsv(EME_genes, file = "script_outputs/EME_genes_lfc1.tsv")
write_tsv(DDB_genes, file = "script_outputs/DDB_genes_lfc1.tsv")

# Used DeepVenn webtool to visualize data

#### Generate top 50 genes for Broad CLUE querying ---------

write_tsv(EME_results$drug_inhibitor_vs_untreated_inhibitor %>% select(external_gene_name, log2FoldChange) %>% arrange(desc(log2FoldChange)) %>% filter(is.na(external_gene_name) != TRUE) %>% slice_max(100), "EME+ZAKi_top_100.tsv")
write_tsv(EME_results$drug_inhibitor_vs_untreated_inhibitor %>% select(external_gene_name, log2FoldChange) %>% arrange(desc(log2FoldChange)) %>% filter(is.na(external_gene_name) != TRUE) %>% slice_min(order_by=log2FoldChange, n = 100), "EME+ZAKi_bottom_100.tsv")

#### Generate Broad GSEA input dataframe for Kate's code ---------

write_tsv(ANS_results$drug_vs_untreated %>% 
            select(ensembl_gene_id, log2FoldChange, pvalue, external_gene_name) %>% mutate(pvalue = pvalue + .Machine$double.xmin), "20251211_ANS-condition_GSEA_input.tsv")

#### Generate ZAK-independent genes in ANS+ZAKi vs EME+ZAKi vs DDB+ZAKi for transcription factor analysis for reviewers ---------

ANS_genes <- ANS_results$drug_vs_untreated %>% filter((padj < 0.01) & (log2FoldChange > 1)) %>% arrange(desc(stat))
EME_genes <- EME_results$drug_vs_untreated %>% filter((padj < 0.01) & (log2FoldChange > 1)) %>% arrange(desc(stat))
DDB_genes <- DDB_results$drug_vs_untreated %>% filter((padj < 0.01) & (log2FoldChange > 1)) %>% arrange(desc(stat))

write_tsv(ANS_genes, file = "script_outputs/TFenrichment/ANS_genes_lfc1.tsv")
write_tsv(EME_genes, file = "script_outputs/TFenrichment/EME_genes_lfc1.tsv")
write_tsv(DDB_genes, file = "script_outputs/TFenrichment/DDB_genes_lfc1.tsv")

ANSZAKi_genes <- ANS_results$drug_inhibitor_vs_untreated_inhibitor %>% filter((padj < 0.01) & (log2FoldChange > 1)) %>% arrange(desc(stat))
EMEZAKi_genes <- EME_results$drug_inhibitor_vs_untreated_inhibitor %>% filter((padj < 0.01) & (log2FoldChange > 1)) %>% arrange(desc(stat))
DDBZAKi_genes <- DDB_results$drug_inhibitor_vs_untreated_inhibitor %>% filter((padj < 0.01) & (log2FoldChange > 1)) %>% arrange(desc(stat))

write_tsv(ANSZAKi_genes, file = "script_outputs/TFenrichment/ANS+ZAKi_genes_lfc1.tsv")
write_tsv(EMEZAKi_genes, file = "script_outputs/TFenrichment/EME+ZAKi_genes_lfc1.tsv")
write_tsv(DDBZAKi_genes, file = "script_outputs/TFenrichment/DDB+ZAKi_genes_lfc1.tsv")

#### determine what is "ZAK-independent" by taking the DEGs and calculating the difference between noZAKi + _ZAKi first

ANS_genes <- ANS_genes %>% select(external_gene_name, log2FoldChange) %>% drop_na(external_gene_name)
ANSZAKi_genes <- ANSZAKi_genes %>% select(external_gene_name, log2FoldChange)%>% drop_na(external_gene_name)

combined <- full_join(ANS_genes, ANSZAKi_genes, by = "external_gene_name", suffix = c("_noZAKi", "_ZAKi")) %>% replace_na(list(log2FoldChange_noZAKi = 0, log2FoldChange_ZAKi = 0)) %>% mutate(difference = log2FoldChange_noZAKi - log2FoldChange_ZAKi)

combined_ZAK_indep <- combined %>% filter(difference < 1)
write_tsv(combined_ZAK_indep, file = "script_outputs/TFenrichment/ANS_ZAK_indep_genes.tsv")

combined_ZAK_dep <- combined %>% filter(difference >= 1)
write_tsv(combined_ZAK_dep, file = "script_outputs/TFenrichment/ANS_ZAK_dep_genes.tsv")

