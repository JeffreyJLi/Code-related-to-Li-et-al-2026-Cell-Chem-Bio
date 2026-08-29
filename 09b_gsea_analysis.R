library("here")
library("tidyverse")
library("enrichplot")
library("ggplot2")
library("ggrepel")
library("msigdbr")

# Define default themes for plotting
theme_greenlab <- function() {
  theme(axis.text = element_text(size = 14), axis.title = element_text(size = 20, face="bold"), plot.title = element_text(size = 24, face = "bold"),
        plot.subtitle = element_text(size = 18), legend.title = element_text(size = 14, face = "bold"), legend.text = element_text(size = 14), legend.key.spacing.y = unit(1, "mm"))
}
greenlab_colorpalette <- c("#ff0c68", "#008E91","#981cd7")

#### HALLMARK GSEA:

ANS_gsea <- read_tsv(file = here("script_outputs/reviewer_analysis/gsea/Hallmark_sets/ANS_gsea_results_hallmark.tsv"))
EME_gsea <- read_tsv(file = here("script_outputs/reviewer_analysis/gsea/Hallmark_sets/EME_gsea_results_hallmark.tsv"))
DDB_gsea <- read_tsv(file = here("script_outputs/reviewer_analysis/gsea/Hallmark_sets/DDB_gsea_results_hallmark.tsv"))
ANSZAKi_gsea <- read_tsv(file = here("script_outputs/reviewer_analysis/gsea/Hallmark_sets/ANSZAKi_gsea_results_hallmark.tsv"))
EMEZAKi_gsea <- read_tsv(file = here("script_outputs/reviewer_analysis/gsea/Hallmark_sets/EMEZAKi_gsea_results_hallmark.tsv"))
DDBZAKi_gsea <- read_tsv(file = here("script_outputs/reviewer_analysis/gsea/Hallmark_sets/DDBZAKi_gsea_results_hallmark.tsv"))

##### Generate summary NES table of all conditions --------------------------------------------------------------------------------------------------------
summarize_gsea_NES <- function(gsea_result, drug) {
  summarized <- gsea_result %>% select(ID, NES) %>% rename_with(~ paste0(., "_", drug), .cols ="NES")
  return(summarized)
}

ANS_hallmark_NES <- summarize_gsea_NES(ANS_gsea, "ANS")
EME_hallmark_NES <- summarize_gsea_NES(EME_gsea, "EME")
DDB_hallmark_NES <- summarize_gsea_NES(DDB_gsea, "DDB")
ANSZAKi_hallmark_NES <- summarize_gsea_NES(ANSZAKi_gsea, "ANSZAKi")
EMEZAKi_hallmark_NES <- summarize_gsea_NES(EMEZAKi_gsea, "EMEZAKi")
DDBZAKi_hallmark_NES <- summarize_gsea_NES(DDBZAKi_gsea, "DDBZAKi")

hallmark_NES <- list(ANS_hallmark_NES, EME_hallmark_NES, DDB_hallmark_NES, ANSZAKi_hallmark_NES, EMEZAKi_hallmark_NES,DDBZAKi_hallmark_NES)

combined_hallmark_NES <- hallmark_NES %>% reduce(full_join, by = "ID")

if (!file.exists(here("script_outputs/reviewer_analysis/gsea/Hallmark_sets/combined_hallmark_NES_table.tsv"))) {
  write_tsv(combined_hallmark_NES, file = "script_outputs/reviewer_analysis/gsea/Hallmark_sets/combined_hallmark_NES_table.tsv")
} else {
  print("File already exists.")
}

##### Plot NES scores of each conditions --------------------------------------------------------------------------------------------------------
# Clean up data for plotting
ANS_hallmark_NES <- ANS_gsea %>% select(ID, NES, p.adjust) %>% mutate(drug = "ANS")
EME_hallmark_NES <- EME_gsea %>% select(ID, NES, p.adjust) %>% mutate(drug = "EME")
DDB_hallmark_NES <- DDB_gsea %>% select(ID, NES, p.adjust) %>% mutate(drug = "DDB")
ANSZAKi_hallmark_NES <- ANSZAKi_gsea %>% select(ID, NES, p.adjust) %>% mutate(drug = "ANSZAKi")
EMEZAKi_hallmark_NES <- EMEZAKi_gsea %>% select(ID, NES, p.adjust) %>% mutate(drug = "EMEZAKi")
DDBZAKi_hallmark_NES <- DDBZAKi_gsea %>% select(ID, NES,p.adjust) %>% mutate(drug = "DDBZAKi")

# Define helper functions for readability
reorder_df_for_plotting <- function(df) {
  # Define rank order of IDs for plotting later
  ord <- df %>% 
    mutate(ID = str_remove(ID, "HALLMARK_")) %>% 
    group_by(ID) %>% 
    summarise(m = max(NES)) %>% 
    arrange(m) %>% 
    pull(ID)
  
  # Generate plotting df, filling in missing values for better plotting.
  # Otherwise, some rows with missing values will have weird bar ordering...
  plot_df <- df %>% 
    mutate(ID = str_remove(ID, "HALLMARK_")) %>% 
    mutate(drug = factor(drug)) %>% 
    complete(ID, drug, fill = list(NES = 0, p.adjust = 1), explicit = FALSE) %>% 
    mutate(ID = factor(ID, levels = ord))
  return(plot_df)
}

save_plot <- function(plot, filename){
  if (!file.exists(here("script_outputs/reviewer_analysis/gsea/plots/", filename))) {
    ggsave(plot = plot, filename = here("script_outputs/reviewer_analysis/gsea/plots/", filename), width = 14, height = 12.3, units = "in")
  } else {
    print("File already exists.")
  }
}

#### 
Add the ANS+EME and -/+ ZAKi GSEA plots and write short paragraph citing the supplemental. EME shuts off metabolism; ANS hits UV
ADD top 10 CHEA3 analysis to rebuttal. Lots of the same transcription factors. GSEA is more helpful to think about these pathways.
ggsave(here("script_outputs/reviewer_analysis/gsea/plots/ANS-ANS+ZAKi_GSEA.pdf"), plot = plot_ANS_ANSZAKi_combinedNES, width = 5, height = 8)

# Plot ANS vs ANSZAKi
ANS_ANSZAKi_combinedNES <- bind_rows(ANS_hallmark_NES, ANSZAKi_hallmark_NES)

plot_ANS_ANSZAKi_combinedNES <- ggplot(reorder_df_for_plotting(ANS_ANSZAKi_combinedNES), aes(x = NES, y = ID, color = drug, fill = drug, size = -log10(p.adjust))) +
  geom_point(shape = 21, stroke = 1)+
  labs(x = "Normalized Enrichment Score (NES)", y = "GSEA Hallmark Gene Set", 
       title = "GSEA: ANS and ANS+ZAKi") +
  geom_vline(xintercept = 0, color = "black", linewidth = 0.6) +
  theme_bw() +
  scale_color_manual(values = c("ANS" = greenlab_colorpalette[1], "ANSZAKi" = greenlab_colorpalette[1]),
                     labels = c("ANSZAKi" = "ANS+ZAKi", "EMEZAKi" = "EME+ZAKi", "DDBZAKi" = "DDB+ZAKi"))+
  scale_fill_manual(
    values = c("ANS" = greenlab_colorpalette[1], "ANSZAKi" = NA),
    labels = c("ANSZAKi" = "ANS+ZAKi", "EMEZAKi" = "EME+ZAKi", "DDBZAKi" = "DDB+ZAKi")) +
  theme_greenlab() +
  scale_size(range = c(3, 8))

plot_ANS_ANSZAKi_combinedNES
save_plot(plot_ANS_ANSZAKi_combinedNES, "ANSvsANS+ZAKi_GSEA.pdf")

# Plot EME vs EMEZAKi
EME_EMEZAKi_combinedNES <- bind_rows(EME_hallmark_NES, EMEZAKi_hallmark_NES)

plot_EME_EMEZAKi_combinedNES <- ggplot(reorder_df_for_plotting(EME_EMEZAKi_combinedNES), aes(x = NES, y = ID, color = drug, fill = drug, size = -log10(p.adjust))) +
  geom_point(shape = 21, stroke = 1)+
  labs(x = "Normalized Enrichment Score (NES)", y = "GSEA Hallmark Gene Set", 
       title = "GSEA: EME and EME+ZAKi") +
  geom_vline(xintercept = 0, color = "black", linewidth = 0.6) +
  theme_bw() +
  scale_color_manual(values = c("EME" = greenlab_colorpalette[2], "EMEZAKi" = greenlab_colorpalette[2]),
                     labels = c("ANSZAKi" = "ANS+ZAKi", "EMEZAKi" = "EME+ZAKi", "DDBZAKi" = "DDB+ZAKi"))+
  scale_fill_manual(
    values = c("EME" = greenlab_colorpalette[2], "EMEZAKi" = NA),
    labels = c("ANSZAKi" = "ANS+ZAKi", "EMEZAKi" = "EME+ZAKi", "DDBZAKi" = "DDB+ZAKi")) +
  theme_greenlab() +
  scale_size(range = c(3, 8))
plot_EME_EMEZAKi_combinedNES
save_plot(plot_EME_EMEZAKi_combinedNES, "EMEvsEME+ZAKi_GSEA.pdf")


# Plot DDB vs DDBZAKi
DDB_DDBZAKi_combinedNES <- bind_rows(DDB_hallmark_NES, DDBZAKi_hallmark_NES)

plot_DDB_DDBZAKi_combinedNES <- ggplot(reorder_df_for_plotting(DDB_DDBZAKi_combinedNES), aes(x = NES, y = ID, color = drug, fill = drug, size = -log10(p.adjust))) +
  geom_point(shape = 21, stroke = 1)+
  labs(x = "Normalized Enrichment Score (NES)", y = "GSEA Hallmark Gene Set", 
       title = "GSEA: DDB and DDB+ZAKi") +
  geom_vline(xintercept = 0, color = "black", linewidth = 0.6) +
  theme_bw() +
  scale_color_manual(values = c("DDB" = greenlab_colorpalette[3], "DDBZAKi" = greenlab_colorpalette[3]),
                     labels = c("ANSZAKi" = "ANS+ZAKi", "EMEZAKi" = "EME+ZAKi", "DDBZAKi" = "DDB+ZAKi"))+
  scale_fill_manual(
    values = c("DDB" = greenlab_colorpalette[3], "DDBZAKi" = NA),
    labels = c("ANSZAKi" = "ANS+ZAKi", "EMEZAKi" = "EME+ZAKi", "DDBZAKi" = "DDB+ZAKi")) +
  theme_greenlab() +
  scale_size(range = c(3, 8))
plot_DDB_DDBZAKi_combinedNES
save_plot(plot_DDB_DDBZAKi_combinedNES, "DDBvsDDB+ZAKi_GSEA.pdf")

# Plot all drugs -ZAKi on same plot
alldrug_combinedNES <- bind_rows(ANS_hallmark_NES, EME_hallmark_NES, DDB_hallmark_NES)
alldrug_combinedNES <- alldrug_combinedNES %>%  mutate(drug = factor(drug, levels = c("ANS", "EME", "DDB")))

plot_alldrug_combinedNES <- ggplot(reorder_df_for_plotting(alldrug_combinedNES), aes(x = NES, y = ID, color = drug, fill = drug, size = -log10(p.adjust))) +
  geom_point(shape = 21, stroke = 1, alpha = 0.5)+
  labs(x = "Normalized Enrichment Score (NES)", y = "GSEA Hallmark Gene Set", 
       title = "GSEA: ANS/EME/DDB -ZAKi") +
  geom_vline(xintercept = 0, color = "black", linewidth = 0.6) +
  theme_bw() +
  scale_color_manual(values = c("ANS" = greenlab_colorpalette[1], "EME" = greenlab_colorpalette[2], "DDB" = greenlab_colorpalette[3]),
                     labels = c("ANSZAKi" = "ANS+ZAKi", "EMEZAKi" = "EME+ZAKi", "DDBZAKi" = "DDB+ZAKi"))+
  scale_fill_manual(
    values = c("ANS" = greenlab_colorpalette[1], "EME" = greenlab_colorpalette[2],"DDB" = greenlab_colorpalette[3]),
    labels = c("ANSZAKi" = "ANS+ZAKi", "EMEZAKi" = "EME+ZAKi", "DDBZAKi" = "DDB+ZAKi")) +
  theme_greenlab() +
  scale_size(range = c(3, 8))
plot_alldrug_combinedNES
save_plot(plot_alldrug_combinedNES, "alldrug-ZAKi_GSEA.pdf")

# Plot all drugs +ZAKi on same plot
alldrugZAKi_combinedNES <- bind_rows(ANSZAKi_hallmark_NES, EMEZAKi_hallmark_NES, DDBZAKi_hallmark_NES)
alldrugZAKi_combinedNES <- alldrugZAKi_combinedNES %>%  mutate(drug = factor(drug, levels = c("ANSZAKi", "EMEZAKi", "DDBZAKi")))

plot_alldrugZAKi_combinedNES <- ggplot(reorder_df_for_plotting(alldrugZAKi_combinedNES), aes(x = NES, y = ID, color = drug, fill = drug, size = -log10(p.adjust))) +
  geom_point(shape = 21, stroke = 1, alpha = 0.5)+
  labs(x = "Normalized Enrichment Score (NES)", y = "GSEA Hallmark Gene Set", 
       title = "GSEA: ANS/EME/DDB +ZAKi") +
  geom_vline(xintercept = 0, color = "black", linewidth = 0.6) +
  theme_bw() +
  scale_color_manual(values = c("ANSZAKi" = greenlab_colorpalette[1], "EMEZAKi" = greenlab_colorpalette[2], "DDBZAKi" = greenlab_colorpalette[3]),
                     labels = c("ANSZAKi" = "ANS+ZAKi", "EMEZAKi" = "EME+ZAKi", "DDBZAKi" = "DDB+ZAKi"))+
  scale_fill_manual(
    values = c("ANSZAKi" = greenlab_colorpalette[1], "EMEZAKi" = greenlab_colorpalette[2],"DDBZAKi" = greenlab_colorpalette[3]),
    labels = c("ANSZAKi" = "ANS+ZAKi", "EMEZAKi" = "EME+ZAKi", "DDBZAKi" = "DDB+ZAKi")) +
  theme_greenlab() +
  scale_size(range = c(3, 8))
plot_alldrugZAKi_combinedNES
save_plot(plot_alldrugZAKi_combinedNES, "alldrug+ZAKi_GSEA.pdf")



h <- msigdbr(species = "Homo sapiens", collection = "H")

tnfa_genes <- h |>
  filter(gs_name == "HALLMARK_UV_RESPONSE_UP")

ANSZAKi_gsea$ID

ANS_gsea %>% filter(ID == "HALLMARK_UV_RESPONSE_UP") %>% select(core_enrichment)
ANS_gsea %>% filter(ID == "HALLMARK_ESTROGEN_RESPONSE_LATE") %>% select(core_enrichment)
ANSZAKi_gsea %>% filter(ID == "HALLMARK_OXIDATIVE_PHOSPHORYLATION") %>% select(core_enrichment) %>% toString()

##### --------------------------------------------------------------------------------------------------------
#### C5 GSEA:

ANS_gsea <- read_tsv(file = here("C5_sets//ANS_gsea_results_C5collection.tsv"))
EME_gsea <- read_tsv(file = here("C5_sets/EME_gsea_results_C5collection.tsv"))
DDB_gsea <- read_tsv(file = here("C5_sets/DDB_gsea_results_C5collection.tsv"))
ANSZAKi_gsea <- read_tsv(file = here("C5_sets/ANSZAKi_gsea_results_C5collection.tsv"))
EMEZAKi_gsea <- read_tsv(file = here("C5_sets/EMEZAKi_gsea_results_C5collection.tsv"))
DDBZAKi_gsea <- read_tsv(file = here("C5_sets/DDBZAKi_gsea_results_C5collection.tsv"))

summarize_gsea_NES <- function(gsea_result, drug) {
  summarized <- gsea_result %>% select(ID, NES) %>% rename_with(~ paste0(., "_", drug), .cols ="NES")
  return(summarized)
}

ANS_C5_NES <- summarize_gsea_NES(ANS_gsea, "ANS")
EME_C5_NES <- summarize_gsea_NES(EME_gsea, "EME")
DDB_C5_NES <- summarize_gsea_NES(DDB_gsea, "DDB")
ANSZAKi_C5_NES <- summarize_gsea_NES(ANSZAKi_gsea, "ANSZAKi")
EMEZAKi_C5_NES <- summarize_gsea_NES(EMEZAKi_gsea, "EMEZAKi")
DDBZAKi_C5_NES <- summarize_gsea_NES(DDBZAKi_gsea, "DDBZAKi")

C5_NES <- list(ANS_C5_NES, EME_C5_NES, DDB_C5_NES, ANSZAKi_C5_NES, EMEZAKi_C5_NES,DDBZAKi_C5_NES)

combined_C5_NES <- C5_NES %>% reduce(full_join, by = "ID")

if (!file.exists(here("C5_sets/combined_C5_NES_table.tsv"))) {
  write_tsv(combined_hallmark_NES, file = "C5_sets/combined_C5_NES_table.tsv")
} else {
  print("File already exists.")
}


