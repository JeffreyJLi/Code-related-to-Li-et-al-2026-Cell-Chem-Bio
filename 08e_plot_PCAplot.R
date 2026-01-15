plot_PCAplot <- function(PCAplot_data, percentVariance) {
  PCAplot <- ggplot(PCAplot_data, aes(x = PC1, y = PC2, color = drug, shape = inhibitor)) +
    geom_point(size = 4, stroke = 2) +
    theme_bw(base_family = "Avenir") +
    theme(legend.title = element_text(size = 20, face="bold"),
          legend.text = element_text(size = 18),
          axis.text = element_text(size = 18),
          axis.title = element_text(size = 20, face="bold"),
          plot.title = element_text(size = 20, face = "bold"),
          panel.grid = element_line(color = "dark grey"),
          legend.key.spacing.y = unit(0.5, "lines")) + 
    coord_fixed() +
    ggtitle("RNA-seq PCA") +
    xlab(paste0("PC1: ",percentVariance[1],"% variance")) +
    ylab(paste0("PC2: ",percentVariance[2],"% variance")) +
    scale_color_manual(name = "group", 
                       labels = c("Untreated","ANS", "EME", "DDB"),
                       values = c("#afafaf", "#b100ff","#084152","#62d5b4")) +
    scale_shape_manual(name = "inhibitor", 
                       labels = c("no ZAKi","+ZAKi"),
                       values = c(19,1))
  return(PCAplot)
}

