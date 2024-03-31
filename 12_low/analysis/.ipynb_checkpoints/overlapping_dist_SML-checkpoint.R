library(ggplot2)
library(reshape2)

data <- small_scores
data$individuals_per_taxa <- as.numeric(data$individuals_per_taxa)

data$vcfbulk_seq_score <- 1 - data$vcfbulk_seq_score
data$hapcon_seq_score <-  1 - data$hapcon_seq_score
odensity_data <- data[data$vcfbulk_seq_score != 2, ]
odensity_data <- odensity_data[odensity_data$hapcon_seq_score != 2, ]

x <- data.frame(hapcon_seq_score = odensity_data$hapcon_seq_score,
                vcfbulk_seq_score = odensity_data$vcfbulk_seq_score)

melted_data <- melt(x)

custom_font_size <- 6  
custom_legend_key_size <- 8 

gg <- ggplot(melted_data, aes(x = value, fill = variable)) +
  geom_density(alpha = 0.5) +
  labs(
    x = "Phylogenetic Accuracy",
    y = "Density",
    title = "Phylogenetic Accuracy of Pooling Methods (SBR-Low): H-CS vs. VCF-PB",
    fill = "Sequencing Method"
  ) +
  theme_minimal() +
  scale_fill_manual(values = c("orangered1", "deepskyblue"),labels = c("H-CS", "VCF-PB"), name = "Sequencing Method") +
  theme(
    plot.title = element_text(size = custom_font_size, family = "serif", color = "black", hjust = .75),
    axis.title = element_text(size = custom_font_size, family = "serif", color = "black"),
    legend.text = element_text(size = custom_font_size, family = "serif", color = "black"),
    legend.title = element_text(size = custom_font_size, family = "serif", color = "black"),
    legend.key.size = unit(custom_legend_key_size, "pt"),
    axis.text = element_text(size = custom_font_size, family = "serif", color = "black"),
    axis.line = element_line(colour = "black", size = 0.1),
    axis.ticks = element_line(colour = "black", size = 0.1),
    panel.grid = element_blank(),
    plot.margin = margin(10, 10, 10, 10),
  )

save_path <- "12_low/results/small_scores/z_overlapping_dist_vcf_hap.png"
ggsave(filename = save_path, device = "png", bg = "white", width = 5.25, height = 3, dpi = 320)
