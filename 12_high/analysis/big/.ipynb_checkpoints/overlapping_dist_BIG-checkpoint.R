library(ggplot2)
library(reshape2)

data <- big_scores
data$individuals_per_taxa <- as.numeric(data$individuals_per_taxa)
data <- na.omit(data)

data$ASTRAL_score <- 1 - data$ASTRAL_score
data$concat_seq_score <- 1 - data$concat_seq_score
data$hap_FastRFS_score <- 1 - data$hap_FastRFS_score
data$vcf_FastRFS_score <- 1 - data$vcf_FastRFS_score
data$vcf_merge_seq_score <- 1 - data$vcf_merge_seq_score

odensity_data <- data[data$ASTRAL_score != 2, ]
odensity_data <- odensity_data[odensity_data$concat_seq_score != 2, ]
odensity_data <- odensity_data[odensity_data$hap_FastRFS_score != 2, ]
odensity_data <- odensity_data[odensity_data$vcf_FastRFS_score != 2, ]
odensity_data <- odensity_data[odensity_data$vcf_merge_seq_score != 2, ]

x <- data.frame(
  ASTRAL_score = odensity_data$ASTRAL_score,
  concat_seq_score = odensity_data$concat_seq_score,
  hap_FastRFS_score = odensity_data$hap_FastRFS_score,
  vcf_merge_seq_score = odensity_data$vcf_merge_seq_score,
  vcf_FastRFS_score = odensity_data$vcf_FastRFS_score
)

melted_data <- melt(x)

custom_font_size <- 6
custom_legend_key_size <- 8

gg <- ggplot(melted_data, aes(x = value, fill = variable)) +
  geom_density(alpha = 0.4, position = "identity") +
  labs(
    x = "Phylogenetic Accuracy",
    y = "Density",
    title = "Distribution of Phylogenetic Accuracy of Various Methods (MR-High)",
    fill = "Sequencing Method"
  ) +
  theme_minimal() +
  scale_fill_discrete(
    name = "Sequencing Method",
    labels = c("ASTRAL", "Concat-H-CS", "FastRFS-H-CS", "Concat-VCF-PB", "FastRFS-VCF-PB")
  ) +
  theme(
    plot.title = element_text(size = custom_font_size, family = "serif", color = "black"),
    axis.title = element_text(size = custom_font_size, family = "serif", color = "black"),
    legend.text = element_text(size = custom_font_size, family = "serif", color = "black"),
    legend.title = element_text(size = custom_font_size, family = "serif", color = "black"),
    legend.key.size = unit(custom_legend_key_size, "pt"),
    axis.text = element_text(size = custom_font_size, family = "serif", color = "black"),
    axis.line = element_line(colour = "black", size = 0.1),
    axis.ticks = element_line(colour = "black", size = 0.1),
    panel.grid = element_blank(),
    plot.margin = margin(10, 10, 10, 10)  
  )

save_path <- "12_high/results/big_scores/z_overlapping_dist_methods.png"
ggsave(filename = save_path, device = "png", bg = "white", width = 6, height = 2.5, dpi = 320)
