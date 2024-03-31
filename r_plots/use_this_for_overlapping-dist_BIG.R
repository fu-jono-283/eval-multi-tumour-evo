library(ggplot2)
library(reshape2)
library(gridExtra)
library(dplyr)
library(grid)

data <- big_scores_full_BIG
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

summary_stats <- melted_data %>%
  group_by(variable) %>%
  summarize(
    Mean = sprintf("%.2f", mean(value)),
    SD = sprintf("%.2f", sd(value)),
    Median = sprintf("%.2f", median(value))
  ) %>%
  mutate(
    Method = case_when(
      variable == "ASTRAL_score" ~ "scSEQ-ASTRAL",
      variable == "concat_seq_score" ~ "H-CS-CONCAT",
      variable == "hap_FastRFS_score" ~ "H-CS-FastRFS",
      variable == "vcf_merge_seq_score" ~ "VCF-PB-CONCAT",
      variable == "vcf_FastRFS_score" ~ "VCF-PB-FASTRFS",
      TRUE ~ "Other"
    )
  ) %>%
  select(Method, Mean, SD, Median)

# Custom font size for table
custom_table_font_size <- 15
custom_font_size <- 15
custom_legend_key_size <- 15

gg <- ggplot(melted_data, aes(x = value, fill = variable)) +
  geom_density(alpha = 0.5, size = 1.5,  position = "identity") +
  labs(
    x = "Phylogenetic Accuracy",
    y = "Density",
    fill = "Sequencing Method"
  ) +
  theme_minimal() +
  scale_fill_discrete(
    name = "Method",
    labels = c("scSEQ-ASTRAL", "H-CS-CONCAT", "H-CS-FastRFS", "VCF-PB-CONCAT", "VCF-PB-FastRFS")
  ) +
  scale_x_continuous(breaks = seq(0, 1, 0.25)) +  # Adjust x-axis label intervals
  theme(
    plot.title = element_text(size = custom_font_size, color = "black"),
    axis.title = element_text(size = custom_font_size, color = "black"),
    legend.text = element_text(size = custom_font_size, color = "black"),
    legend.title = element_text(size = custom_font_size, color = "black"),
    legend.key.size = unit(custom_legend_key_size, "pt"),
    axis.text = element_text(size = custom_font_size, color = "black"),
    axis.line = element_line(colour = "black", size = 0.1),
    axis.ticks = element_line(colour = "black", size = 0.1),
    panel.grid = element_blank(),
    plot.margin = margin(10, 10, 10, 10)  
  )
gg

# Create the table
table_grob <- tableGrob(summary_stats, rows = NULL)

# Combine plot and table using grid.arrange
grid.arrange(ggplotGrob(gg), table_grob, nrow = 2, heights = c(3, 1))
