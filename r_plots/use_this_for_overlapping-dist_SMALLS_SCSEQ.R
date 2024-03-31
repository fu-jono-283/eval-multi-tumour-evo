library(ggplot2)
library(reshape2)
library(dplyr)

data <- small_scores_high_SMALL
data$individuals_per_taxa <- as.numeric(data$individuals_per_taxa)

data$vcfbulk_seq_score <- 1 - data$vcfbulk_seq_score
data$hapcon_seq_score <-  1 - data$hapcon_seq_score
odensity_data <- data[data$vcfbulk_seq_score != 2, ]
odensity_data <- odensity_data[odensity_data$hapcon_seq_score != 2, ]

x <- data.frame(hapcon_seq_score = odensity_data$hapcon_seq_score,
                vcfbulk_seq_score = odensity_data$vcfbulk_seq_score)
# Create melted data for plotting
melted_data <- melt(x)

# Calculate summary statistics
summary_stats <- melted_data %>%
  mutate(Method = ifelse(variable == "hapcon_seq_score", "H-CS", "VCF-PB")) %>%
  group_by(Method) %>%
  summarize(
    Mean = sprintf("%.2f", mean(value)),
    SD = sprintf("%.2f", sd(value)),
    Median = sprintf("%.2f", median(value))
  )

# Custom font size for table
custom_table_font_size <- 15
custom_font_size <- 15
custom_legend_key_size <- 15

# Plot
gg <- ggplot(melted_data, aes(x = value, fill = variable)) +
  geom_density(alpha = 0.5, size = 1.5) +
  labs(
    x = "Phylogenetic Accuracy",
    y = "Density",
    fill = "Sequencing Method"
  ) +
  theme_minimal() +
  scale_fill_manual(values = c("chartreuse1", "deepskyblue"), labels = c("H-CS", "VCF-PB"), name = "Sequencing Method") +
  theme(
    plot.title = element_text(size = custom_font_size, color = "black", hjust = 0),
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

# Create the table
table_grob <- tableGrob(summary_stats, rows = NULL)

# Combine plot and table
library(gridExtra)
grid.arrange(gg, table_grob, nrow = 2, heights = c(3, 1))
