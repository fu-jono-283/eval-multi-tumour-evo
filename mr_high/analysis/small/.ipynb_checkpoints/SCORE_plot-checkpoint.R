library(ggplot2)
library(dplyr)
library(tidyr)

data <- small_scores
data$individuals_per_taxa <- as.numeric(data$individuals_per_taxa)
data <- na.omit(data)

data$scSEQ_score <- 1 - data$scSEQ_score
clean_data <- data[data$scSEQ_score != 2,]
data$hapcon_seq_score <- 1 - data$hapcon_seq_score 
clean_data <- data[data$hapcon_seq_score != 2,]
data$vcfbulk_seq_score <- 1 - data$vcfbulk_seq_score 
clean_data <- data[data$vcfbulk_seq_score != 2,]

long_data <- clean_data %>%
  pivot_longer(cols = c(scSEQ_score, hapcon_seq_score, vcfbulk_seq_score), names_to = "Method", values_to = "Score")

long_data$Method <- factor(long_data$Method, levels = c("scSEQ_score", "hapcon_seq_score", "vcfbulk_seq_score"))

plot <- ggplot(data = long_data, aes(x = Score, y = species_birth_rate, color = Method)) +
  stat_summary(
    geom = "point",
    fun.y = "mean",
    size = 0.1,
    shape = 1,
    fill = "chocolate1",
    alpha = 0.5  
  ) +
  labs(x = "Phylogenetic Accuracy Score", y = "Metastatic Rate",
       color = "Method of Reconstruction") +  
  scale_color_manual(values = c("orangered1", "chartreuse3", "deepskyblue"),
                     labels = c("scSEQ", "H-CS", "VCF-PB")) +  
  ggtitle("Phylogenetic Accuracy vs. Metastatic Rate: Comparison of Three Methods (MR-High)") +  
  theme_minimal() + 
  theme(
    text = element_text(family = "serif", color = "black"),  
    plot.title = element_text(size = 5, color = "black"),  
    axis.title = element_text(size = 5, color = "black"),  
    legend.text = element_text(size = 4.5, color = "black"),
    axis.text = element_text(size = 5, color = "black"),
    axis.line = element_line(colour = "black", size = 0.1),  
    axis.ticks = element_line(colour = "black", size = 0.1),  
    panel.grid = element_blank(),  
    legend.title = element_text(size = 5, hjust = -0.5),  
    legend.spacing.x = unit(0.1, "cm"),
    legend.key.size = unit(0.5, "lines")
  )

ggsave("12_high/results/small_scores/z_single_plot.png", plot, width = 3.5, height = 2.5, dpi = 320, bg = "white")
