library(ggplot2)
library(dplyr)
library(tidyr)
library(gridExtra)

# Load data and preprocess for MR-Full
data <- small_scores_27_  # Use small_scores_27_ dataset for MR-Full
data$individuals_per_taxa <- as.numeric(data$individuals_per_taxa)
data <- na.omit(data)

data$scSEQ_score <- 1 - data$scSEQ_score
clean_data <- data[data$scSEQ_score != 2,]
data$hapcon_seq_score <- 1 - data$hapcon_seq_score 
clean_data <- data[data$hapcon_seq_score != 2,]
data$vcfbulk_seq_score <- 1 - data$vcfbulk_seq_score 
clean_data <- data[data$vcfbulk_seq_score != 2,]

# Pivot longer
long_data <- clean_data %>%
  pivot_longer(cols = c(scSEQ_score, hapcon_seq_score, vcfbulk_seq_score), names_to = "Method", values_to = "Score")

long_data$Method <- factor(long_data$Method, levels = c("scSEQ_score"))

# Filter the data for only scSEQ method for MR-Full
scSEQ_data_full <- filter(long_data, Method == "scSEQ_score")

# Plot for MR-Full
plot_full <- ggplot(data = scSEQ_data_full, aes(x = Score, y = species_birth_rate)) +
  geom_point(size = 2, alpha = 0.2, color = "orangered1") +
  labs(x = "Phylogenetic Accuracy", y = "Metastatic Rate") +  
  ggtitle("Performance of scSEQ (MR-Full)") +  
  theme_minimal() + 
  theme(
    plot.title = element_text(size = 14, color = "black", hjust = 0.5), # Center the title
    axis.title = element_text(size = 14, color = "black"),  
    legend.text = element_text(size = 14, color = "black"),
    axis.text = element_text(size = 14, color = "black"),
    axis.line = element_line(colour = "black", size = 0.1),  
    axis.ticks = element_line(colour = "black", size = 0.1),  
    panel.grid = element_blank(),  
    legend.title = element_text(size = 14, hjust = -0.5),  
    legend.spacing.x = unit(0.1, "cm"),
    legend.key.size = unit(0.5, "lines")
  ) +
  guides(color = FALSE)+  # Remove the legend
  xlim(0.4, 1) +
  ylim(0, 0.0001)  # Set x-axis

# Load data and preprocess for MR-Low
data <- small_scores_28_  # Use small_scores_28_ dataset for MR-Low
data$individuals_per_taxa <- as.numeric(data$individuals_per_taxa)
data <- na.omit(data)

data$scSEQ_score <- 1 - data$scSEQ_score
clean_data <- data[data$scSEQ_score != 2,]
data$hapcon_seq_score <- 1 - data$hapcon_seq_score 
clean_data <- data[data$hapcon_seq_score != 2,]
data$vcfbulk_seq_score <- 1 - data$vcfbulk_seq_score 
clean_data <- data[data$vcfbulk_seq_score != 2,]

# Pivot longer
long_data <- clean_data %>%
  pivot_longer(cols = c(scSEQ_score, hapcon_seq_score, vcfbulk_seq_score), names_to = "Method", values_to = "Score")

long_data$Method <- factor(long_data$Method, levels = c("scSEQ_score"))

# Filter the data for only scSEQ method for MR-Low
scSEQ_data_low <- filter(long_data, Method == "scSEQ_score")

# Plot for MR-Low
plot_low <- ggplot(data = scSEQ_data_low, aes(x = Score, y = species_birth_rate)) +
  geom_point(size = 2, alpha = 0.2, color = "orangered1") +
  labs(x = "Phylogenetic Accuracy", y = "Metastatic Rate") +  
  ggtitle("Performance of scSEQ (MR-Low)") +  
  theme_minimal() + 
  theme(
    plot.title = element_text(size = 14, color = "black", hjust = 0.5), # Center the title
    axis.title = element_text(size = 14, color = "black"),  
    legend.text = element_text(size = 14, color = "black"),
    axis.text = element_text(size = 14, color = "black"),
    axis.line = element_line(colour = "black", size = 0.1),  
    axis.ticks = element_line(colour = "black", size = 0.1),  
    panel.grid = element_blank(),  
    legend.title = element_text(size = 14, hjust = -0.5),  
    legend.spacing.x = unit(0.1, "cm"),
    legend.key.size = unit(0.5, "lines")
  ) +
  guides(color = FALSE) +  # Remove the legend
  xlim(0.4, 1)+ # Set x-axis+
  ylim(0, 0.0001)
  
data <- small_scores_29_
data$individuals_per_taxa <- as.numeric(data$individuals_per_taxa)
data <- na.omit(data)

data$scSEQ_score <- 1 - data$scSEQ_score
clean_data <- data[data$scSEQ_score != 2,]
data$hapcon_seq_score <- 1 - data$hapcon_seq_score 
clean_data <- data[data$hapcon_seq_score != 2,]
data$vcfbulk_seq_score <- 1 - data$vcfbulk_seq_score 
clean_data <- data[data$vcfbulk_seq_score != 2,]

# Pivot longer
long_data <- clean_data %>%
  pivot_longer(cols = c(scSEQ_score, hapcon_seq_score, vcfbulk_seq_score), names_to = "Method", values_to = "Score")

long_data$Method <- factor(long_data$Method, levels = c("scSEQ_score"))

# Filter the data for only scSEQ method for MR-High
scSEQ_data_high <- filter(long_data, Method == "scSEQ_score")

# Plot for MR-High
plot_high <- ggplot(data = scSEQ_data_high, aes(x = Score, y = species_birth_rate)) +
  geom_point(size = 2, alpha = 0.2, color = "orangered1") +
  labs(x = "Phylogenetic Accuracy", y = "Metastatic Rate") +  
  ggtitle("Performance of scSEQ (MR-High)") +  
  theme_minimal() + 
  theme(
    plot.title = element_text(size = 14, color = "black", hjust = 0.5), # Center the title
    axis.title = element_text(size = 14, color = "black"),  
    legend.text = element_text(size = 14, color = "black"),
    axis.text = element_text(size = 14, color = "black"),
    axis.line = element_line(colour = "black", size = 0.1),  
    axis.ticks = element_line(colour = "black", size = 0.1),  
    panel.grid = element_blank(),  
    legend.title = element_text(size = 14, hjust = -0.5),  
    legend.spacing.x = unit(0.1, "cm"),
    legend.key.size = unit(0.5, "lines")
  ) +
  guides(color = FALSE) +  # Remove the legend
  xlim(0.4, 1) +
  ylim(0, 0.0001)# Set x-axis

plot(plot_full)
plot(plot_low)

# Save plot for MR-Full
ggsave("plot_full.pdf", plot_full, width = 1600/300, height = 1200/300, units = "in")

# Save plot for MR-Low
ggsave("plot_low.pdf", plot_low, width = 1600/300, height = 1200/300, units = "in")

# Save plot for MR-High
ggsave("plot_high.pdf", plot_high, width = 1600/300, height = 1200/300, units = "in")

