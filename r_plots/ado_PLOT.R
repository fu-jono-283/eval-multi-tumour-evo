# Plot for MR-Full with ADO
plot_fullADO <- ggplot(data = scSEQ_data_full, aes(x = ado, y = Score)) +
  geom_point(size = 2, alpha = 0.2, color = "cyan4") +
  labs(x = "ADO", y = "Phylogenetic Accuracy") +  
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
  ylim(0.4, 1) +
  xlim(0, 0.2)  # Set y-axis to ADO range

# Plot for MR-Low with ADO
plot_lowADO <- ggplot(data = scSEQ_data_low, aes(x = ado, y = Score)) +
  geom_point(size = 2, alpha = 0.2, color = "cyan4") +
  labs(x = "ADO", y = "Phylogenetic Accuracy") + 
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
  ylim(0.4, 1)+ # Set x-axis
  xlim(0, 0.2)  # Set y-axis to ADO range

# Plot for MR-High with ADO
plot_highADO <- ggplot(data = scSEQ_data_high, aes(x = ado, y = Score)) +
  geom_point(size = 2, alpha = 0.2, color = "cyan4") +
  labs(x = "ADO", y = "Phylogenetic Accuracy") + 
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
  ylim(0.4, 1) +
  xlim(0, 0.2)  # Set y-axis to ADO range


plot((plot_fullADO))
plot((plot_highADO))
plot((plot_lowADO))


# Save plot for MR-Full
ggsave("plot_fullADO.pdf", plot_fullADO, width = 1600/300, height = 1200/300, units = "in")

# Save plot for MR-Low
ggsave("plot_lowADO.pdf", plot_lowADO, width = 1600/300, height = 1200/300, units = "in")

# Save plot for MR-High
ggsave("plot_highADO.pdf", plot_highADO, width = 1600/300, height = 1200/300, units = "in")
