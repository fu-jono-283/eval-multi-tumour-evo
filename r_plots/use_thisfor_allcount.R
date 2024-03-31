  library(dplyr)
  library(ggplot2)
  
  allcount <- allcountdata_3_
  
  # Calculate Proportion_N
  allcount$total_sites <- allcount$no_of_taxa * allcount$individuals_per_taxa * allcount$no_of_sites
  
  allcount$Proportion_ALL <- allcount$ALL_count / allcount$total_sites
  
  # Calculate mean for each interval
  mean_data <- allcount %>%
    group_by(CutoffValue) %>%
    summarize(mean_Proportion_ALL = mean(Proportion_ALL, na.rm = TRUE))
  
  # Create the plot
  all_plot <- ggplot(data = allcount, aes(x = CutoffValue, y = Proportion_ALL)) +
    geom_point(color = "green3", alpha = 0.05) +
    geom_point(data = mean_data, aes(y = mean_Proportion_ALL), color = "orange", size = 3) +
    geom_text(data = mean_data, aes(label = round(mean_Proportion_ALL, 2), y = mean_Proportion_ALL + 0.01), color = "black", vjust = -0.5, size = 3.5) +
    labs(x = "Cutoff Value", y = "Proportion of IUPAC Ambiguity Bases") +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 10),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_line(color = "black"),
      axis.ticks = element_line(colour = "black", size = 0.1),
      axis.text = element_text(size = 15),
      axis.title = element_text(size = 15)
    )
  
  all_plot
  
  # Save the plot as a PDF
  ggsave("all_count.pdf", my_plot, width = 1600/300, height = 1200/300, units = "in")
