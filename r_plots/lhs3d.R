library(plot3D)

# Create a 3D scatter plot
scatter3D(small_scores_27_$effective_population_size, 
          small_scores_27_$tree_wide_sub_rate, 
          small_scores_27_$ado, 
          colvar = small_scores_27_$no_of_taxa,  # Color by the "no_of_taxa" variable
          colkey = list(at = unique(small_scores_27_$no_of_taxa), 
                        col = heat.colors(length(unique(small_scores_27_$no_of_taxa))), 
                        length = 0.5, 
                        width = 0.5, 
                        cex.clab = 1),  # Set color key
          pch = 20,  # Set point character
          cex = 1,  # Set point size
          xlab = "Effective Population Size",
          ylab = "Mutation Rate",
          zlab = "ADO",
          main = "",
          col.axis = "black",  # Set axis color
          grid = TRUE,  # Show grid
          box = TRUE,  # Show box around plot
          phi = 5,  # Set the perspective angle
          theta = 45,  # Set the rotation angle
          bty = "g",  # Set box type
          ticktype ="detailed", 
          clab = "No. of Tumours")
