library(ggtree)
library(cowplot)

# Read the tree from the file
tree <- read.tree("10397_S.trees")

# Group tip nodes based on the first number in their labels
grouped_tips <- split(tree$tip.label, substr(tree$tip.label, 1, 1))

# Convert the tree object to a ggtree object
z <- ggtree(tree)

# Add colored points to the tips of the tree based on grouping
z <- z + geom_tippoint(aes(color = factor(substr(label, 1, 1))), size = 3)

# Modify the tip labels to add "M" in front of the numerical value
z <- z + geom_tiplab(aes(label = paste0("M", label)))

# Remove the legend
z <- z + theme(legend.position = "none")+ geom_treescale()

# Plot the tree
print(z)
