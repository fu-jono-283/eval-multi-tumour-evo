library(ggtree)
library(cowplot)

# Read the tree from the file
tree <- read.tree("10397r2_G.trees")

# Group tip nodes based on the first number in their labels
grouped_tips <- split(tree$tip.label, substr(tree$tip.label, 1, 1))

# Convert the tree object to a ggtree object
p <- ggtree(tree)

# Add colored points to the tips of the tree based on grouping
p <- p + geom_tippoint(aes(color = factor(substr(label, 1, 1))), size = 3)

# Modify the tip labels to the desired format (e.g., "M7C8")
p <- p + geom_tiplab(aes(label = gsub("^(\\d+)_(\\d+)_(\\d+)$", "M\\1C\\3", label)))

# Remove the legend
p <- p + theme(legend.position = "none")+ geom_treescale()

# Plot the tree
print(p)
