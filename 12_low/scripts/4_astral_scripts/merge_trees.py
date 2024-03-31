import dendropy
import os

root_directory = "data/4_astral"

for subdir, dirs, files in os.walk(root_directory):
    iteration_folder = os.path.basename(subdir)
    
    tree_files = [filename for filename in files if filename.startswith("r") and filename.endswith(".raxml.bestTree")]
    
    tree_files.sort()

    if tree_files:
        tree_list = dendropy.TreeList()
        
        for tree_file in tree_files:
            tree_path = os.path.join(subdir, tree_file)
            trees = dendropy.TreeList.get_from_path(tree_path, "newick")
            tree_list.extend(trees)

        output_file = os.path.join(subdir, "merged_trees.tre")
        tree_list.write_to_path(output_file, "newick")

        print("Merged trees for", iteration_folder, "written to:", output_file)
