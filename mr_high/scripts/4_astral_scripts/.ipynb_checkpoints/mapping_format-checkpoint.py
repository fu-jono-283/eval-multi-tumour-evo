import dendropy
import os

root_directory = "data/4_astral"

for subdir, dirs, files in os.walk(root_directory):
    iteration_folder = os.path.basename(subdir)
    tree_file_path = os.path.join(subdir, "merged_trees.tre")

    if os.path.isfile(tree_file_path):
        mapping = {}

        with open(tree_file_path, 'r') as tree_file:
            tree_string = tree_file.read()

        tree_cells = tree_string.split(',')

        for cell in tree_cells:
            label = cell.split(':')[0].replace('(', '').replace(')', '').split('_')
            species_label = label[0]
            individual_identifier = "_".join(label[1:])  
            individual_key = f"{species_label}_{individual_identifier}"
            mapping[individual_key] = species_label

        mapping_file_path = os.path.join(subdir, "mapping.txt")
        with open(mapping_file_path, 'w') as f:
            for key, value in mapping.items():
                f.write(f"{key} {value}\n")

        print("Processed and saved mapping file for", iteration_folder, "in:", mapping_file_path)
