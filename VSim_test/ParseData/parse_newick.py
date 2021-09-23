from parse_newick_functions import parse_newick_string, TreeArrayToTreeClass

newick_file = open('/Users/LAB-SCG-125/Desktop/SarsData/uncondensed-final-tree.unamb_rus.CorrectedBL.nwk','r')
newick_string = newick_file.read()
newick_file.close()

tree_in_array_format = parse_newick_string(newick_string)
tree_in_tree_format = TreeArrayToTreeClass(tree_in_array_format)

tree_size, haps = 0, []
for node in tree_in_tree_format.all_nodes():
    if node.data['haplotype'] != 'Unknown':
        tree_size += 1
        haps.append(node.data['haplotype'])
number_of_haps = len(set(haps))
print(tree_size, number_of_haps)