from parse_newick_function import parse

newick_file = open('/Users/LAB-SCG-125/Desktop/SarsData/uncondensed-final-tree.unamb_rus.CorrectedBL.nwk','r')
newick_string = newick_file.read()
newick_file.close()

tree_in_array_format = parse(newick_string)
a=5