from treelib import Tree
from build_raw_nodes_from_newick import nodes_from_newick_file
from VSim_test.tree_functions import DataOnNode


def MakeAnyTree(datafile):

    raw_nodes = nodes_from_newick_file(datafile)

    tree = Tree()

    print('STARTED TREE BUILDING')

    def add_children(some_tree, node):  # adds all children to the tree by checking all parent ids
        if node["parentid"] != -1:
            some_tree.create_node(node["name"], node["id"], parent=node["parentid"],
                                  data = DataOnNode(mutation_name='1', old_nucleotyde='A',
                                                    new_nucleotyde='A', time_of_birth= tree.get_node(node["parentid"]).data.time_of_birth + node["length"]))
        else:
            some_tree.create_node(node["name"], node["id"], parent = None,
                                  data=DataOnNode(mutation_name='1', old_nucleotyde='A',
                                                  new_nucleotyde='A',
                                                  time_of_birth=0))
        for child_node in node["children"]:
            add_children(some_tree, child_node)

    add_children(tree, raw_nodes)

    print('ENDED TREE BUILDING')

    # tree.save2file('tree.txt')
    return tree

builded_tree = MakeAnyTree(datafile="newick_test/data10.txt")
builded_tree.show()

