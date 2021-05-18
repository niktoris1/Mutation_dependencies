from treelib import Tree
from build_tree_from_newick import nodes_from_newick_file


def MakeAnyTree(datafile):

    raw_nodes = nodes_from_newick_file(datafile)

    name_id_dict = {}  # dictionary in type id: name

    def add_children(some_tree, node):  # adds all children to the tree by checking all parent ids
        if node["parentid"] is not None:
            some_tree.create_node(node["name"], node["id"], parent=node["parentid"])
        else:
            some_tree.create_node(node["name"], node["id"])
        name_id_dict[node["name"]] = node["id"]
        for child_node in node["children"]:
            add_children(some_tree, child_node)

    tree = Tree()
    print('STARTED TREE BUILDING')

    tree.create_node(raw_nodes["name"], raw_nodes["id"])  # we build a tree from newick here
    for children_node in raw_nodes["children"]:
        add_children(tree, children_node)

    print('ENDED TREE BUILDING')

    # tree.save2file('tree.txt')
    return tree

tree = MakeAnyTree(datafile="newick_test/data.txt")
print(tree)