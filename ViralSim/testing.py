from treelib import Tree

def ArrayTreeToTreeClass(array_tree):

    tree = Tree()

    tree.create_node(0, 0, data=None)

    for i in range(1, len(array_tree)):
        tree.create_node(i, i, parent=array_tree[i], data=None)

    return tree

