from treelib import Tree

from tree_functions import MutationOnNode



def ArrayTreeToTreeClass(array_tree, times):

    tree = Tree()

    tree.create_node("0", "0", data=MutationOnNode(mutation_name="0", old_nucleotyde="G", new_nucleotyde="T", time_of_birth=times[0]))

    for i in range(1, len(array_tree)):
        tree.create_node(str(i), str(i), parent=str(array_tree[i]), data=MutationOnNode(mutation_name=str(i), old_nucleotyde="X", new_nucleotyde="W", time_of_birth=times[i]))

    return tree




