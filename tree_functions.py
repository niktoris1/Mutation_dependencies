from treelib import Node, Tree
import re

def getcite(mutation_name): # gets a cite name from the name of mutation
    return int(re.findall('\d+', mutation_name)[0])

class MutationOnNode: # defines a mutation on specific node
    def __init__(self, mutation_name, old_nucleotyde, new_nucleotyde, time_of_birth):
        self.mutation_name = mutation_name
        self.mutation_cite = getcite(mutation_name)
        self.old_nucleotyde = old_nucleotyde
        self.new_nucleotyde = new_nucleotyde
        self.time_of_birth = time_of_birth

def ArrayTreeToTreeClass(array_tree):

    tree = Tree()

    tree.create_node(0, 0, data=None)

    for i in range(1, len(array_tree)):
        tree.create_node(i, i, parent=array_tree[i], data=None)

    return tree