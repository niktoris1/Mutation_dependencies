from treelib import Node, Tree
from initialize_tree import  raw_nodes
from parse_mutations import mutations
import re


def getcite(mutation_name): # gets a cite name from the name of mutation
    return re.findall('\d+', mutation_name)[0]

class MutationOnNode: # defines mutztion on specific node
    def __init__(self, mutation_name, old_nucleotyde, new_nucleotyde):
        self.mutation_name = mutation_name
        self.mutation_cite = getcite(mutation_name)
        self.old_nucleotyde = old_nucleotyde
        self.new_nucleotyde = new_nucleotyde


name_id_dict = {} # dictionary in type id: name

def add_children(tree, node):
    if node["parentid"] is not None:
        tree.create_node(node["name"], node["id"], parent=node["parentid"])
    else:
        tree.create_node(node["name"], node["id"])
    name_id_dict[node["name"]] = node["id"]
    for child_node in node["children"]:
        add_children(tree, child_node)


tree = Tree()
print('STARTED TREE BUILDING')

tree.create_node(raw_nodes["name"], raw_nodes["id"])
for children_node in raw_nodes["children"]:
    add_children(tree, children_node)


#tree.save2file('tree.txt')

for mutation in mutations: # here for every node we write a mutation in it
    #if len(mutation.mutation_nodes) > 5:
    #    print("Mutation ", mutation.mutation_name, " has ", len(mutation.mutation_nodes), " mutations")
    #    for mutation_node in mutation.mutation_nodes:
    #        print(mutation_node.old_nucleotyde, " maps to ", mutation_node.new_nucleotyde)
    for mutation_node in mutation.mutation_nodes: # we take all possible couples of old-new nycleotydes
        for node_name in mutation_node.nodes_names:
            if node_name[-3:] == '[F]' or node_name[-3:] == '[B]':
                node_name = node_name[:-3]
            chosen_node = tree.get_node(name_id_dict[node_name])
            chosen_node.data = MutationOnNode(mutation.mutation_name, mutation_node.old_nucleotyde, mutation_node.new_nucleotyde)
            # Костыль. Мы оставляем на ноде в дереве только самые важные данные - название мутации, место мутции, мутации, а также стратовый и новый нуклеотид



