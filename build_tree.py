from treelib import Node, Tree
from initialize_tree import raw_nodes
from parse_mutations import mutations
from tree_functions import add_children, name_id_dict, MutationOnNode
import re




covid_tree = Tree()
print('STARTED TREE BUILDING')

covid_tree.create_node(raw_nodes["name"], raw_nodes["id"]) # we build a tree from newick here
for children_node in raw_nodes["children"]:
    add_children(covid_tree, children_node)

print('ENDED TREE BUILDING')


#tree.save2file('tree.txt')

for mutation in mutations: # here for every node we write a mutation in it
    for mutation_node in mutation.mutation_nodes: # we take all possible couples of old-new nycleotydes
        for node_name in mutation_node.nodes_names:
            if node_name[-3:] == '[F]' or node_name[-3:] == '[B]': # remove forward-backward notation
                node_name = node_name[:-3]
            chosen_node = covid_tree.get_node(name_id_dict[node_name])
            chosen_node.data = MutationOnNode(mutation_name=mutation.mutation_name, old_nucleotyde=mutation_node.old_nucleotyde, new_nucleotyde=mutation_node.new_nucleotyde, time_of_birth='Unknown')
            # Костыль. Мы оставляем на ноде в дереве только самые важные данные - название мутации, место мутaции, стратовый и новый нуклеотид и время до родителя




