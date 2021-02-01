import math
from build_tree import covid_tree
from build_tree import name_id_dict

import treelib

class SubtreeCreation:

    def ProcessNode(self, root_node, A_nucleotyde, A_cite, B_nucleotyde, B_cite, current_A_state, current_B_state):
        if root_node.data is not None:
            if root_node.data[0].mutation_cite == A_cite and root_node.data[0].new_nucleotyde == A_nucleotyde:
                root_node.data[1] = 1
                current_A_state = 1
            if root_node.data[0].mutation_cite == B_cite and root_node.data[0].new_nucleotyde == B_nucleotyde:
                root_node.data[2] = 1
                current_B_state = 1
            if root_node.data[0].mutation_cite == A_cite and root_node.data[0].new_nucleotyde != A_nucleotyde:
                root_node.data[1] = 0
                current_A_state = 0
            if root_node.data[0].mutation_cite == B_cite and root_node.data[0].new_nucleotyde != B_nucleotyde:
                root_node.data[2] = 0
                current_B_state = 0


        for child_node in covid_tree.children(root_node.identifier):
            SubtreeCreation.ProcessNode(self, child_node, A_nucleotyde, A_cite, B_nucleotyde, B_cite, current_A_state, current_B_state)

    def CheckIfABHaplotype(self, node):

        if node.data is not None:
            if node.data[1] == 1 and node.data[2] == 1:
                SubtreeCreation.AB_haplotype_node_list.append(node)
        for child in covid_tree.children(node.identifier):
            SubtreeCreation.CheckIfABHaplotype(self, child)

        return SubtreeCreation.AB_haplotype_node_list

    AB_haplotype_node_list = [] # a local variable for the class

    def GetABsubtrees(self, A_nucleotyde, A_cite, B_nucleotyde, B_cite, base_tree): # here we do not look at what the old haplotype is
        for node in base_tree.all_nodes():
            if node.data is not None:
                node.data = [node.data] + [0] + [0] # here first element corresponds to node having an A_nucleotyde and second - to node having a B_nucleotyde in B_cite

        root_node = base_tree.get_node(base_tree.root)


        if root_node.data is None:
            starting_A_state = 0
            starting_B_state = 0
        else:
            if root_node.data.mutation_cite is not A_cite:
                starting_A_state = 0
            else:
                starting_A_state = 1

            if root_node.data.mutation_cite is not B_cite:
                starting_B_state = 0
            else:
                starting_B_state = 1

        SubtreeCreation.ProcessNode(self, root_node, A_nucleotyde, A_cite, B_nucleotyde, B_cite, starting_A_state, starting_B_state) # at the root there are no mutations

        list_of_AB_roots = SubtreeCreation.CheckIfABHaplotype(self, root_node)

        for node in base_tree.all_nodes():
            if node.data is not None:
                node.data = node.data[0] # clean the data

        return list_of_AB_roots




SomeCreation = SubtreeCreation()

a = SomeCreation.GetABsubtrees('A', 27502, 'A', 27502, covid_tree)
print(a)
print(a[0].data.mutation_cite)