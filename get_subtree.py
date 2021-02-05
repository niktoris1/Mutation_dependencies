import math
from build_tree import covid_tree
from build_tree import name_id_dict

import treelib

from treelib import tree

class SubtreeCreation:


    AB_haplotype_subtree_roots = []
    not_AB_haplotype_subtree_roots = []


    def ProcessNode(self, root_node, A_nucleotyde, A_cite, B_nucleotyde, B_cite, current_A_state, current_B_state, base_tree):
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


        for child_node in base_tree.children(root_node.identifier):
            SubtreeCreation.ProcessNode(self, child_node, A_nucleotyde, A_cite, B_nucleotyde, B_cite, current_A_state, current_B_state, base_tree)

    def ABHaplotypeIsPresent(self, node): #checks if both haplotypes are pressent
        if node.data is None:
            return 0
        else:
            if node.data[1] == 1 and node.data[2] == 1:
                return 1
            else:
                return 0

    def GetParentNode(self, node, base_tree):
        if node.is_root():
            return node

        # noinspection INSPECTION_NAME
        return base_tree.get_node(node.bpointer)

    def CheckIfABHaplotype(self, node, base_tree):

        if node.is_root():
            if SubtreeCreation.ABHaplotypeIsPresent(self, node) == 1:
                SubtreeCreation.AB_haplotype_subtree_roots.append(node)
            else:
                SubtreeCreation.not_AB_haplotype_subtree_roots.append(node)

        else:
            if SubtreeCreation.ABHaplotypeIsPresent(self, node) and not SubtreeCreation.ABHaplotypeIsPresent(self, SubtreeCreation.GetParentNode(self, node, base_tree)):
                SubtreeCreation.AB_haplotype_subtree_roots.append(node)
            if not SubtreeCreation.ABHaplotypeIsPresent(self, node) and SubtreeCreation.ABHaplotypeIsPresent(self, SubtreeCreation.GetParentNode(self, node, base_tree)):
                SubtreeCreation.not_AB_haplotype_subtree_roots.append(node)



        for child in base_tree.children(node.identifier):
            SubtreeCreation.CheckIfABHaplotype(self, child, base_tree)

        return [SubtreeCreation.AB_haplotype_subtree_roots, SubtreeCreation.not_AB_haplotype_subtree_roots] # we return two arrays - the roots of defining subtrees



    def GetABsubtrees(self, A_nucleotyde, A_cite, B_nucleotyde, B_cite, base_tree): # here we do not look at what the old haplotype is

        # also, this part of code can be done better by simply getting said mutations from data. Still, no harm done.

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

        SubtreeCreation.ProcessNode(self, root_node, A_nucleotyde, A_cite, B_nucleotyde, B_cite, starting_A_state, starting_B_state, base_tree) # at the root there are no mutations

        list_of_AB_roots = SubtreeCreation.CheckIfABHaplotype(self, root_node, base_tree)

        for node in base_tree.all_nodes():
            if node.data is not None:
                node.data = node.data[0] # clean the data

        return list_of_AB_roots

    def SubtreesFromRoots(self, AB_haplotype_subtree_roots, not_AB_haplotype_subtree_roots, base_tree):

        subtrees = []
        for good_root in AB_haplotype_subtree_roots:
            subtrees.append(base_tree.subtree(good_root.id))

        for bad_root in not_AB_haplotype_subtree_roots:
            for subtree in subtrees:
                subtree.remove_node(bad_root.id)

        return subtrees




SomeCreation = SubtreeCreation()

a = SomeCreation.GetABsubtrees('A', 27502, 'A', 27502, covid_tree)
print(a)
print(a[0][0].data.mutation_cite)