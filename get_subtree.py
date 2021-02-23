from treelib import Node, Tree
from build_tree import MutationOnNode


class SubtreeCreation:

    AB_haplotype_subtree_roots = []
    not_AB_haplotype_subtree_roots = []

    def ProcessNode(self, root_node, A_nucleotyde, A_cite, B_nucleotyde, B_cite, current_A_state, current_B_state, base_tree):
        # here we mark all nodes as [0, 0] if no A_nucleotyde are in A_cite and no A_nycleotyde are in B_cite

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

            root_node.data[1] = current_A_state
            root_node.data[2] = current_B_state

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
        #format is [[AB_roots], [not_AB_roots]]

        AB_subtrees = SubtreeCreation.SubtreesFromRoots(self, list_of_AB_roots[0], list_of_AB_roots[1], base_tree)

        for node in base_tree.all_nodes():
            if node.data is not None:
                node.data = node.data[0] # clean the data

        return AB_subtrees

    def SubtreesFromRoots(self, AB_haplotype_subtree_roots, not_AB_haplotype_subtree_roots, base_tree):

        subtrees = []
        for good_root in AB_haplotype_subtree_roots:
            subtrees.append(base_tree.subtree(good_root.identifier))

        for bad_root in not_AB_haplotype_subtree_roots:
            for subtree in subtrees:
                if subtree.contains(bad_root.identifier):
                    subtree.remove_node(bad_root.identifier)

        return subtrees


SomeCreation = SubtreeCreation()

test_tree = Tree()

test_tree.create_node(1, 1, data=None)
test_tree.create_node(2, 2, parent = 1, data=MutationOnNode(mutation_name="2A", old_nucleotyde="G", new_nucleotyde="T", time_of_birth=1))
test_tree.create_node(3, 3, parent = 1, data=MutationOnNode(mutation_name="3", old_nucleotyde="G", new_nucleotyde="T", time_of_birth=2))
test_tree.create_node(4, 4, parent = 2, data=MutationOnNode(mutation_name="4", old_nucleotyde="G", new_nucleotyde="T", time_of_birth=3))
test_tree.create_node(5, 5, parent = 2, data=MutationOnNode(mutation_name="5", old_nucleotyde="G", new_nucleotyde="T", time_of_birth=4))
test_tree.create_node(6, 6, parent = 4, data=MutationOnNode(mutation_name="2B", old_nucleotyde="T", new_nucleotyde="G", time_of_birth=5))
test_tree.create_node(7, 7, parent = 3, data=MutationOnNode(mutation_name="2C", old_nucleotyde="G", new_nucleotyde="T", time_of_birth=11))
test_tree.create_node(8, 8, parent = 4, data=MutationOnNode(mutation_name="6", old_nucleotyde="T", new_nucleotyde="G", time_of_birth=10))
test_tree.create_node(9, 9, parent = 3, data=MutationOnNode(mutation_name="7", old_nucleotyde="G", new_nucleotyde="T", time_of_birth=2))


#test_trees = SomeCreation.GetABsubtrees('T', 2, 'T', 2, test_tree)

#for tree in test_trees:
#    tree.show()
