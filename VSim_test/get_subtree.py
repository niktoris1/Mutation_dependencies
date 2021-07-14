class SubtreeCreation:


    def __init__(self, A_nucleotyde, A_cite, B_nucleotyde, B_cite, tree):
        # problem is clear. We have nucleotyde as number in code and as letter here
        def LetterToNumber(number):
            if number == 'A':
                return 0
            if number == 'T':
                return 1
            if number == 'C':
                return 2
            if number == 'G':
                return 3

        self.A_nucleotyde = LetterToNumber(A_nucleotyde)
        self.A_cite = A_cite
        self.B_nucleotyde = LetterToNumber(B_nucleotyde)
        self.B_cite = B_cite
        self.tree = tree


    def GetABsubtrees(self):
        A_roots = []
        B_roots = []
        not_AB_roots = []
        AB_roots = []

        def CheckNode(node):
            if node.data.is_a_mutation == False:
                return 'Nothing'
            elif node.data.mutation_cite == self.A_cite and node.data.new_nucleotyde == self.A_nucleotyde:
                return 'A'
            elif node.data.mutation_cite == self.A_cite and node.data.new_nucleotyde != self.A_nucleotyde:
                return 'not A'
            elif node.data.mutation_cite == self.B_cite and node.data.new_nucleotyde == self.B_nucleotyde:
                return 'B'
            elif node.data.mutation_cite == self.B_cite and node.data.new_nucleotyde != self.B_nucleotyde:
                return 'not B'
            else:
                return 'Nothing'

        def CheckAncestors(needed_type, node): # returns the first ancestor with the needed type of mutation - or false
            if node.identifier == self.tree.root:
                return False
            parent = self.tree.get_node(node.bpointer)
            if CheckNode(parent) == needed_type:
                return parent
            else:
                return CheckAncestors(needed_type, parent)


        for node in self.tree.all_nodes():
            if node.data.mutation_cite == 0 and node.data.new_nucleotyde == 1:
                continue
            if CheckNode(node) == 'A':
                A_roots.append(node)
            elif CheckNode(node) == 'B':
                B_roots.append(node)
            elif CheckNode(node) == 'not A':
                not_AB_roots.append(node)
            elif CheckNode(node) == 'not B':
                not_AB_roots.append(node)


        #need to make a change and identify A_roots, B_roots and AB_roots correctly

        for A_root in A_roots:
            parent_AB_root = CheckAncestors('B', A_root)
            if parent_AB_root != False:
                AB_roots.append(A_root)

        for B_root in B_roots:
            parent_AB_root = CheckAncestors('A', B_root)
            if parent_AB_root != False:
                AB_roots.append(parent_AB_root)

        AB_roots = list(set(AB_roots))
        not_AB_roots = list(set(not_AB_roots))

        subtrees = []
        for good_root in AB_roots:
            subtrees.append(self.tree.subtree(good_root.identifier))

        for bad_root in not_AB_roots:
            for subtree in subtrees:
                if subtree.contains(bad_root.identifier):
                    subtree.remove_node(bad_root.identifier)


        return subtrees





