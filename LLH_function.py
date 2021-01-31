import math
from build_tree import tree
from build_tree import name_id_dict

import treelib

def LLH_function(iteration, coal_rate, coal_iteration, number_of_lineages, event_probability, event_type):
    return LLH_function(iteration-1) * math.exp(coal_rate * (TimeFromIteration(coal_iteration) - TimeFromIteration(coal_iteration - 1))) \
           * math.comb(DistinctLineages(TimeFromIteration(coal_iteration - 1)), 2) * EventProbability(event_type, coal_rate, DistinctLineages(TimeFromIteration(iteration)))


def TimeFromIteration(iteration):
    time = iteration ## to do
    return time

def DistinctLineages(time):
    number_of_lineages = time
    return number_of_lineages


def EventProbability(event_type, coal_rate, distinct_lineages):
    if event_type == 1:
        return 1
    else:
        probability = 1
        for i in range(0, distinct_lineages):
            probability = probability * coal_rate * math.comb(distinct_lineages - i + 1, 2)
        return probability

def GetABsubtrees(A_haplotype, A_cite, B_haplotype, B_cite, covid_tree):
    for node in covid_tree.all_nodes():
        if node.data is not None:
            node.data = [node.data] + [0] + [0] # here first element corresponds to node having an A_haplotype and second - to node having a B_haplotype in B_cite

    root_node = tree.get_node(covid_tree.root)
    ProcessNode(root_node, A_haplotype, A_cite, B_haplotype, B_cite)

    A_B_haplotype_list = CheckIfABHaplotype(root_node)

    for node in covid_tree.all_nodes():
        if node.data is not None:
            node.data = node.data[0] # clean the data

    return A_B_haplotype_list

def ProcessNode(node, A_haplotype, A_cite, B_haplotype, B_cite):
    if node.data is not None:
        if node.data[0].mutation_cite == A_cite and node.data[0].new_haplotype == A_haplotype:
            node.data[1] = 1
        if node.data[0].mutation_cite == B_cite and node.data[0].new_haplotype == B_haplotype:
            node.data[2] = 1
        if node.data[0].mutation_cite == A_cite and node.data[0].new_haplotype != A_haplotype:
            node.data[1] = 0
        if node.data[0].mutation_cite == B_cite and node.data[0].new_haplotype != B_haplotype:
            node.data[2] = 0

    for child_node in tree.children(node.identifier):
        ProcessNode(child_node, A_haplotype, A_cite, B_haplotype, B_cite)

    return 0

def CheckIfABHaplotype(node):
    if 'ABHaplotypeNodeList' not in locals():
        ABHaplotypeNodeList = []

    if node.data is not None:
        if node.data[1] == 1 and node.data[2] == 1:
            ABHaplotypeNodeList = ABHaplotypeNodeList + node
        else:
            for child in tree.children(node.identifier):
                CheckIfABHaplotype(child)

    if 'ABHaplotypeNodeList' not in locals():
        return []
    else:
        return ABHaplotypeNodeList



class TestHaplotypeAB:
    def __init__(self, A_haplotype, A_cite, B_haplotype, B_cite, phylogenetic_tree):
        self.A_haplotype = A_haplotype
        self.A_cite = A_cite
        self.B_haplotype = B_haplotype
        self.B_cite = B_cite
        self.phylogenetic_tree = phylogenetic_tree

    def getLikelyhoodAB(self):
        return 0

a = GetABsubtrees('G', 100, 'G', 200, tree)

print(a)



