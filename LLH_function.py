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



def GetABsubtrees(A_nucleotyde, A_cite, B_nucleotyde, B_cite, covid_tree):
    for node in covid_tree.all_nodes():
        if node.data is not None:
            node.data = [node.data] + [0] + [0] # here first element corresponds to node having an A_nucleotyde and second - to node having a B_nucleotyde in B_cite

    root_node = tree.get_node(covid_tree.root)
    ProcessNode(root_node, A_nucleotyde, A_cite, B_nucleotyde, B_cite)

    AB_haplotype_node_list = []
    list_of_AB_roots = CheckIfABHaplotype(root_node, AB_haplotype_node_list)

    for node in covid_tree.all_nodes():
        if node.data is not None:
            node.data = node.data[0] # clean the data

    return list_of_AB_roots

def ProcessNode(node, A_nucleotyde, A_cite, B_nucleotyde, B_cite):
    if node.data is not None:
        if node.data[0].mutation_cite == A_cite and node.data[0].new_nucleotyde == A_nucleotyde:
            node.data[1] = 1
        if node.data[0].mutation_cite == B_cite and node.data[0].new_nucleotyde == B_nucleotyde:
            node.data[2] = 1
        if node.data[0].mutation_cite == A_cite and node.data[0].new_nucleotyde != A_nucleotyde:
            node.data[1] = 0
        if node.data[0].mutation_cite == B_cite and node.data[0].new_nucleotyde != B_nucleotyde:
            node.data[2] = 0

    for child_node in tree.children(node.identifier):
        ProcessNode(child_node, A_nucleotyde, A_cite, B_nucleotyde, B_cite)


def CheckIfABHaplotype(node, AB_haplotype_node_list):
    if node.data is not None:
        if node.data[1] == 1 and node.data[2] == 1:
            AB_haplotype_node_list = AB_haplotype_node_list + [node]
        else:
            for child in tree.children(node.identifier):
                CheckIfABHaplotype(child, AB_haplotype_node_list)
    else:
        for child in tree.children(node.identifier):
            CheckIfABHaplotype(child, AB_haplotype_node_list)

    return AB_haplotype_node_list




class TestHaplotypeAB:
    def __init__(self, A_nucleotyde, A_cite, B_nucleotyde, B_cite, phylogenetic_tree):
        self.A_nucleotyde = A_nucleotyde
        self.A_cite = A_cite
        self.B_nucleotyde = B_nucleotyde
        self.B_cite = B_cite
        self.phylogenetic_tree = phylogenetic_tree

    def getLikelyhoodAB(self):
        return 0

a = GetABsubtrees('G', 270, 'G', 270, tree)
print(a)

for node in tree.all_nodes():
    if node.data is not None:
        a = node.data.mutation_cite








