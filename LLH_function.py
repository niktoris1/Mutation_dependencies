import math
from build_tree import tree
from build_tree import name_id_dict

def LLH_function(iteration, coal_rate, coal_iteration, number_of_lineages, event_probability, event_type):
    return LLH_function(iteration-1) * math.exp(coal_rate * (TimeFromIteration(coal_iteration) - TimeFromIteration(coal_iteration - 1))) \
           * math.comb(DistinctLineages(TimeFromIteration(coal_iteration - 1)), 2) * EventProbability(event_type, coal_rate, DistinctLineages(TimeFromIteration(iteration)))




def TimeFromIteration(iteration)
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

class SubTree:
    def __init__(self, tree_origin): ## subtree is defined by its origin node, old and new nucleotyde variants
        self.origin = tree_origin

    def GetSubtree(self):
        return tree.subtree(name_id_dict[self.tree_origin])

    def GetLeaves(self):
        return tree.leaves(SubTree.GetSubtree(self))



