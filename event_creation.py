import math
from build_tree import covid_tree
from build_tree import name_id_dict

import treelib

class Event:
    def __init__(self, event_type, number_of_lineages, lineages):
        self.event_type = event_type
        self.number_of_lineages = number_of_lineages
        self.lineages = lineages

        if event_type is not 'adding_lineage' or 'coalescence':
            raise('Error - wrong event type')
        if event_type is 'adding_lineage' and len(lineages) is not 1:
            raise('Error - trying to add more or less than 1 lineage')

class EventSequence:
    def __init__(self, event_sequence):
        self.event_sequence = event_sequence





















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












