import math
from get_subtree import test_trees
from events_from_tree import GetEventsFromTree, events_sequence

import treelib

class LikelyhoodEstimation:

    def __init__(self, estimated_tree, events_sequence):
        self.estimated_tree = estimated_tree
        self.events_sequence = events_sequence



    def LLH_function(self, iteration, coal_rate):
        if iteration == 0:
            return 1
        else:
            event = LikelyhoodEstimation.EventFromIteration(self, iteration)
            event_type = event.event_type
            time = LikelyhoodEstimation.TimeFromIteration(self, iteration)
            distinct_lineages = LikelyhoodEstimation.DistinctLineages(self, time)
            event_probability = LikelyhoodEstimation.EventProbability(self, event_type, coal_rate, distinct_lineages)

            return LikelyhoodEstimation.LLH_function(self, iteration - 1, coal_rate, coal_iteration) * \
               math.exp(coal_rate * (LikelyhoodEstimation.TimeFromIteration(self, iteration) - LikelyhoodEstimation.TimeFromIteration(self, coal_iteration - 1))) \
               * math.comb(LikelyhoodEstimation.DistinctLineages(self, LikelyhoodEstimation.TimeFromIteration(self, iteration - 1)), 2) * \
               event_probability


    def TimeFromIteration(self, iteration):
        time = events_sequence[iteration-1].event_time
        return time

    def IterationFromTime(self, time):
        for iteration_number in range(0, len(events_sequence)):
            if time >= events_sequence[iteration_number].event_time and time < events_sequence[iteration_number+1].event_time:
                return iteration_number

    def DistinctLineages(self, time, events_sequence):

        iteration = LikelyhoodEstimation.IterationFromTime(self, time)

        dl = 0
        for iteration_number in range(0, iteration):
            if events_sequence[iteration_number].event_type == "adding_lineage":
                dl = dl+1
            if events_sequence[iteration_number].event_type == "coalescence":
                for individual_tree in self.estimated_tree:
                    if individual_tree.contains(events_sequence[iteration_number].vertex_tag):
                        number_of_children = len(individual_tree.get_node(events_sequence[iteration_number].vertex_tag).fpointer)
                        dl = dl - number_of_children + 1
        return dl

    def EventFromIteration(self, iteration):
        return self.events_sequence[iteration-1]

    def EventProbability(self, event_type, coal_rate, distinct_lineages):
        if event_type == 1:
            return 1
        else:
            probability = 1
            for i in range(0, distinct_lineages):
                probability = probability * coal_rate * math.comb(distinct_lineages - i + 1, 2)


print(events_sequence)

es = events_sequence.event_sequence

le = LikelyhoodEstimation(test_trees, es)
iteration=len(es)-1
coal_rate=1 # will have to estimate
coal_iteration=len(es)
number_of_lineages = len(es)
event_type = es[len(es) - 1].event_type
event_probability = le.EventProbability(le, event_type, coal_rate)

LLH = le.LLH_function(iteration=iteration, coal_rate=coal_rate)


