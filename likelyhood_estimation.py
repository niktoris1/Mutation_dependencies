import math
from get_subtree import test_trees
from events_from_tree import GetEventsFromTree, events_sequence

class LikelyhoodEstimation:

    def __init__(self, estimated_tree):
        self.estimated_tree = estimated_tree
        self.lineages_array = [i+1 for i in range(0, len(estimated_tree.all_nodes()))] # array [1, 2, 3, 4, 5...]




    def LLH_function(self, iteration, coal_rate, coal_iteration, number_of_lineages, event_probability, event_type):
        return LikelyhoodEstimation.LLH_function(self, iteration - 1, coal_rate, coal_iteration, number_of_lineages, event_probability, event_type) * \
               math.exp(coal_rate * (LikelyhoodEstimation.TimeFromIteration(coal_iteration) - LikelyhoodEstimation.TimeFromIteration(coal_iteration - 1))) \
               * math.comb(LikelyhoodEstimation.DistinctLineages(LikelyhoodEstimation.TimeFromIteration(coal_iteration - 1)), 2) * \
               LikelyhoodEstimation.EventProbability(event_type, coal_rate, LikelyhoodEstimation.DistinctLineages(LikelyhoodEstimation.TimeFromIteration(iteration)))


    def TimeFromIteration(self, iteration):
        time = events_sequence[iteration].event_time
        return time

    def DistinctLineages(self, iteration):
        return self.lineages_array[iteration]

    def EventFromIteration(self, iteration):
        return events_sequence[iteration]


    def EventProbability(self, event_type, coal_rate, distinct_lineages):
        if event_type == 1:
            return 1
        else:
            probability = 1
            for i in range(0, distinct_lineages):
                probability = probability * coal_rate * math.comb(distinct_lineages - i + 1, 2)


print(events_sequence)

le = LikelyhoodEstimation(test_trees)

iteration=len(events_sequence)
coal_rate=1
coal_iteration=len(events_sequence)
number_of_lineages = len(events_sequence)
event_type = events_sequence[len(events_sequence)].event_type
event_probability = le.EventProbability(le, event_type, coal_rate, number_of_lineages)

LLH = le.LLH_function(le, iteration=iteration, coal_rate=coal_rate, coal_iteration=coal_iteration, number_of_lineages = number_of_lineages, \
                      event_probability = le.EventProbability(le, events_sequence[len(events_sequence)]))


