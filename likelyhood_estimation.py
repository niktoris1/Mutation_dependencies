import math
from scipy import optimize
#from get_subtree import test_trees
from events_from_tree import GetEventsFromTree, test_events_sequence, test_tree
from treelib import Tree

import matplotlib.pyplot as plt

import treelib

class LikelyhoodEstimation:

    def __init__(self, estimated_tree, events_sequence):
        if isinstance(estimated_tree, list):
            self.estimated_tree = estimated_tree
        elif isinstance(estimated_tree, Tree):
            self.estimated_tree = [estimated_tree]
        else:
            raise("Error - we have on the input nor tree or the list of trees. Check the input, please.")

        self.events_sequence = events_sequence



    def LLH_function(self, iteration, coal_rate, events_sequence):
        if iteration > 0:
            event = LikelyhoodEstimation.EventFromIteration(self, iteration, events_sequence)
            event_type = event.event_type
            time = LikelyhoodEstimation.TimeFromIteration(self, iteration, events_sequence)
            current_time = LikelyhoodEstimation.TimeFromIteration(self, iteration, events_sequence)
            previous_time = LikelyhoodEstimation.TimeFromIteration(self, iteration - 1, events_sequence)
            distinct_lineages = LikelyhoodEstimation.DistinctLineages(self, time, events_sequence)
            event_probability = LikelyhoodEstimation.EventProbability(self, iteration, event_type, coal_rate, events_sequence)


            value_on_this_step = LikelyhoodEstimation.LLH_function(self, iteration - 1, coal_rate, events_sequence) * \
               math.exp(- coal_rate * (current_time - previous_time) * \
               math.comb(distinct_lineages, 2)) * event_probability

            #print("distinct_lineages = ", distinct_lineages)
            #print("event_probability = ", event_probability)
            #print("value_on_this_step = ", value_on_this_step)
            #print("event_type = ", event_type)
            #print("-----------------")

            return value_on_this_step
        if iteration == 0:
            return 1



    def TimeFromIteration(self, iteration, events_sequence):
        if iteration < 0:
            return 0
        time = events_sequence[iteration].event_time
        return time

    def IterationFromTime(self, time, events_sequence):
        for iteration_number in range(0, len(events_sequence) - 1):
            if time >= events_sequence[iteration_number].event_time and time < events_sequence[iteration_number+1].event_time:
                return iteration_number
        return len(events_sequence) - 1 #if the last iteration

    def DistinctLineages(self, time, events_sequence):

        iteration = LikelyhoodEstimation.IterationFromTime(self, time, events_sequence)

        dl = 1
        for iteration_number in range(0, iteration+1):
            if events_sequence[iteration_number].event_type == "adding_lineage":
                dl = dl - 1
            if events_sequence[iteration_number].event_type == "coalescence":
                for individual_tree in self.estimated_tree:
                    if individual_tree.contains(events_sequence[iteration_number].vertex_tag):
                        number_of_children = len(individual_tree.get_node(events_sequence[iteration_number].vertex_tag).fpointer)
                        dl = dl + number_of_children - 1
        return dl

    def EventFromIteration(self, iteration, events_sequence):
        return events_sequence[iteration]

    def EventProbability(self, iteration, event_type, coal_rate, events_sequence):

        time = LikelyhoodEstimation.TimeFromIteration(self, iteration, events_sequence)
        distinct_lineages = LikelyhoodEstimation.DistinctLineages(self, time, events_sequence)

        if event_type == "adding_lineage":
            return 1
        else:
            probability = 1
            for i in range(0, distinct_lineages):
                probability = probability * coal_rate * math.comb(distinct_lineages - i + 1, 2)
            return probability


number_of_events = len(test_events_sequence)

le = LikelyhoodEstimation(test_tree, test_events_sequence)

test_tree.show()

iteration=number_of_events - 1
coal_rate=0.5 # will have to estimate

LLH = le.LLH_function(iteration=iteration, coal_rate=coal_rate, events_sequence=test_events_sequence)

wrapper = lambda coal_rate, iteration, test_events_sequence: - le.LLH_function(iteration=iteration, coal_rate=coal_rate, events_sequence=test_events_sequence)

LLH_optimised = optimize.minimize(fun=wrapper, x0=coal_rate, args=(iteration, test_events_sequence), method='Nelder-Mead')

print(LLH_optimised.fun, LLH_optimised.x)

x = [x / 100.0 for x in range(1, 200, 1)]
y = [le.LLH_function(iteration=iteration, coal_rate=i, events_sequence=test_events_sequence) for i in x]

plt.plot(x, y)

plt.show()

print(LLH_optimised.fun, LLH_optimised.x)




