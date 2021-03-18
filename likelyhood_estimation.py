import math
from scipy import optimize
#from get_subtree import test_trees
from tree_functions import GetEventsFromTree
from treelib import Tree

import matplotlib.pyplot as plt

import treelib

import logging # a workaround to kill warnings
logging.captureWarnings(True)

class LikelyhoodEstimation:

    def __init__(self, estimated_tree):
        if isinstance(estimated_tree, list):
            if len(estimated_tree) == 0:
                raise ('Trying to estimate empty subree - undefined. Seems that this mutation simply does not exist.')
            self.estimated_tree = estimated_tree
        elif isinstance(estimated_tree, Tree):
            self.estimated_tree = [estimated_tree]
        else:
            raise("Error - we have on the input neither tree or the list of trees. Check the input, please.")

        self.events_sequence = GetEventsFromTree(self.estimated_tree)

        self.number_of_events = len(self.events_sequence)



    def LLH_function(self, coal_rate):

        LHH_values = [0] * self.number_of_events
        events = [0] * self.number_of_events
        event_type_array = [0] * self.number_of_events
        time_array = [0] * self.number_of_events
        current_time_array = [0] * self.number_of_events
        previous_time_array = [0] * self.number_of_events
        distinct_lineages_array = [0] * self.number_of_events
        event_probability_array = [0] * self.number_of_events

        for iteration in range(0, self.number_of_events):
            events[iteration] = LikelyhoodEstimation.EventFromIteration(self, iteration)
            event_type_array[iteration] = events[iteration].event_type
            time_array[iteration] = LikelyhoodEstimation.TimeFromIteration(self, iteration)
            current_time_array[iteration] = LikelyhoodEstimation.TimeFromIteration(self, iteration)
            previous_time_array[iteration] = LikelyhoodEstimation.TimeFromIteration(self, iteration - 1)
            distinct_lineages_array[iteration] = LikelyhoodEstimation.DistinctLineages(self, current_time_array[iteration])
            event_probability_array[iteration] = LikelyhoodEstimation.EventProbability(self, current_time_array[iteration], event_type_array[iteration], coal_rate)

        for iteration in range(1, self.number_of_events):
            LHH_values[iteration] = LHH_values[iteration - 1] + \
               (- coal_rate * (current_time_array[iteration] - previous_time_array[iteration]) * \
               math.comb(distinct_lineages_array[iteration], 2)) + math.log( event_probability_array[iteration])

            #print("distinct_lineages = ", distinct_lineages)
            #print("event_probability = ", event_probability)
            #print("value_on_this_step = ", value_on_this_step)
            #print("event_type = ", event_type)
            #print("-----------------")

            #print(value_on_this_step)

        return LHH_values[-1]


    def TimeFromIteration(self, iteration):
        if iteration < 0:
            return 0
        return self.events_sequence[iteration].event_time

    def IterationFromTime(self, time):
        for iteration_number in range(0, len(self.events_sequence) - 1):
            if time >= self.events_sequence[iteration_number].event_time and time < self.events_sequence[iteration_number+1].event_time:
                return iteration_number
        return len(self.events_sequence) - 1 #if the last iteration

    def DistinctLineages(self, time):

        iteration = LikelyhoodEstimation.IterationFromTime(self, time)

        dl = 0 # Число деревьев?
        for iteration_number in range(0, iteration):
            event = self.events_sequence[iteration_number]

            for individual_tree in self.estimated_tree:
                if individual_tree.contains(event.vertex_tag):
                    if event.event_type == "adding_lineage":
                        if individual_tree.contains(individual_tree.get_node(event.vertex_tag).bpointer): # not root
                            dl = dl - 1
                        else:
                            pass
                        break

                    if event.event_type == "coalescence":
                        number_of_children = len(individual_tree.get_node(event.vertex_tag).fpointer)
                        if individual_tree.contains(individual_tree.get_node(event.vertex_tag).bpointer): # not root
                            dl = dl + number_of_children - 1
                        else:
                            dl = dl + number_of_children # if the tree is root - we add all of the children as a new lineage. If not - we add children and remove ancestor lineage
                        break
        return dl

    def EventFromIteration(self, iteration):
        return self.events_sequence[iteration]

    def EventProbability(self, time, event_type, coal_rate):

        distinct_lineages = LikelyhoodEstimation.DistinctLineages(self, time)

        if event_type == "adding_lineage":
            return 1
        else:
            probability = 1
            for i in range(0, distinct_lineages):
                probability = probability * coal_rate * math.comb(distinct_lineages - i + 1, 2)
            return probability


    def GetEstimation(self):

        print("GETTING ESTIMATION")
        start_coal_rate=20 # will have to estimate
        wrapper = lambda coal_rate: - self.LLH_function(coal_rate=coal_rate)

        LLH_optimised = optimize.minimize_scalar(fun=wrapper, bounds = (0.01, 1000), bracket = (10, 100), method='Bounded', tol=1e-1)
        #LLH_optimised = optimize.minimize(fun=wrapper, x0=20, method='Nelder-Mead', tol=1e-2)

        #x = [x / 100.0 for x in range(1, 100, 1)]
        #y = [LikelyhoodEstimation.LLH_function(self, iteration=iteration, coal_rate=i) for i in x]

        #plt.plot(x, y)
        #plt.show()

        LLH_result = LLH_optimised.fun
        LLH_point = LLH_optimised.x

        print(LLH_result, LLH_point)

        return [LLH_result, LLH_point]




#le = LikelyhoodEstimation(test_tree)
#test_tree.show()

#le.GetEstimation()





