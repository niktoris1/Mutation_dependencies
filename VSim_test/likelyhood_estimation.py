import math
from scipy import optimize
from VSim_test.tree_functions import GetEventsFromTree, IterationFromTime
from treelib import Tree

import matplotlib.pyplot as plt

import logging # a workaround to kill warnings
logging.captureWarnings(True)

class LikelyhoodEstimation:

    def __init__(self, estimated_tree=None, es=None):
        if isinstance(estimated_tree, list):
            if len(estimated_tree) == 0:
                raise ('Trying to estimate empty subtree - undefined. Seems that this mutation simply does not exist.')
            self.estimated_tree = estimated_tree
        elif isinstance(estimated_tree, Tree):
            self.estimated_tree = [estimated_tree] # making a single tree a list for the sake of consitency
        else:
            raise("Error - we have on the input neither tree or the list of trees. Check the input, please.")

        if estimated_tree != None:
            self.es = GetEventsFromTree(self.estimated_tree)
        else:
            self.es = es
        self.number_of_events = len(self.es.sequence)
        self.distinct_lineages = self.BuildDistinctLineages(self.es.sequence)



    def LLH_function(self, coal_rate):

        LLH_values = [0] * self.number_of_events
        events = [0] * self.number_of_events
        event_type_array = [0] * self.number_of_events
        current_time_array = [0] * self.number_of_events
        distinct_lineages_array = [0] * self.number_of_events
        event_probability_array = [0] * self.number_of_events
        addition = [0] * self.number_of_events


        for iteration in range(self.number_of_events):
            events[iteration] = LikelyhoodEstimation.EventFromIteration(self, iteration)
            event_type_array[iteration] = events[iteration].event_type
            current_time_array[iteration] = self.es.TimeFromIteration(iteration)
            distinct_lineages_array[iteration] = LikelyhoodEstimation.DistinctLineages(self, current_time_array[iteration])
            event_probability_array[iteration] = LikelyhoodEstimation.EventProbability(self, current_time_array[iteration], events[iteration], coal_rate)

        for iteration in range(1, self.number_of_events):
            addition[iteration] = (- coal_rate * (current_time_array[iteration] - current_time_array[iteration - 1]) *
               math.comb(distinct_lineages_array[iteration], 2)) + math.log(event_probability_array[iteration])

            LLH_values[iteration] = LLH_values[iteration - 1] + addition[iteration]

        return LLH_values[-1]


    def BuildDistinctLineages(self, event_sequence):
        # returns number of distinct lineages which corresponds to the times at the event sequence

        distinct_lineages = [0] * len(event_sequence)

        for iteration in range(len(event_sequence)):
            if iteration > 0:
                distinct_lineages[iteration] = distinct_lineages[iteration-1] - 1
            else:
                distinct_lineages[iteration] = -1

            for individual_tree in self.estimated_tree:
                if individual_tree.contains(event_sequence[iteration].vertex_id):
                    if event_sequence[iteration].event_type == "coalescence":
                        number_of_children = len(individual_tree.get_node(event_sequence[iteration].vertex_id).fpointer)
                        distinct_lineages[iteration] = distinct_lineages[iteration] + number_of_children
                            #print('Added ', number_of_children, 'children')
                    if individual_tree.root == event_sequence[iteration].vertex_id:
                        distinct_lineages[iteration] = distinct_lineages[iteration] + 1
                            #print('Added ', 1, 'for root')
                    break

            if distinct_lineages[iteration] < 0:
                raise Exception('Error, less than zero lineages!')

                #print('Currently', dl, 'distinct lineages')

        return distinct_lineages

    def DistinctLineages(self, time):
        iteration = IterationFromTime(time, es=self.es)
        return self.distinct_lineages[iteration]

    def EventFromIteration(self, iteration):
        return self.es.sequence[iteration]

    def EventProbability(self, time, event, coal_rate):

        distinct_lineages = LikelyhoodEstimation.DistinctLineages(self, time)

        if event.event_type == "adding_lineage":
            return 1
        else:
            probability = 1
            for i in range(0, event.number_of_children - 1): # ???
                probability = probability * coal_rate * math.comb(distinct_lineages - i + 1, 2)
            return probability


    def GetEstimation(self):

        print("GETTING ESTIMATION")
        #start_coal_rate=20 # will have to estimate
        wrapper = lambda coal_rate: - self.LLH_function(coal_rate=coal_rate)

        LLH_optimised = optimize.minimize_scalar(fun=wrapper, bounds = (0.0001, 1000), bracket = (0.0001, 10), method='Bounded', tol=1e-5)
        #LLH_optimised = optimize.minimize(fun=wrapper, x0=20, method='Nelder-Mead', tol=1e-2)

        x = [x /1000.0 for x in range(1, 1000, 1)]
        y = [LikelyhoodEstimation.LLH_function(self, coal_rate=i) for i in x]

        plt.plot(x, y)
        plt.show()



        LLH_result = LLH_optimised.fun
        LLH_point = LLH_optimised.x

        print(LLH_result, LLH_point)

        return [LLH_result, LLH_point]




#le = LikelyhoodEstimation(test_tree)
#test_tree.show()

#le.GetEstimation()





