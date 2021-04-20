import math
from scipy import optimize
#from get_subtree import test_trees
from tree_functions import GetEventsFromTree
from treelib import Tree

import matplotlib.pyplot as plt

import logging # a workaround to kill warnings
logging.captureWarnings(True)

class LikelyhoodEstimation:

    def __init__(self, estimated_tree):
        if isinstance(estimated_tree, list):
            if len(estimated_tree) == 0:
                raise ('Trying to estimate empty subtree - undefined. Seems that this mutation simply does not exist.')
            self.estimated_tree = estimated_tree
        elif isinstance(estimated_tree, Tree):
            self.estimated_tree = [estimated_tree]
        else:
            raise("Error - we have on the input neither tree or the list of trees. Check the input, please.")

        self.events_sequence = GetEventsFromTree(self.estimated_tree)
        self.number_of_events = len(self.events_sequence)
        self.distinct_lineages = self.BuildDistinctLineages(self.events_sequence)



    def LLH_function(self, coal_rate):

        LLH_values = [0] * self.number_of_events
        events = [0] * self.number_of_events
        event_type_array = [0] * self.number_of_events
        time_array = [0] * self.number_of_events
        current_time_array = [0] * self.number_of_events
        #previous_time_array = [0] * self.number_of_events
        distinct_lineages_array = [0] * self.number_of_events
        event_probability_array = [0] * self.number_of_events
        addition = [0] * self.number_of_events


        for iteration in range(0, self.number_of_events):
            events[iteration] = LikelyhoodEstimation.EventFromIteration(self, iteration)
            event_type_array[iteration] = events[iteration].event_type
            time_array[iteration] = LikelyhoodEstimation.TimeFromIteration(self, iteration)
            current_time_array[iteration] = LikelyhoodEstimation.TimeFromIteration(self, iteration)
            #previous_time_array[iteration] = LikelyhoodEstimation.TimeFromIteration(self, iteration - 1)
            distinct_lineages_array[iteration] = LikelyhoodEstimation.DistinctLineages(self, current_time_array[iteration])
            event_probability_array[iteration] = LikelyhoodEstimation.EventProbability(self, current_time_array[iteration], events[iteration], coal_rate)

        for iteration in range(1, self.number_of_events):
            addition[iteration] = (- coal_rate * (current_time_array[iteration] - current_time_array[iteration - 1]) * \
               math.comb(distinct_lineages_array[iteration], 2)) + math.log(event_probability_array[iteration])

            LLH_values[iteration] = LLH_values[iteration - 1] + addition[iteration]

            #print("distinct_lineages = ", distinct_lineages)
            #print("event_probability = ", event_probability)
            #print("value_on_this_step = ", value_on_this_step)
            #print("event_type = ", event_type)
            #print("-----------------")

            #print(value_on_this_step)

        #plt.show()


        return LLH_values[-1]


    def TimeFromIteration(self, iteration):
        if iteration < 0:
            return 0
        return self.events_sequence[iteration].event_time

    def IterationFromTime(self, time):

        def IterationFromTimeStartFinish(time, start, finish):
            middle = (start + finish) // 2

            if time <= self.events_sequence[start].event_time:
                return start
            if time >= self.events_sequence[finish].event_time:
                return finish

            if time == self.events_sequence[middle].event_time:
                return middle
            if time > self.events_sequence[middle].event_time:
                return IterationFromTimeStartFinish(time, middle + 1, finish)
            else:
                return IterationFromTimeStartFinish(time, start, middle - 1)

        return IterationFromTimeStartFinish(time, 0, len(self.events_sequence) - 1)


    def BuildDistinctLineages(self, event_sequence): #returns number of distinct lineages which corresponds to the event sequence

        distinct_lineages = [0] * len(event_sequence)

        for iteration in range(0, len(event_sequence)):
            if iteration > 0:
                distinct_lineages[iteration] = distinct_lineages[iteration-1] - 1
            else:
                distinct_lineages[iteration] = -1

            for individual_tree in self.estimated_tree:
                if individual_tree.contains(event_sequence[iteration].vertex_tag):
                    if event_sequence[iteration].event_type == "coalescence":
                        number_of_children = len(individual_tree.get_node(event_sequence[iteration].vertex_tag).fpointer)
                        distinct_lineages[iteration] = distinct_lineages[iteration] + number_of_children
                            #print('Added ', number_of_children, 'children')
                    if individual_tree.root == event_sequence[iteration].vertex_tag:
                        distinct_lineages[iteration] = distinct_lineages[iteration] + 1
                            #print('Added ', 1, 'for root')
                    break

            if distinct_lineages[iteration] < 0:
                raise Exception('Error, less than zero lineages!')

                #print('Currently', dl, 'distinct lineages')

        return distinct_lineages

    def DistinctLineages(self, time):
        iteration = self.IterationFromTime(time)
        return self.distinct_lineages[iteration]

    def EventFromIteration(self, iteration):
        return self.events_sequence[iteration]

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





