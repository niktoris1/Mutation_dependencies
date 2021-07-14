import math
from scipy import optimize
from VSim_test.tree_functions import GetEventsFromTree, IterationFromTime
from treelib import Tree
from Simulator.TreeDismember import TreeDismemberIO, TreeDismember
from VSim_test.VGsim_test import simulation

import matplotlib.pyplot as plt

import logging # a workaround to kill warnings
logging.captureWarnings(True)



tdm = simulation.gettdm() #get tdm object
trees_funct, trees_neutral = tdm.Dismember() #перед получением таблиц, нужно разчленить дерево
#получение таблиц
event_table_funct, event_table_neutral = tdm.getEventTable() #[{time: [n_samples, n_coals]}]

sample_fraction_table = tdm.getSampleFracTable([0.3, 0.6, 0.9, 1.2, 1.5, 1.8, 2.1, 2.4, 2.7])


class LikelyhoodEstimationDismembered:
    class Event:
        def __init__(self, time, is_sample, is_coal, timestamp, fraction):
            self.time = time
            self.is_sample = is_sample
            self.is_coal = is_coal
            self.timestamp = timestamp
            self.fraction = fraction

    def __init__(self, event_table_funct=None, event_table_neutral=None, sample_fraction_table=None): # here we have not tree, but tables

        # for now we do not care if the tree contains mutation or not
        event_table_funct[0].update(event_table_neutral[0]) # merges dictionaries in-place
        event_table = list(event_table_funct[0].items())

        #format of data is [left timestamp, number of samples, number_of coals, fraction]
        self.bracket_data = []
        for bracket in sample_fraction_table:
            self.bracket_data.append([[], [], [], [bracket, sample_fraction_table[bracket]] ])

        for i in range(len(self.bracket_data)):
            if i < len(self.bracket_data) - 1:
                time_start = self.bracket_data[i][3][0]
                time_finish = self.bracket_data[i+1][3][0]
            else:
                time_start = self.bracket_data[i][3][0]
                time_finish = 999
            for event in event_table:
                if event[0] >= time_start and event[0] < time_finish:
                    self.bracket_data[i][0].append(event[0])
                    self.bracket_data[i][1].append(event[1][0])
                    self.bracket_data[i][2].append(event[1][1])

        x = 5
        # bracket_data is a list of brackets Each correspons to a timeframe
        # each bracket is a list of 4 values
        # 1) list of times
        # 2) list of indicators if the event is sample
        # 3) list of indicators if the event is coalescence



        self.event_array = []
        for i in range(len(self.bracket_data)):
            data_piece = self.bracket_data[i]
            for j in range(len(data_piece[0])):
                event = self.Event(time = data_piece[0][j], is_sample=data_piece[1][j],
                    is_coal = data_piece[2][j], timestamp=data_piece[3][0],
                    fraction=data_piece[3][1])

                self.event_array.append(event)

        self.number_of_events = len(self.event_array)





    def LLH_function(self, coal_rate, needed_timestamp): # this must be rewritten from the tables perspective

        all_events = [self.Event(0, 0, 0, 0, 0)] * self.number_of_events
        # Event(time, is_sample, is_coal, timestamp, fraction)

        events = []
        for i in range(len(all_events)):
            if self.event_array[i].timestamp == needed_timestamp:
                events.append(self.event_array[i])

        def takeTime(elem):
            return elem.time

        events.sort(key=takeTime)

        LLH_values = [0] * len(events)
        distinct_lineages_array = [0] * len(events)
        event_probability_array = [0] * len(events)
        addition = [0] * len(events)

        for i in range(len(events)):
            if i > 0:
                distinct_lineages_array[i] = distinct_lineages_array[i - 1]
            if events[i].is_coal == 1:
                if distinct_lineages_array[i] == 0:
                    distinct_lineages_array[i] = distinct_lineages_array[i] + 2
                else:
                    distinct_lineages_array[i] = distinct_lineages_array[i] + 1






        for i in range(len(events)):
            event_probability_array[i] = LikelyhoodEstimationDismembered.EventProbability(self, coal_rate, i, distinct_lineages_array)

        for i in range(1, len(events)):
            addition[i] = (- coal_rate * (events[i].time - events[i].time) *
               math.comb(distinct_lineages_array[i], 2)) + \
                          math.log(event_probability_array[i])

            LLH_values[i] = LLH_values[i - 1] + addition[i]

        return LLH_values[-1]


    def BuildDistinctLineages(self, event_sequence): # same here, must go from tables
        # returns number of distinct lineages which corresponds to the times at the event sequence

        distinct_lineages = [0] * len(event_sequence)

        for iteration in range(len(event_sequence)):
            if iteration > 0:
                distinct_lineages[iteration] = distinct_lineages[iteration-1] - 1
            else:
                distinct_lineages[iteration] = -1

            for individual_tree in self.estimated_tree:
                if individual_tree.contains(event_sequence[iteration].node_id):
                    if event_sequence[iteration].tree_type == "coalescence":
                        number_of_children = len(individual_tree.get_node(event_sequence[iteration].node_id).fpointer)
                        distinct_lineages[iteration] = distinct_lineages[iteration] + number_of_children
                            #print('Added ', number_of_children, 'children')
                    if individual_tree.root == event_sequence[iteration].node_id:
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
        return self.es.tree_sequence[iteration]

    def EventProbability(self, coal_rate, iteration, distinct_lineages):

        if self.event_array[iteration].is_sample == 1:
            return 1
        else:
            probability = 1
            for i in range(0, 2): # Need to have a number of children
                probability = probability * coal_rate * math.comb(distinct_lineages[iteration] - i + 1, 2)
            return probability


    def GetEstimation(self, needed_timestamp): # this needs no real change, just redo the funcction

        print("GETTING ESTIMATION")
        #start_coal_rate=20 # will have to estimate
        wrapper = lambda coal_rate: - self.LLH_function(coal_rate=coal_rate, needed_timestamp=needed_timestamp)

        LLH_optimised = optimize.minimize_scalar(fun=wrapper, bounds = (0.0001, 1000), bracket = (0.0001, 10), method='Bounded', tol=1e-5)
        #LLH_optimised = optimize.minimize(fun=wrapper, x0=20, method='Nelder-Mead', tol=1e-2)

        # plotting

        x = [x /1000.0 for x in range(1, 1000, 1)]
        y = [LikelyhoodEstimationDismembered.LLH_function(self, coal_rate=i, needed_timestamp=needed_timestamp) for i in x]

        plt.plot(x, y)
        plt.show()



        LLH_result = LLH_optimised.fun
        LLH_point = LLH_optimised.x

        print(LLH_result, LLH_point)

        return [LLH_result, LLH_point]


LED = LikelyhoodEstimationDismembered(event_table_funct, event_table_neutral, sample_fraction_table)
a = LED.GetEstimation(1.2)
print(a)




