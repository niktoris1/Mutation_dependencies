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
brackets = [0.3, 0.6, 0.9, 1.2, 1.5, 1.8, 2.1, 2.4, 2.7]
sample_fraction_table = tdm.getSampleFracTable(brackets)


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
        event_table = event_table_funct[0] + event_table_neutral[0] # merges lists

        #format of data is [left timestamp, number of samples, number_of coals, fraction]
        self.bracket_data = []
        for bracket in sample_fraction_table:
            self.bracket_data.append([[], [], [], bracket[0], bracket[1]])


        self.timestamps = [self.bracket_data[i][4] for i in range(len(self.bracket_data))]

        # bracket_data is a list of brackets Each corresponds to a timeframe
        # each bracket is a list of 4 values
        # 1) list of times
        # 2) list of indicators if the event is sample
        # 3) list of indicators if the event is coalescence
        # 4) bracket timstamp
        # 5) bracket fraction


        for i in range(len(self.bracket_data)):
            if i < len(self.bracket_data) - 1:
                time_start = self.bracket_data[i][3]
                time_finish = self.bracket_data[i+1][3]
            else:
                time_start = self.bracket_data[i][3]
                time_finish = 999
            for event in event_table:
                if event[0] >= time_start and event[0] < time_finish:
                    self.bracket_data[i][0].append(event[0])
                    self.bracket_data[i][1].append(event[1])
                    self.bracket_data[i][2].append(event[2])


        self.event_array = [[] for _ in range(len(self.bracket_data))]
        for bracket_num in range(len(self.bracket_data)):
            current_bracket = self.bracket_data[bracket_num]
            for j in range(len(current_bracket[0])):
                event = self.Event(time = current_bracket[0][j], is_sample=current_bracket[1][j],
                    is_coal = current_bracket[2][j], timestamp=current_bracket[3],
                    fraction=current_bracket[4])

                self.event_array[bracket_num].append(event)

            def takeTime(elem):
                return elem.time

            self.event_array[bracket_num].sort(key=takeTime)

        self.number_of_events = len(self.event_array)

    def LLH_function(self, coal_rate):

        # returns LLH function for all time brackets

        # Event(time, is_sample, is_coal, timestamp, fraction)



        distinct_lineages_array = [[] for _ in range(len(self.event_array)) ]
        event_probability_array = [[] for _ in range(len(self.event_array)) ]

        for i in range(len(self.event_array)):
            for j in range(len(self.event_array[i])):
                distinct_lineages_array[i].append(0)
                event_probability_array[i].append(0)



        for i in range(len(self.event_array)):
            for j in range(len(self.event_array[i])):
                if i == 0 and j == 0:
                    distinct_lineages_array[i][j] = 0
                elif i > 0 and j == 0:
                    distinct_lineages_array[i][j] = distinct_lineages_array[i - 1][len(distinct_lineages_array[i - 1]) - 1]
                elif j > 0:
                    distinct_lineages_array[i][j] = distinct_lineages_array[i][j - 1]

                if self.event_array[i][j].is_coal == 1:
                    distinct_lineages_array[i][j] = distinct_lineages_array[i][j] + 1
                if self.event_array[i][j].is_sample == 1:
                    distinct_lineages_array[i][j] = distinct_lineages_array[i][j] - 1

        for i in range(len(self.event_array)):
            for j in range(len(self.event_array[i])):
                if self.event_array[i][j].is_sample == 1:
                    event_probability_array[i][j] = 1
                else:
                    probability = 1
                    for iter in range(0, 1): # TODO - update for a bigger number of children
                        probability = probability * coal_rate * math.comb(distinct_lineages_array[i][j] - iter + 1, 2)
                    event_probability_array[i][j] = probability


        for i in range(self.number_of_events):
            event_probability_array[i] = LikelyhoodEstimationDismembered.EventProbability(self, coal_rate, i, distinct_lineages_array)

        for timestamp_number in range(len(self.timestamps)):
            LLH_value_by_timestamp = [[0] for _ in range(len(self.timestamps))]


        LLH_values = [[0] for _ in range(self.number_of_events)]
        addition = [[0] for _ in range(self.number_of_events)]

        for i in range(1, self.number_of_events):
            addition[i] = (- coal_rate * (self.event_array[i].time - self.event_array[i-1].time) *
               math.comb(distinct_lineages_array[i], 2)) + \
                          math.log(event_probability_array[i])

            LLH_values[i] = LLH_values[i - 1] + addition[i]

        return LLH_values[-1]

    def EventProbability(self, coal_rate, iteration, distinct_lineages_array):

        if self.event_array[iteration].is_coal == 0:
            return 1
        else:
            probability = 1
            for i in range(0, 1): # TODO - update for a bigger number of children
                probability = probability * coal_rate * math.comb(distinct_lineages_array[iteration] - i + 1, 2)
            return probability


    def GetEstimation(self, needed_timestamp): # this needs no real change, just redo the funcction

        print("GETTING ESTIMATION")
        #start_coal_rate=20 # will have to estimate
        wrapper = lambda coal_rate: - self.LLH_function(coal_rate=coal_rate)

        LLH_optimised = optimize.minimize_scalar(fun=wrapper, bounds = (0.0001, 1000), bracket = (0.0001, 10), method='Bounded', tol=1e-5)
        #LLH_optimised = optimize.minimize(fun=wrapper, x0=20, method='Nelder-Mead', tol=1e-2)

        # plotting

        x = [x /10.0 for x in range(1, 1000, 1)]
        y = [LikelyhoodEstimationDismembered.LLH_function(self, coal_rate=i, needed_timestamp=needed_timestamp) for i in x]

        plt.plot(x, y)
        plt.show()



        LLH_result = LLH_optimised.fun
        LLH_point = LLH_optimised.x

        print(LLH_result, LLH_point)

        return [LLH_result, LLH_point]


LED = LikelyhoodEstimationDismembered(event_table_funct, event_table_neutral, sample_fraction_table)
a = LED.GetEstimation(1.5)
print(a)




