import math
import statistics

from scipy import optimize
from Simulator.VGsim._BirthDeath import BirthDeathModel
import matplotlib.pyplot as plt
from simulation import iterations, bRate, dRate, sRate, mRate, popModel, susceptible, lockdownModel, rndseed

import logging # a workaround to kill warnings
logging.captureWarnings(True)

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
        if event_table_funct != [] and event_table_neutral != []:
            event_table = event_table_funct[0] + event_table_neutral[0] # merges lists
        if event_table_funct == []:
            event_table = event_table_neutral[0]
        if event_table_neutral == []:
            event_table = event_table_funct[0]


        #format of data is [left timestamp, number of samples, number_of coals, fraction]
        self.bracket_data = []
        for bracket in sample_fraction_table:
            self.bracket_data.append([[], [], [], bracket[0], bracket[1]])



        self.timestamps = [sample_fraction_table[i][0] for i in range(len(sample_fraction_table))]

        # bracket_data is a list of brackets Each corresponds to a timeframe
        # each bracket is a list of 4 values
        # 0) list of times
        # 1) list of indicators if the event is sample
        # 2) list of indicators if the event is coalescence
        # 3) bracket timstamp
        # 4) bracket fraction


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

        self.distinct_lineages_array = [[] for _ in range(len(self.event_array))]
        for i in range(len(self.event_array)):
            for j in range(len(self.event_array[i])):
                self.distinct_lineages_array[i].append(0)

        current_lineages = 1
        for i in range(len(self.event_array)):
            for j in range(len(self.event_array[i])):
                if self.event_array[i][j].is_coal == 1:
                    if current_lineages == 0:
                        self.distinct_lineages_array[i][j] = current_lineages + 2
                    else:
                        self.distinct_lineages_array[i][j] = current_lineages + 1
                    current_lineages = self.distinct_lineages_array[i][j]
                if self.event_array[i][j].is_sample == 1:
                    self.distinct_lineages_array[i][j] = current_lineages - 1
                    current_lineages = self.distinct_lineages_array[i][j]

    # we do a preprocessing of values for LLH
    # LLH = -coal_rate*coal_rate_multiplier + sum_of_logs


        self.coal_rate_multipliers = [0 for _ in range(len(self.bracket_data))]
        self.number_of_coals = [0 for _ in range(len(self.bracket_data))]
        self.number_of_samples = [0 for _ in range(len(self.bracket_data))]


        for i in range(len(self.coal_rate_multipliers)):
            for j in range(1, len(self.event_array[i])):
                self.coal_rate_multipliers[i] = self.coal_rate_multipliers[i] + \
                - (self.event_array[i][j].time - self.event_array[i][j-1].time) * \
                    math.comb(self.distinct_lineages_array[i][j], 2)

        # sums_of_logs equals number of coalescent events

        for i in range(len(self.number_of_coals)):
            self.number_of_coals[i] = sum(self.bracket_data[i][2])
        for i in range(len(self.number_of_samples)):
            self.number_of_samples[i] = sum(self.bracket_data[i][1])


    def GetEstimationAnalytics(self):
        print("GETTING ESTIMATION")
        results = []
        for timestamp_num in range(len(self.timestamps)):
            if self.coal_rate_multipliers[timestamp_num] == 0:
                result = 0
            else:
                result = (-1) * self.number_of_coals[timestamp_num] / self.coal_rate_multipliers[timestamp_num]
            results.append(result)
        print("ESTIMATED")

        return results











