import math
import sys

import scipy
import numpy as np
import statistics

from scipy import optimize
from scipy.optimize import Bounds
from Simulator.VGsim._BirthDeath import BirthDeathModel
import matplotlib.pyplot as plt
from simulation import iterations, bRate, dRate, sRate, mRate, popModel, susceptible, lockdownModel, rndseed
from scipy.optimize import fsolve

import logging # a workaround to kill warnings
logging.captureWarnings(True)

class LikelyhoodEstimationDismembered:
    def __init__(self, event_table_funct=None, event_table_neutral=None, sample_fraction_table=None): # here we have not tree, but tables

        def takeTime(event):
            return event[0]

        self.roots_neutral = []
        for event_tree in event_table_neutral:
            event_tree.sort(key=takeTime)
            self.roots_neutral.append(event_tree[0])

        self.roots_funct = []
        for event_tree in event_table_funct:
            event_tree.sort(key=takeTime)
            self.roots_funct.append(event_tree[0])

        #if event_table_funct != []:
        #    for timestamp_num in range(1, len(event_table_funct)):
        #        event_table_funct[0] = event_table_funct[0] + event_table_funct[timestamp_num]
        #   event_table_funct = event_table_funct[0]
        #if event_table_neutral != []
        #    event_table_neutral = event_table_neutral[0]


        self.timestamps = [0.0] + [sample_fraction_table[i][0] for i in range(1, len(sample_fraction_table))]
        # TODO - after Shishkin makes a fix - delete 0.0
        self.number_of_timestamps = len(self.timestamps) - 1

        #format of data is [left timestamp, number of samples, number_of coals, fraction]

        self.bracket_data_neutral = [[] for _ in range(self.number_of_timestamps)]
        self.bracket_data_funct = [[] for _ in range(self.number_of_timestamps)]

        # bracket_data is a list of brackets Each corresponds to a timeframe
        # each bracket is a list of 4 values
        # 0) list of times
        # 1) list of indicators if the event is sample
        # 2) list of indicators if the event is coalescence

        for timestamp_num in range(self.number_of_timestamps):
            if timestamp_num < self.number_of_timestamps - 1:
                time_start = self.timestamps[timestamp_num]
                time_finish = self.timestamps[timestamp_num + 1]
            else:
                time_start = self.timestamps[timestamp_num]
                time_finish = sys.maxsize
            for event_tree_neutral in event_table_neutral:
                for event_neutral in event_tree_neutral:
                    if event_neutral[0] >= time_start and event_neutral[0] < time_finish:
                        self.bracket_data_neutral[timestamp_num].append(event_neutral)
            for event_tree_funct in event_table_funct:
                for event_funct in event_tree_funct:
                    if event_funct[0] >= time_start and event_funct[0] < time_finish:
                        self.bracket_data_funct[timestamp_num].append(event_funct)




        for bracket in self.bracket_data_neutral:
            bracket.sort(key=takeTime)
        for bracket in self.bracket_data_funct:
            bracket.sort(key=takeTime)

        self.distinct_lineages_array_neutral = [[] for _ in range(self.number_of_timestamps)]
        self.distinct_lineages_array_funct = [[] for _ in range(self.number_of_timestamps)]

        for timestamp_num in range(self.number_of_timestamps):
            for sample_num in range(len(self.bracket_data_neutral[timestamp_num])):
                self.distinct_lineages_array_neutral[timestamp_num].append(0)
            for sample_num in range(len(self.bracket_data_funct[timestamp_num])):
                self.distinct_lineages_array_funct[timestamp_num].append(0)


        current_lineages_neutral = 0
        current_lineages_funct = 0
        for timestamp_num in range(self.number_of_timestamps):
            for sample_num in range(len(self.bracket_data_neutral[timestamp_num])):
                if self.bracket_data_neutral[timestamp_num][sample_num] in self.roots_neutral:
                    if self.bracket_data_neutral[timestamp_num][sample_num][2] == 1: # is coal
                        self.distinct_lineages_array_neutral[timestamp_num][sample_num] = current_lineages_neutral + 2
                    else: # is sample
                        continue
                else:
                    if self.bracket_data_neutral[timestamp_num][sample_num][2] == 1: # is coal
                        self.distinct_lineages_array_neutral[timestamp_num][sample_num] = current_lineages_neutral + 1
                    else: # is sample
                        self.distinct_lineages_array_neutral[timestamp_num][sample_num] = current_lineages_neutral - 1

                current_lineages_neutral = self.distinct_lineages_array_neutral[timestamp_num][sample_num]

            for sample_num in range(len(self.bracket_data_funct[timestamp_num])):
                if self.bracket_data_funct[timestamp_num][sample_num] in self.roots_funct:
                    if self.bracket_data_funct[timestamp_num][sample_num][2] == 1:  # is coal
                        self.distinct_lineages_array_funct[timestamp_num][sample_num] = current_lineages_funct + 2
                    else:  # is sample
                        continue
                else:
                    if self.bracket_data_funct[timestamp_num][sample_num][2] == 1:  # is coal
                        self.distinct_lineages_array_funct[timestamp_num][sample_num] = current_lineages_funct + 1
                    else:  # is sample
                        self.distinct_lineages_array_funct[timestamp_num][sample_num] = current_lineages_funct - 1

                current_lineages_funct = self.distinct_lineages_array_funct[timestamp_num][sample_num]


        for timestamp_num in range(self.number_of_timestamps):
            for sample_num in range(len(self.bracket_data_neutral[timestamp_num])):
                if self.distinct_lineages_array_neutral[timestamp_num][sample_num] < 0:
                    raise(ValueError)
        for timestamp_num in range(self.number_of_timestamps):
            for sample_num in range(len(self.bracket_data_funct[timestamp_num])):
                if self.distinct_lineages_array_funct[timestamp_num][sample_num] < 0:
                    raise(ValueError)


        # we do a preprocessing of values for LLH
        # LLH = -coal_rate*coal_rate_multiplier + sum_of_logs

        # sums_of_logs equals number of coalescent events

        self.number_of_coals_neutral = [0 for _ in range(self.number_of_timestamps)]
        self.number_of_coals_funct = [0 for _ in range(self.number_of_timestamps)]
        self.number_of_samples_neutral = [0 for _ in range(self.number_of_timestamps)]
        self.number_of_samples_funct = [0 for _ in range(self.number_of_timestamps)]

        def SumCoals(dataEvents):
            number_of_coals = 0
            for i in range(len(dataEvents)):
                number_of_coals = number_of_coals + dataEvents[i][2]
            return number_of_coals

        def SumSamples(dataEvents):
            number_of_samples = 0
            for i in range(len(dataEvents)):
                number_of_samples = number_of_samples + dataEvents[i][1]
            return number_of_samples


        for timestamp_num in range(self.number_of_timestamps):
            if self.bracket_data_neutral[timestamp_num] == []:
                self.number_of_coals_neutral[timestamp_num] = 0
                self.number_of_samples_neutral[timestamp_num] = 0
            else:
                self.number_of_coals_neutral[timestamp_num] = SumCoals(self.bracket_data_neutral[timestamp_num])
                self.number_of_samples_neutral[timestamp_num] = SumSamples(self.bracket_data_neutral[timestamp_num])

            if self.bracket_data_funct[timestamp_num] == []:
                self.number_of_coals_funct[timestamp_num] = 0
                self.number_of_samples_funct[timestamp_num] = 0
            else:
                self.number_of_coals_funct[timestamp_num] = SumCoals(self.bracket_data_funct[timestamp_num])
                self.number_of_samples_funct[timestamp_num] = SumSamples(self.bracket_data_funct[timestamp_num])



    def GetEstimationConstantsGivenRho(self, rho): # returns an estimation of the s_i with respect to rho

        # here we define constants in the LLH_1 and LLH_2 as follows
        # the LLH results in the formula for neutral and funct cases
        # c_1 * \lambda - number_of_coals_neutral * log \lambda + c_2
        # c_3 * (\rho * \lambda) - number_of_coals_funct * log (\lambda * \rho) + c_4

        c1s = [0 for _ in range(self.number_of_timestamps)]
        c2s = [0 for _ in range(self.number_of_timestamps)]
        c3s = [0 for _ in range(self.number_of_timestamps)]
        c4s = [0 for _ in range(self.number_of_timestamps)]

        for timestamp_num in range(self.number_of_timestamps):
            for j in range(len(self.bracket_data_neutral[timestamp_num]) - 1):
                c1s[timestamp_num] = c1s[timestamp_num] + \
                                          (self.bracket_data_neutral[timestamp_num][j+1][0] -
                                           self.bracket_data_neutral[timestamp_num][j][0]) * \
                                          math.comb(self.distinct_lineages_array_neutral[timestamp_num][j], 2)


        for timestamp_num in range(self.number_of_timestamps):
            for j in range(len(self.bracket_data_neutral[timestamp_num]) - 1):
                if (self.distinct_lineages_array_neutral[timestamp_num][j] >= 2):
                    c2s[timestamp_num] = c2s[timestamp_num] + self.number_of_coals_neutral[timestamp_num]* \
                                           math.log(math.comb(self.distinct_lineages_array_neutral[timestamp_num][j], 2))
                else:
                    # a workaround since log of zero is not defined
                    continue


        for timestamp_num in range(self.number_of_timestamps):
            for j in range(len(self.bracket_data_funct[timestamp_num]) - 1):
                c3s[timestamp_num] = c3s[timestamp_num] + \
                                          (self.bracket_data_funct[timestamp_num][j+1][0] -
                                           self.bracket_data_funct[timestamp_num][j][0]) * \
                                          math.comb(self.distinct_lineages_array_funct[timestamp_num][j], 2)

        for timestamp_num in range(self.number_of_timestamps):
            for j in range(len(self.bracket_data_funct[timestamp_num]) - 1):
                c4s[timestamp_num] = c4s[timestamp_num] +\
                                          (self.bracket_data_funct[timestamp_num][j+1][0] -
                                           self.bracket_data_funct[timestamp_num][j][0]) * \
                                          math.comb(self.distinct_lineages_array_funct[timestamp_num][j], 2)

                if (self.distinct_lineages_array_funct[timestamp_num][j] >= 2):
                    c4s[timestamp_num] = c4s[timestamp_num] + \
                                     math.log(math.comb(self.distinct_lineages_array_funct[timestamp_num][j], 2))
                else:
                    # a workaround since log of zero is not defined
                    continue

        return c1s, c2s, c3s, c4s


    def GetLHHOptimumsOnSamples(self, rho):

        #returns estimations for lambdas

        c1s, c2s, c3s, c4s = self.GetEstimationConstantsGivenRho(rho)

        lambdas = [0 for _ in range(self.number_of_timestamps)]
        LLHOptimumResultsNoConstantTerm = [0 for _ in range(self.number_of_timestamps)]

        for timestamp_num in range(self.number_of_timestamps):
            if (self.number_of_coals_neutral[timestamp_num] + self.number_of_coals_funct[timestamp_num] == 0) or \
                (c1s[timestamp_num] + c3s[timestamp_num] == 0):
                lambdas[timestamp_num] = 0
                LLHOptimumResultsNoConstantTerm[timestamp_num] = 0
                # since we know nothing, it doesn't influence the LLH
            else:
                lambdas[timestamp_num] = (self.number_of_coals_neutral[timestamp_num] + self.number_of_coals_funct[timestamp_num]) / (c1s[timestamp_num] + rho * c3s[timestamp_num])
                LLHOptimumResultsNoConstantTerm[timestamp_num] = (c1s[timestamp_num] + rho * c3s[timestamp_num]) * lambdas[timestamp_num] - \
                                                                 (self.number_of_coals_neutral[timestamp_num]+self.number_of_coals_funct[timestamp_num])*math.log(lambdas[timestamp_num]) - \
                    self.number_of_coals_funct[timestamp_num]*rho

        return LLHOptimumResultsNoConstantTerm


    def GetLLHOptimumTotal(self, rho):
        LLHOptimumResultsNoConstantTerm = self.GetLHHOptimumsOnSamples(rho)
        result = sum(LLHOptimumResultsNoConstantTerm)
        # we use an addition, since we work with the logarithms

        return result

    def OptimiseLLH(self):
        overall_optimizer = lambda rho: - self.GetLLHOptimumTotal(rho)
        overall_result = scipy.optimize.minimize(fun=overall_optimizer, x0=np.array([3]), bounds=Bounds(0.001, 1000, keep_feasible=True)).x
        return overall_result

    def PlotLLH(self):
        results = [0 for _ in range(0, 10)]

        for i in range(len(results)):
            results[i] = self.GetLLHOptimumTotal(0.2*i+0.2)

        plt.plot(results)
        plt.show()
        return 0


