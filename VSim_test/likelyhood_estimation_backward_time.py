import math
import sys

import scipy
import numpy as np
import scipy.stats
from scipy.stats.distributions import chi2
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
    def __init__(self, event_table_funct=None, event_table_neutral=None, number_of_brackets=None, simulation=None): # here we have not tree, but tables


        for event_tree_num in range(len(event_table_neutral)): # ugly way of time reversal
            for event_num in range(len(event_table_neutral[event_tree_num])):
                time =  event_table_neutral[event_tree_num][event_num][0]
                sample = event_table_neutral[event_tree_num][event_num][1]
                coal = event_table_neutral[event_tree_num][event_num][2]
                event_table_neutral[event_tree_num][event_num] = [-time, sample, coal]


        for event_tree_num in range(len(event_table_funct)): # ugly way of time reversal
            for event_num in range(len(event_table_funct[event_tree_num])):
                time =  event_table_funct[event_tree_num][event_num][0]
                sample = event_table_funct[event_tree_num][event_num][1]
                coal = event_table_funct[event_tree_num][event_num][2]
                event_table_funct[event_tree_num][event_num] = [-time, sample, coal]

        self.simulation = simulation

        def TakeEventTime(event):
            return event[0]

        self.roots_neutral = []
        for event_tree in event_table_neutral:
            event_tree.sort(key=TakeEventTime, reverse=False)
            self.roots_neutral.append(event_tree[-1])

        self.roots_funct = []
        for event_tree in event_table_funct:
            event_tree.sort(key=TakeEventTime, reverse=False)
            self.roots_funct.append(event_tree[-1])

        #if event_table_funct != []:
        #    for timestamp_num in range(1, len(event_table_funct)):
        #        event_table_funct[0] = event_table_funct[0] + event_table_funct[timestamp_num]
        #   event_table_funct = event_table_funct[0]
        #if event_table_neutral != []
        #    event_table_neutral = event_table_neutral[0]

        min_time = 0
        for neutral_tree in event_table_neutral:
            for neutral_event in neutral_tree:
                if neutral_event[0] < min_time:
                    min_time = neutral_event[0]
        for funct_tree in event_table_funct:
            for funct_event in funct_tree:
                if funct_event[0] < min_time:
                    min_time = funct_event[0]


        self.timestamps = [_ for _ in np.linspace(min_time, 0, number_of_brackets + 1, endpoint=True)]
        self.number_of_brackets = len(self.timestamps) - 1

        #format of data is [left timestamp, number of samples, number_of coals, fraction]

        self.bracket_data_neutral = [[] for _ in range(self.number_of_brackets)]
        self.bracket_data_funct = [[] for _ in range(self.number_of_brackets)]


        # bracket_data is a list of brackets Each corresponds to a timeframe
        # each bracket is a list of 4 values
        # 0) list of times
        # 1) list of indicators if the event is sample
        # 2) list of indicators if the event is coalescence

        def MapTimeToBracket(time, brackets, bracket_start_num, bracket_finish_num):
            test_bracket_num = (bracket_start_num + bracket_finish_num) // 2
            if (time < brackets[test_bracket_num][0]):
                return MapTimeToBracket(time, brackets, bracket_start_num, test_bracket_num - 1)
            elif (time > brackets[test_bracket_num][1]):
                return MapTimeToBracket(time, brackets, test_bracket_num + 1, bracket_finish_num)
            else:
                return test_bracket_num


        starts_and_finishes = []
        for bracket_num in range(self.number_of_brackets):
            time_start = self.timestamps[bracket_num]
            time_finish = self.timestamps[bracket_num + 1]
            starts_and_finishes.append([time_start, time_finish])


        for event_tree_neutral in event_table_neutral:
            for event_neutral in event_tree_neutral:
                bracket_num = MapTimeToBracket(event_neutral[0], starts_and_finishes, 0, self.number_of_brackets)
                self.bracket_data_neutral[bracket_num].append(event_neutral)
        for event_tree_funct in event_table_funct:
            for event_funct in event_tree_funct:
                bracket_num = MapTimeToBracket(event_funct[0], starts_and_finishes, 0, self.number_of_brackets)
                self.bracket_data_funct[bracket_num].append(event_funct)

        for bracket in self.bracket_data_neutral:
            bracket.sort(key=TakeEventTime, reverse=False)
        for bracket in self.bracket_data_funct:
            bracket.sort(key=TakeEventTime, reverse=False)

        self.distinct_lineages_array_neutral = [[0 for _ in range(len(self.bracket_data_neutral[timestamp_num]))]
                                                for timestamp_num in range(self.number_of_brackets)]
        self.distinct_lineages_array_funct = [[0 for _ in range(len(self.bracket_data_funct[timestamp_num]))]
                                                for timestamp_num in range(self.number_of_brackets)]

        current_lineages_neutral = 0
        current_lineages_funct = 0

        root_coals = 0
        non_root_coals = 0
        root_samples = 0
        non_root_samples = 0

        for bracket_num in range(self.number_of_brackets-1, -1, -1):
            for sample_num in range(len(self.bracket_data_neutral[bracket_num])-1, -1, -1):
                if len(self.bracket_data_neutral[bracket_num]) > 0:
                    # if there are many roots - next line adds much running time
                    if self.bracket_data_neutral[bracket_num][sample_num] in self.roots_neutral: # is root
                        if self.bracket_data_neutral[bracket_num][sample_num][2] == 1: # is coal
                            self.distinct_lineages_array_neutral[bracket_num][sample_num] =  current_lineages_neutral + 2
                        else: # is sample
                            self.distinct_lineages_array_neutral[bracket_num][sample_num] =  current_lineages_neutral
                    else: # is not root
                        if self.bracket_data_neutral[bracket_num][sample_num][2] == 1: # is coal
                            self.distinct_lineages_array_neutral[bracket_num][sample_num] = current_lineages_neutral + 1
                        else: # is sample
                            self.distinct_lineages_array_neutral[bracket_num][sample_num] = current_lineages_neutral - 1
                    current_lineages_neutral = self.distinct_lineages_array_neutral[bracket_num][sample_num]
                else:
                    continue

            for sample_num in range(len(self.bracket_data_funct[bracket_num])-1, -1, -1):
                if len(self.bracket_data_funct[bracket_num]) > 0:
                    if self.bracket_data_funct[bracket_num][sample_num] in self.roots_funct: # is root
                        if self.bracket_data_funct[bracket_num][sample_num][2] == 1: # is coal
                            self.distinct_lineages_array_funct[bracket_num][sample_num] = current_lineages_funct + 2
                            root_coals +=1
                        else: # is sample
                            self.distinct_lineages_array_funct[bracket_num][sample_num] = current_lineages_funct
                            root_samples +=1
                    else: # is not root
                        if self.bracket_data_funct[bracket_num][sample_num][2] == 1: # is coal
                            self.distinct_lineages_array_funct[bracket_num][sample_num] = current_lineages_funct + 1
                            non_root_coals +=1
                        else: # is sample
                            self.distinct_lineages_array_funct[bracket_num][sample_num] = current_lineages_funct - 1
                            non_root_samples +=1
                    current_lineages_funct = self.distinct_lineages_array_funct[bracket_num][sample_num]
                else:
                    continue

        print("Root coals", root_coals)
        print("Root samples", root_samples)
        print("Non Root coals", non_root_coals)
        print("Non Root samples", non_root_samples)


        for bracket_num in range(self.number_of_brackets):
            for sample_num in range(len(self.bracket_data_neutral[bracket_num])):
                if self.distinct_lineages_array_neutral[bracket_num][sample_num] < 0:
                    raise(ValueError)
        for bracket_num in range(self.number_of_brackets):
            for sample_num in range(len(self.bracket_data_funct[bracket_num])):
                if self.distinct_lineages_array_funct[bracket_num][sample_num] < 0:
                    raise(ValueError)



        # we do a preprocessing of values for LLH
        # LLH = -coal_rate*coal_rate_multiplier + sum_of_logs

        # sums_of_logs equals number of coalescent events

        self.number_of_coals_neutral = [0 for _ in range(self.number_of_brackets)]
        self.number_of_coals_funct = [0 for _ in range(self.number_of_brackets)]
        self.number_of_samples_neutral = [0 for _ in range(self.number_of_brackets)]
        self.number_of_samples_funct = [0 for _ in range(self.number_of_brackets)]

        def SumCoals(dataEvents):
            number_of_coals = 0
            for i in range(len(dataEvents)):
                number_of_coals += dataEvents[i][2]
            return number_of_coals

        def SumSamples(dataEvents):
            number_of_samples = 0
            for i in range(len(dataEvents)):
                number_of_samples += dataEvents[i][1]
            return number_of_samples


        for bracket_num in range(self.number_of_brackets):
            if self.bracket_data_neutral[bracket_num] == []:
                self.number_of_coals_neutral[bracket_num] = 0
                self.number_of_samples_neutral[bracket_num] = 0
            else:
                self.number_of_coals_neutral[bracket_num] = SumCoals(self.bracket_data_neutral[bracket_num])
                self.number_of_samples_neutral[bracket_num] = SumSamples(self.bracket_data_neutral[bracket_num])

            if self.bracket_data_funct[bracket_num] == []:
                self.number_of_coals_funct[bracket_num] = 0
                self.number_of_samples_funct[bracket_num] = 0
            else:
                self.number_of_coals_funct[bracket_num] = SumCoals(self.bracket_data_funct[bracket_num])
                self.number_of_samples_funct[bracket_num] = SumSamples(self.bracket_data_funct[bracket_num])

        self.number_of_neutral_vertices = sum(self.number_of_coals_neutral) + sum(self.number_of_samples_neutral)
        self.number_of_funct_vertices = sum(self.number_of_coals_funct) + sum(self.number_of_samples_funct)
        self.number_of_overall_vertices = sum(self.number_of_coals_neutral) + sum(self.number_of_samples_neutral) + \
                                          sum(self.number_of_coals_funct) + sum(self.number_of_samples_funct)

        print("There are", self.number_of_neutral_vertices, "vertices with a neutral haplotype out of",
              self.number_of_overall_vertices)
        print("Overall", sum(self.number_of_samples_neutral) + sum(self.number_of_samples_funct), "vertices were sampled out of", self.number_of_overall_vertices)




    def GetEstimationConstants(self): # returns an estimation of the s_i with respect to rho

        # here we define constants in the LLH_1 and LLH_2 as follows
        # the LLH results in the formula for neutral and funct cases
        # c_1 * \lambda - number_of_coals_neutral * log \lambda + c_2
        # c_3 * (\rho * \lambda) - number_of_coals_funct * log (\lambda * \rho) + c_4

        c1s = [0 for _ in range(self.number_of_brackets)]
        c3s = [0 for _ in range(self.number_of_brackets)]

        current_neutral_lineages = 0
        current_neutral_time = 0
        for timestamp_num in range(self.number_of_brackets-1, -1, -1):
            if len(self.bracket_data_neutral[timestamp_num]) > 0:
                if timestamp_num == self.number_of_brackets-1:
                    for j in range(len(self.bracket_data_neutral[timestamp_num]) - 1, 0, -1):
                        c1s[timestamp_num] += (self.bracket_data_neutral[timestamp_num][j-1][0] -
                                               self.bracket_data_neutral[timestamp_num][j][0]) * \
                                              math.comb(self.distinct_lineages_array_neutral[timestamp_num][j], 2)
                else:
                    c1s[timestamp_num] += (self.bracket_data_neutral[timestamp_num][-1][0] -
                                               current_neutral_time) * \
                                              math.comb(current_neutral_lineages, 2)
                    for j in range(len(self.bracket_data_neutral[timestamp_num]) - 2, 0, -1):
                        c1s[timestamp_num] += (self.bracket_data_neutral[timestamp_num][j-1][0] -
                                               self.bracket_data_neutral[timestamp_num][j][0]) * \
                                              math.comb(self.distinct_lineages_array_neutral[timestamp_num][j], 2)

                current_neutral_time = self.bracket_data_neutral[timestamp_num][0][0]
                current_neutral_lineages = self.distinct_lineages_array_neutral[timestamp_num][0]
            else:
                continue



        current_funct_lineages = 0
        current_funct_time = 0
        for timestamp_num in range(self.number_of_brackets-1, -1, -1):
            if len(self.bracket_data_funct[timestamp_num]) > 0:
                if timestamp_num == self.number_of_brackets-1:
                    for j in range(len(self.bracket_data_funct[timestamp_num]) - 1, 0, -1):
                        c3s[timestamp_num] += (self.bracket_data_funct[timestamp_num][j-1][0] -
                                               self.bracket_data_funct[timestamp_num][j][0]) * \
                                              math.comb(self.distinct_lineages_array_funct[timestamp_num][j], 2)
                else:
                    c3s[timestamp_num] += (self.bracket_data_funct[timestamp_num][-1][0] -
                                               current_funct_time) * \
                                              math.comb(current_funct_lineages, 2)
                    for j in range(len(self.bracket_data_funct[timestamp_num]) - 2, 0, -1):
                        c3s[timestamp_num] += (self.bracket_data_funct[timestamp_num][j-1][0] -
                                               self.bracket_data_funct[timestamp_num][j][0]) * \
                                              math.comb(self.distinct_lineages_array_funct[timestamp_num][j], 2)

                current_funct_time = self.bracket_data_funct[timestamp_num][0][0]
                current_funct_lineages = self.distinct_lineages_array_funct[timestamp_num][0]
            else:
                continue

        return c1s, c3s


    def GetLLHOptimumTotal(self, rho):

        #returns estimations for lambdas and for the sum

        c1s, c3s = self.GetEstimationConstants()

        lambdas = [0 for _ in range(self.number_of_brackets)]
        LLHOptimumResultsNoConstantTerm = [0 for _ in range(self.number_of_brackets)]
        estimated_infected_ratio = [0 for _ in range(self.number_of_brackets)]
        #true_infected_ratio = [0 for _ in range(self.number_of_brackets)]
        #hd = self.simulation.GetHaplotypeDynamics(2*self.number_of_brackets)[1::2] # we don't take a 0.0 timestamp
        #ld = self.simulation.LogDynamics(2*self.number_of_brackets)[1::2]

        #for bracket_num in range(self.number_of_brackets):
        #    if hd[bracket_num][3] != 0:
        #        true_infected_ratio[bracket_num] = hd[bracket_num][0]/hd[bracket_num][3]
        # here 0 means that we cannot get any info

        for timestamp_num in range(self.number_of_brackets-1, -1, -1):
            if (self.number_of_coals_neutral[timestamp_num] == 0) or (self.number_of_coals_funct[timestamp_num] == 0) or \
                (self.number_of_samples_neutral[timestamp_num] == 0) or (self.number_of_samples_funct[timestamp_num] == 0)\
                    or (c1s[timestamp_num] + c3s[timestamp_num] == 0):
                LLHOptimumResultsNoConstantTerm[timestamp_num] = 0
                # since we know nothing, it doesn't influence the LLH
            else:
                estimated_infected_ratio[timestamp_num] = self.number_of_samples_neutral[timestamp_num] / self.number_of_samples_funct[timestamp_num]
                lambdas[timestamp_num] = -(self.number_of_coals_neutral[timestamp_num] + self.number_of_coals_funct[timestamp_num]) / \
                                         (c1s[timestamp_num] + c3s[timestamp_num] * estimated_infected_ratio[timestamp_num] * rho)

                #experimental - we eliminate a constant term (c1s[timestamp_num] + c3s[timestamp_num] * estimated_infected_ratio[timestamp_num] * rho) * lambdas[timestamp_num]
                LLHOptimumResultsNoConstantTerm[timestamp_num] = - (self.number_of_coals_neutral[timestamp_num]
                                                                    + self.number_of_coals_funct[timestamp_num]) * \
                                                                 math.log(lambdas[timestamp_num]) - \
                    self.number_of_coals_funct[timestamp_num] * math.log(rho)
        #print("Estimated infected ratio", estimated_infected_ratio)

        #print("Estimated:", estimated_infected_ratio)
        #print("True:", true_infected_ratio)
        result = sum(LLHOptimumResultsNoConstantTerm)
        # we use an addition, since we work with the logarithms

        #result = LLHOptimumResultsNoConstantTerm[2]
        return result

    def OptimiseLLH(self):
        overall_optimizer = lambda rho: self.GetLLHOptimumTotal(rho)
        optimum = scipy.optimize.minimize_scalar(fun=overall_optimizer, bracket=(0.001, 10), bounds=(0.001, 1000000), method='Bounded')
        return optimum

    def PlotLLH(self):
        results = [0 for _ in range(0, 40)]
        # we need a minimum gere

        for i in range(len(results)):
            results[i] = - self.GetLLHOptimumTotal(0.01*i+0.9)

        plt.plot(results)
        plt.show()
        return 0

    def ConductLikelyhoodRatioTest(self, resulting_LLH, hypothesis_value):

        lr = 2 * (hypothesis_value - resulting_LLH)

        p = chi2.sf(lr, 0)

        if p > 0.9772:
            print("Likelyhood ratio test has passed")
        else:
            print("WARNING, likelyhood ratio test has failed")

        return chi2


