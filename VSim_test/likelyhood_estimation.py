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

        self.simulation = simulation
        self.event_table_funct = event_table_funct
        self.event_table_neutral = event_table_neutral

        def TakeEventTime(event):
            return event[0]

        self.roots_neutral = []
        for event_tree in event_table_neutral:
            event_tree.sort(key=TakeEventTime)
            self.roots_neutral.append(event_tree[0])

        self.roots_funct = []
        for event_tree in event_table_funct:
            event_tree.sort(key=TakeEventTime)
            self.roots_funct.append(event_tree[0])

        max_time = 0
        for neutral_tree in event_table_neutral:
            for neutral_event in neutral_tree:
                if neutral_event[0] > max_time:
                    max_time = neutral_event[0]
        for funct_tree in event_table_funct:
            for funct_event in funct_tree:
                if funct_event[0] > max_time:
                    max_time = funct_event[0]


        self.timestamps = [_ for _ in np.linspace(0, max_time, number_of_brackets + 1, endpoint=True)]
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
            bracket.sort(key=TakeEventTime)
        for bracket in self.bracket_data_funct:
            bracket.sort(key=TakeEventTime)



        self.distinct_lineages_array_neutral = [[0 for _ in range(len(self.bracket_data_neutral[timestamp_num]))]
                                                for timestamp_num in range(self.number_of_brackets)]
        self.distinct_lineages_array_funct = [[0 for _ in range(len(self.bracket_data_funct[timestamp_num]))]
                                                for timestamp_num in range(self.number_of_brackets)]

        current_lineages_neutral = 0
        current_lineages_funct = 0


        for bracket_num in range(self.number_of_brackets):

            for sample_num in range(len(self.bracket_data_neutral[bracket_num])):
                if len(self.bracket_data_neutral[bracket_num]) > 0:
                    # if there are many roots - next line adds much running time
                    if self.bracket_data_neutral[bracket_num][sample_num] in self.roots_neutral: # is root
                        if self.bracket_data_neutral[bracket_num][sample_num][2] == 1: # is coal
                            self.distinct_lineages_array_neutral[bracket_num][sample_num] = current_lineages_neutral + 2
                        else: # is sample
                            self.distinct_lineages_array_neutral[bracket_num] [sample_num] = current_lineages_neutral
                    else: # is not root
                        if self.bracket_data_neutral[bracket_num][sample_num][2] == 1: # is coal
                            self.distinct_lineages_array_neutral[bracket_num][sample_num] = current_lineages_neutral + 1
                        else: # is sample
                            self.distinct_lineages_array_neutral[bracket_num][sample_num] = current_lineages_neutral - 1
                    current_lineages_neutral = self.distinct_lineages_array_neutral[bracket_num][sample_num]
                else:
                    continue

            for sample_num in range(len(self.bracket_data_funct[bracket_num])):
                if len(self.bracket_data_funct[bracket_num]) > 0:
                    if self.bracket_data_funct[bracket_num][sample_num] in self.roots_funct: # is root
                        if self.bracket_data_funct[bracket_num][sample_num][2] == 1: # is coal
                            self.distinct_lineages_array_funct[bracket_num][sample_num] = current_lineages_funct + 2
                        else: # is sample
                            self.distinct_lineages_array_funct[bracket_num][sample_num] = current_lineages_funct
                    else: # is not root
                        if self.bracket_data_funct[bracket_num][sample_num][2] == 1: # is coal
                            self.distinct_lineages_array_funct[bracket_num][sample_num] =  current_lineages_funct + 1
                        else: # is sample
                            self.distinct_lineages_array_funct[bracket_num] [sample_num] =  current_lineages_funct - 1
                    current_lineages_funct = self.distinct_lineages_array_funct[bracket_num][sample_num]
                else:
                    continue


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
        for timestamp_num in range(self.number_of_brackets):
            if len(self.bracket_data_neutral[timestamp_num]) > 0:
                if timestamp_num == 0:
                    for j in range(1, len(self.bracket_data_neutral[timestamp_num])):
                        c1s[timestamp_num] += (self.bracket_data_neutral[timestamp_num][j][0] -
                                               self.bracket_data_neutral[timestamp_num][j - 1][0]) * \
                                              math.comb(self.distinct_lineages_array_neutral[timestamp_num][j - 1], 2)
                else:
                    c1s[timestamp_num] += (self.bracket_data_neutral[timestamp_num][0][0] -
                                               current_neutral_time) * \
                                              math.comb(current_neutral_lineages, 2)
                    for j in range(1, len(self.bracket_data_neutral[timestamp_num])):
                        c1s[timestamp_num] += (self.bracket_data_neutral[timestamp_num][j][0] -
                                               self.bracket_data_neutral[timestamp_num][j - 1][0]) * \
                                              math.comb(self.distinct_lineages_array_neutral[timestamp_num][j - 1], 2)

                current_neutral_time = self.bracket_data_neutral[timestamp_num][-1][0]
                current_neutral_lineages = self.distinct_lineages_array_neutral[timestamp_num][-1]
            else:
                continue

        current_funct_lineages = 0
        current_funct_time = 0
        for timestamp_num in range(self.number_of_brackets):
            if len(self.bracket_data_funct[timestamp_num]) > 0:
                if timestamp_num == 0:
                    for j in range(1, len(self.bracket_data_funct[timestamp_num])):
                        c3s[timestamp_num] += (self.bracket_data_funct[timestamp_num][j][0] -
                                               self.bracket_data_funct[timestamp_num][j - 1][0]) * \
                                              math.comb(self.distinct_lineages_array_funct[timestamp_num][j - 1], 2)
                else:
                    c3s[timestamp_num] += (self.bracket_data_funct[timestamp_num][0][0] -
                                               current_funct_time) * \
                                              math.comb(current_funct_lineages, 2)
                    for j in range(1, len(self.bracket_data_funct[timestamp_num])):
                        c3s[timestamp_num] += (self.bracket_data_funct[timestamp_num][j][0] -
                                               self.bracket_data_funct[timestamp_num][j - 1][0]) * \
                                              math.comb(self.distinct_lineages_array_funct[timestamp_num][j - 1], 2)

                current_funct_time = self.bracket_data_funct[timestamp_num][-1][0]
                current_funct_lineages = self.distinct_lineages_array_funct[timestamp_num][-1]
            else:
                continue

        return c1s, c3s


    def GetLLHOptimumTotal(self, rho):

        #returns estimations for lambdas and for the sum

        c1s, c3s = self.GetEstimationConstants()

        self.lambdas = [0 for _ in range(self.number_of_brackets)]
        LLHOptimumResultsNoConstantTerm = [0 for _ in range(self.number_of_brackets)]
        self.estimated_infected_ratio = [0 for _ in range(self.number_of_brackets)]
        self.true_infected_ratio = [0 for _ in range(self.number_of_brackets)]
        self.average = [0 for _ in range(self.number_of_brackets)]


        self.hd = self.simulation.log_dynamics(step=self.number_of_brackets, output_file=False)
        #S_len = len(hd['P0']['H0'])
        self.hd_a = self.hd['P0']['H0']
        self.hd_g = self.hd['P0']['H3']

        for i in range(len(self.hd_g)):
            if self.hd_g[i] != 0:
                self.true_infected_ratio[i] = self.hd_a[i] / self.hd_g[i]

        #self.ratio = [self.estimated_infected_ratio[i] / self.true_infected_ratio[i] for i in range(self.number_of_brackets)]
        #average = sum(self.ratio) / len(self.ratio)
        #print('Average is', average)
        #ld = self.simulation.LogDynamics(2*self.number_of_brackets)[1::2]

        #for bracket_num in range(self.number_of_brackets):
        #    if hd[bracket_num][3] != 0:
        #        true_infected_ratio[bracket_num] = hd[bracket_num][0]/hd[bracket_num][3]
        # here 0 means that we cannot get any info

        for timestamp_num in range(self.number_of_brackets):
            if (self.number_of_coals_neutral[timestamp_num] == 0) or (self.number_of_coals_funct[timestamp_num] == 0) or \
                (self.number_of_samples_neutral[timestamp_num] == 0) or (self.number_of_samples_funct[timestamp_num] == 0)\
                    or (c1s[timestamp_num] + c3s[timestamp_num] == 0) or (self.true_infected_ratio == 0)\
                    or (self.estimated_infected_ratio == 0):
                # will have to remove true in the final version
                LLHOptimumResultsNoConstantTerm[timestamp_num] = 0
                # since we know nothing, it doesn't influence the LLH
            else:
                self.estimated_infected_ratio[timestamp_num] = self.number_of_samples_neutral[timestamp_num] / self.number_of_samples_funct[timestamp_num]
                self.lambdas[timestamp_num] = (self.number_of_coals_neutral[timestamp_num] + self.number_of_coals_funct[timestamp_num]) / \
                                         (c1s[timestamp_num] + c3s[timestamp_num] * self.estimated_infected_ratio[timestamp_num] * rho)
                #experimental - we eliminate a constant term (c1s[timestamp_num] + c3s[timestamp_num] * estimated_infected_ratio[timestamp_num] * rho) * lambdas[timestamp_num]
                LLHOptimumResultsNoConstantTerm[timestamp_num] = - (self.number_of_coals_neutral[timestamp_num] + self.number_of_coals_funct[timestamp_num]) * math.log(self.lambdas[timestamp_num]) - \
                    self.number_of_coals_funct[timestamp_num] * math.log(rho)

        result = sum(LLHOptimumResultsNoConstantTerm)
        # we use an addition, since we work with the logarithms

        return result

    def OptimiseLLH(self):
        overall_optimizer = lambda rho: self.GetLLHOptimumTotal(rho)
        optimum = scipy.optimize.minimize_scalar(fun=overall_optimizer, bracket=(0.01, 10), bounds=(0.001, 1000000), method='Bounded', tol=0.0001)

        #self.lambda_div_si = [0 for _ in range(self.number_of_brackets)]
        #for timestamp_num in range(len(self.hd['P0']['S0'])):
        #    if self.hd['P0']['H0'][timestamp_num] > 0:
        #        self.lambda_div_si[timestamp_num] = self.lambdas[timestamp_num] / (self.hd['P0']['S0'][timestamp_num] / self.hd['P0']['H0'][timestamp_num])

        #plt.plot(self.lambda_div_si)

        plt.plot(self.true_infected_ratio, color='red')
        #plt.plot(self.estimated_infected_ratio, color='blue')

        #for i in range(self.number_of_brackets):
        #    if self.estimated_infected_ratio[i] != 0:
        #        self.average[i] = self.true_infected_ratio[i] / self.estimated_infected_ratio[i]

        #plt.plot(self.average, color='green')
        #plt.plot(np.ones(len(self.average)), color='black')

        plt.show()


        return optimum

    def PlotLLH(self):
        results = [0 for _ in range(0, 40)]
        # we need a minimum here

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


