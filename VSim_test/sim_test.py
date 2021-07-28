from Simulator.VGsim._BirthDeath import BirthDeathModel
import matplotlib.pyplot as plt
from simulation import iterations, bRate, dRate, sRate, mRate, popModel, susceptible, lockdownModel, rndseed
import statistics
from random import randrange
from likelyhood_estimation_dismembered import LikelyhoodEstimationDismembered

import sys



results = []
sampled = []
for i in range(0, 2):
    rndseed = randrange(sys.maxsize)
    simulation = BirthDeathModel(iterations, bRate, dRate, sRate, mRate, populationModel=popModel,
                                 susceptible=susceptible, lockdownModel=lockdownModel, rndseed=rndseed)
    simulation.SimulatePopulation(iterations)
    simulation.GetGenealogy()
    tdm = simulation.gettdm() #get tdm object
    trees_funct, trees_neutral = tdm.Dismember() #перед получением таблиц, нужно разчленить дерево
            #получение таблиц
    event_table_funct, event_table_neutral = tdm.getEventTable() #[{time: [n_samples, n_coals]}]
    brackets = [_ for _ in range(1, 101, 1)]
    sample_fraction_table = tdm.getSampleFracTable(brackets)
    LED = LikelyhoodEstimationDismembered(event_table_funct, event_table_neutral, sample_fraction_table)
    result = LED.GetEstimationAnalytics()
    results.append(result)
    sampled.append(LED.number_of_samples)

fractions = []
for result_num in range(len(results[0])):
    if results[1][result_num] == 0 or sampled[1][result_num] == 0:
        fraction = 999999
    else:
        fraction = (results[0][result_num] / results[1][result_num]) * (sampled[0][result_num] / sampled[1][result_num])

    fractions.append(fraction)

    fractions = list(filter(lambda a: a != 0 and a!=999999, fractions))

overall_result = statistics.median(fractions)
print(fractions)
print(overall_result)



