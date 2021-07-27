from Simulator.VGsim._BirthDeath import BirthDeathModel
import matplotlib.pyplot as plt
from simulation import iterations, bRate, dRate, sRate, mRate, popModel, susceptible, lockdownModel, rndseed
import statistics
from random import randrange
from likelyhood_estimation_dismembered import LikelyhoodEstimationDismembered

import sys

results = []
for i in range(0, 2):
    simulation = BirthDeathModel(iterations, bRate, dRate, sRate, mRate, populationModel=popModel,
                                 susceptible=susceptible, lockdownModel=lockdownModel, rndseed=randrange(sys.maxsize))
    simulation.SimulatePopulation(iterations)
    simulation.GetGenealogy()
    tdm = simulation.gettdm() #get tdm object
    trees_funct, trees_neutral = tdm.Dismember() #перед получением таблиц, нужно разчленить дерево
            #получение таблиц
    event_table_funct, event_table_neutral = tdm.getEventTable() #[{time: [n_samples, n_coals]}]
    brackets = [_/5 for _ in range(1, 51, 1)]
    sample_fraction_table = tdm.getSampleFracTable(brackets)
    LED = LikelyhoodEstimationDismembered(event_table_funct, event_table_neutral, sample_fraction_table)
    result = LED.GetEstimationAnalytics()
    results.append(result)

fractions = []
for result_num in range(len(results[0])):
    if results[1][result_num] == 0:
        fraction = 0
    else:
        fraction = results[0][result_num] / results[1][result_num]
        fractions.append(fraction)

print(fractions)
print(statistics.mean(fractions))


