from Simulator.VGsim._BirthDeath import BirthDeathModel
import matplotlib.pyplot as plt
from simulation import iterations, bRate, dRate, sRate, mRate, popModel, susceptible, lockdownModel, rndseed
import statistics
from random import randrange
from likelyhood_estimation_dismembered import LikelyhoodEstimationDismembered

import sys

results = []
sampled = []
brackets = [_ for _ in range(1, 101, 1)]

def Simulate(iterations, bRate, dRate, sRate, mRate, popModel,
                                 susceptible, lockdownModel, rndseed):

    simulation = BirthDeathModel(iterations, bRate, dRate, sRate, mRate, populationModel=popModel,
                                     susceptible=susceptible, lockdownModel=lockdownModel, rndseed=rndseed)
    simulation.SimulatePopulation(iterations)
    simulation.GetGenealogy()
    tdm = simulation.gettdm() #get tdm object
    trees_funct, trees_neutral = tdm.Dismember() #перед получением таблиц, нужно разчленить дерево
                #получение таблиц
    event_table_funct, event_table_neutral = tdm.getEventTable() #[{time: [n_samples, n_coals]}]
    sample_fraction_table = tdm.getSampleFracTable(brackets)
    LED = LikelyhoodEstimationDismembered(event_table_funct, event_table_neutral, sample_fraction_table)
    result = LED.GetEstimationAnalytics()

    return result, LED.number_of_samples

result1, sampled1 = Simulate(iterations, bRate, dRate, sRate, mRate, popModel,
                                 susceptible, lockdownModel, randrange(sys.maxsize))
result2, sampled2 = Simulate(iterations, bRate, dRate, sRate, mRate, popModel,
                                 susceptible, lockdownModel, randrange(sys.maxsize))

results.append(result1)
results.apend(result2)
sampled.append(sampled1)
sampled.append(sampled2)

fractions = []
for result_num in range(len(results[0])):
    if results[1][result_num] == 0 or sampled[1][result_num] == 0:
        fraction = 999999
    else:
        fraction = (results[0][result_num] / results[1][result_num]) * (sampled[0][result_num] / sampled[1][result_num])

    fractions.append(fraction)

fractions = list(filter(lambda a: a != 0 and a!=999999, fractions))
med = statistics.median(fractions)
fractions = sorted(fractions, key = lambda a: abs(a-med))
fractions = fractions[: len(fractions) // 2]

print(statistics.mean(fractions))

def EstimateRho(simulation, hap1, hap2):
    return 0