from Simulator.VGsim._BirthDeath import BirthDeathModel
import matplotlib.pyplot as plt
from simulation import iterations, bRate, dRate, sRate, mRate, popModel, susceptible, lockdownModel, rndseed
import statistics
from random import randrange

from likelyhood_estimation_dismembered_new import LikelyhoodEstimationDismembered

import sys

results = []
sampled = []
LLHs = []
brackets = [_*2 for _ in range(50)]

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
    rho = LED.GetLLHOptimumsTotal()

    return rho

rho = Simulate(iterations, bRate, dRate, sRate, mRate, popModel,
                                 susceptible, lockdownModel, 1302640546012562505)
print(rho)


