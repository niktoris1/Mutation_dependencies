import random

from Simulator.VGsim._BirthDeath import BirthDeathModel
import matplotlib.pyplot as plt
from simulation import iterations, bRate, dRate, sRate, mRate, popModel, susceptible, lockdownModel, rndseed
import statistics
from random import randrange
import numpy as np

from likelyhood_estimation_dismembered_new import LikelyhoodEstimationDismembered

import sys



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

    max_time = 0
    for neutral_tree in event_table_neutral:
        for neutral_event in neutral_tree:
            if neutral_event[0] > max_time:
                max_time = neutral_event[0]
    for funct_tree in event_table_funct:
        for funct_event in funct_tree:
            if funct_event[0] > max_time:
                max_time = funct_event[0]

    freq = 50 # how frequent brackets are
    brackets = [_ for _ in np.linspace(max_time/freq, max_time, freq, endpoint=True)]
    sample_fraction_table = tdm.getSampleFracTable(brackets)
    LED = LikelyhoodEstimationDismembered(event_table_funct, event_table_neutral, sample_fraction_table)
    rho = LED.OptimiseLLH()

    return rho

rho = Simulate(iterations, bRate, dRate, sRate, mRate, popModel,
                                 susceptible, lockdownModel, randrange(sys.maxsize))
print(rho)


