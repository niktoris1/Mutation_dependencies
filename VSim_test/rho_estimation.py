import random

from Simulator.VGsim._BirthDeath import BirthDeathModel
import matplotlib.pyplot as plt
from simulation import iterations, bRate, dRate, sRate, mRate, popModel, susceptible, lockdownModel, rndseed
import statistics
from random import randrange
import numpy as np

from likelyhood_estimation_dismembered import LikelyhoodEstimationDismembered

import sys



def Simulate(iterations, bRate, dRate, sRate, mRate, popModel,
                                 susceptible, lockdownModel, rndseed):

    print("Random seed is: ", rndseed)

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

    freq = 20 # how frequent brackets are
    timestamps = [_ for _ in np.linspace(0, max_time, freq + 1, endpoint=True)]
    #sample_fraction_table = tdm.getSampleFracTable(timestamps)
    LED = LikelyhoodEstimationDismembered(event_table_funct, event_table_neutral, timestamps)
    optimum = LED.OptimiseLLH()
    LED.PlotLLH()
    rho = optimum.x
    LLH_observed = optimum.fun
    lb = simulation.LogDynamics()
    print(lb)
    a=3
    #LLH_hypothesis = LED.GetLLHOptimumTotal(1)
    #LED.ConductLikelyhoodRatioTest(LLH_observed, LLH_hypothesis)
    #print("Rho equals:", rho)

    return rho, LLH_observed

#for randomiztion use randrange(sys.maxsize)

rho, LLH_observed = Simulate(iterations, bRate, dRate, sRate, mRate, popModel,
                                 susceptible, lockdownModel, 1834944703824815026)




