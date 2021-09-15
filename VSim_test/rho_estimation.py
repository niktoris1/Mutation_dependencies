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

    nc, ns, fc, fs = 0, 0, 0, 0
    for neutral_tree in event_table_neutral:
        for neutral_event in neutral_tree:
            if neutral_event[2] == 1:
                nc+=1
            else:
                ns+=1
    for funct_tree in event_table_funct:
        for funct_event in funct_tree:
            if funct_event[2] == 1:
                fc+=1
            else:
                fs+=1

    print("Neutral coals:", nc)
    print("Neutral samples:", ns)
    print("Funct coals:", fc)
    print("Funct samples:", fs)



    number_of_timestamps = 50 # how frequent brackets are
    #sample_fraction_table = tdm.getSampleFracTable(timestamps)
    LED = LikelyhoodEstimationDismembered(event_table_funct, event_table_neutral, number_of_timestamps, simulation)
    optimum = LED.OptimiseLLH()
    LED.PlotLLH()
    rho = optimum.x
    LLH_observed = optimum.fun
    #hd = simulation.GetHaplotypeDynamics(number_of_timestamps)
    #LLH_hypothesis = LED.GetLLHOptimumTotal(1)
    #LED.ConductLikelyhoodRatioTest(LLH_observed, LLH_hypothesis)
    LED.PlotLLH()
    print("Rho equals:", rho)

    return rho, LLH_observed

#for randomiztion use randrange(sys.maxsize)

rho, LLH_observed = Simulate(iterations, bRate, dRate, sRate, mRate, popModel,
                                 susceptible, lockdownModel, 1834944703824815026)




