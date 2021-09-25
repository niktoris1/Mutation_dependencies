import random

from Simulator.VGsim._BirthDeath import BirthDeathModel
import matplotlib.pyplot as plt
from simulation import iterations, bRate, dRate, sRate, mRate, popModel, \
    susceptible, lockdownModel, suscepTransition, samplingMultiplier, sampleSize, rndseed
import statistics
from random import randrange
import numpy as np

from likelyhood_estimation_dismembered import LikelyhoodEstimationDismembered

import sys

def Simulate(iterations, bRate, dRate, sRate, mRate, popModel,
                                 susceptible, lockdownModel, suscepTransition, samplingMultiplier, rndseed):

    print("Random seed is: ", rndseed)

    #simulation = BirthDeathModel(iterations, bRate, dRate, sRate, mRate, populationModel=popModel,
    #                                 susceptible=susceptible, lockdownModel=lockdownModel,
    #                             suscepTransition=suscepTransition, samplingMultiplier=samplingMultiplier, rndseed=rndseed)
    #simulation.SimulatePopulation(iterations, sampleSize)
    #simulation.GetGenealogy()
    #tdm = simulation.gettdm() #get tdm object
    #trees_funct, trees_neutral = tdm.Dismember() #перед получением таблиц, нужно разчленить дерево
                #получение таблиц
    #event_table_funct, event_table_neutral = tdm.getEventTable() #[{time: [n_samples, n_coals]}]

    #nc, ns, fc, fs = 0, 0, 0, 0
    #for neutral_tree in event_table_neutral:
    #    for neutral_event in neutral_tree:
    #        if neutral_event[1] == 1:
    #            ns+=1
    #        if neutral_event[2] == 1:
    #            nc+=1
    #for funct_tree in event_table_funct:
    #    for funct_event in funct_tree:
    #        if funct_event[1] == 1:
    #            fs+=1
    #        if funct_event[2] == 1:
    #            fc+=1

    #print("Neutral coals:", nc)
    #print("Neutral samples:", ns)
    #print("Funct coals:", fc)
    #print("Funct samples:", fs)
    #print("Neutral trees", len(event_table_neutral))
    #print("Funct trees", len(event_table_funct))

    number_of_timestamps = 1 # how frequent brackets are
    #sample_fraction_table = tdm.getSampleFracTable(timestamps)

    def GenerateCoals(time, number_of_coals):
        return [[time, 0, 1] for _ in range(number_of_coals)]

    def GenerateSamples(time, number_of_coals):
        return [[time, 1, 0] for _ in range(number_of_coals)]

    def GenerateTree(step_time, number_of_levels):
        tree=[]
        for i in range(number_of_levels - 1):
            tree += GenerateCoals(step_time*i, pow(2, i))
        tree += GenerateSamples(step_time*number_of_levels, pow(2, number_of_levels-1))
        return tree

    def GenerateTreeFamily(step_time, number_of_levels):
        family = []

        for i in range(2, number_of_levels+1):
            family += GenerateTree(step_time, i)

        return family


    event_table_neutral_test=[
        GenerateTreeFamily(1, 9)
        #[[0, 0, 1],
        #[1, 1, 0], [1, 1, 0]]
        #[2, 1, 0], [2, 1, 0], [2, 1, 0], [2, 1, 0]]
        #[3, 1, 0],  [3, 1, 0], [3, 1, 0], [3, 1, 0], [3, 1, 0], [3, 1, 0], [3, 1, 0], [3, 1, 0]]
    ]


    event_table_funct_test=[
        GenerateTreeFamily(0.5, 14)
        #[[0, 0, 1],
        #[1, 0, 1], [1, 0, 1],
        #[2, 1, 0], [2, 1, 0], [2, 1, 0], [2, 1, 0]]
        #[6, 1, 0],  [6, 1, 0], [6, 1, 0], [6, 1, 0], [6, 1, 0], [6, 1, 0], [6, 1, 0], [6, 1, 0]]
    ]


    LED = LikelyhoodEstimationDismembered(event_table_funct=event_table_funct_test,
                                          event_table_neutral=event_table_neutral_test,
                                          number_of_brackets=1,
                                          simulation=None)
    #LED = LikelyhoodEstimationDismembered(et1, et2, 1, None)
    optimum = LED.OptimiseLLH()
    rho = optimum.x
    LLH_observed = optimum.fun
    #hd = simulation.GetHaplotypeDynamics(number_of_timestamps)
    #LLH_hypothesis = LED.GetLLHOptimumTotal(1)
    #LED.ConductLikelyhoodRatioTest(LLH_observed, LLH_hypothesis)
    #LED.PlotLLH()
    print("Rho equals:", rho)

    return rho, LLH_observed

#for randomiztion use randrange(sys.maxsize)
#793948375341945111 and 0.01 -
rho, LLH_observed = Simulate(iterations, bRate, dRate, sRate, mRate, popModel,
                                 susceptible, lockdownModel, suscepTransition, samplingMultiplier, 793948375341945111)




