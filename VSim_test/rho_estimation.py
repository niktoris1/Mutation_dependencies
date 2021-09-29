import random

from Simulator.VGsim._BirthDeath import BirthDeathModel
import matplotlib.pyplot as plt
from simulation import iterations, bRate, dRate, sRate, mRate, popModel, \
    susceptible, lockdownModel, suscepTransition, samplingMultiplier, sampleSize

from likelyhood_estimation import LikelyhoodEstimationDismembered

import sys

def Simulate(iterations, bRate, dRate, sRate, mRate, popModel,
                                 susceptible, lockdownModel, suscepTransition, samplingMultiplier, rndseed):

    print("Random seed is: ", rndseed)

    simulation = BirthDeathModel(bRate, dRate, sRate, mRate, populationModel=popModel,
                                     susceptible=susceptible, lockdownModel=lockdownModel,
                                 suscepTransition=suscepTransition, samplingMultiplier=samplingMultiplier, rndseed=rndseed)
    simulation.SimulatePopulation(iterations, sampleSize, time=1)
    simulation.GetGenealogy()
    tdm = simulation.gettdm() #get tdm object
    trees_funct, trees_neutral = tdm.Dismember() #перед получением таблиц, нужно разчленить дерево
                #получение таблиц
    event_table_funct, event_table_neutral = tdm.getEventTable() #[{time: [n_samples, n_coals]}]


    number_of_timestamps = 1 # how frequent brackets are


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


    #event_table_neutral_test=[
    #    GenerateTreeFamily(1, 12)
    #]


    #event_table_funct_test=[
    #    GenerateTreeFamily(1, 18)
    #]



    #frate1 = 'test/test_tree_1.rt'
    #bRate, dRate, sRate, mRate = ReadRates(frate1)
    #simulation_1 = BirthDeathModel(iterations, bRate, dRate, sRate, mRate, populationModel=popModel,
    #                                 susceptible=susceptible, lockdownModel=lockdownModel,
    #                             suscepTransition=suscepTransition, samplingMultiplier=samplingMultiplier, rndseed=randrange(sys.maxsize))
    #simulation_1.SimulatePopulation(iterations, sampleSize)
    #simulation_1.GetGenealogy()
    #tdm = simulation_1.gettdm() #get tdm object
    #trees_funct_1, trees_neutral_1 = tdm.Dismember() #перед получением таблиц, нужно разчленить дерево
    # получение таблиц
    #event_table_funct_1, event_table_neutral_1 = tdm.getEventTable() #[{time: [n_samples, n_coals]}]



    #frate2 = 'test/test_tree_2.rt'
    #bRate, dRate, sRate, mRate = ReadRates(frate2)
    #simulation_2 = BirthDeathModel(iterations, bRate, dRate, sRate, mRate, populationModel=popModel,
    #                                 susceptible=susceptible, lockdownModel=lockdownModel,
    #                             suscepTransition=suscepTransition, samplingMultiplier=samplingMultiplier, rndseed=randrange(sys.maxsize))
    #simulation_2.SimulatePopulation(iterations, sampleSize)
    #simulation_2.GetGenealogy()
    #tdm = simulation_2.gettdm() #get tdm object
    #trees_funct_2, trees_neutral_2 = tdm.Dismember() #перед получением таблиц, нужно разчленить дерево
    # получение таблиц
    #event_table_funct_2, event_table_neutral_2 = tdm.getEventTable() #[{time: [n_samples, n_coals]}]

    #for event_tree_num in range(len(event_table_neutral_2)): # change of time
    #    for event_num in range(len(event_table_neutral_2[event_tree_num])):
    #        time =  event_table_neutral_2[event_tree_num][event_num][0]
    #        sample = event_table_neutral_2[event_tree_num][event_num][1]
    #        coal = event_table_neutral_2[event_tree_num][event_num][2]
    #        event_table_neutral_2[event_tree_num][event_num] = [time/2+(random.random()-0.5)/1000, sample, coal]




    LED = LikelyhoodEstimationDismembered(event_table_funct=event_table_funct,
                                          event_table_neutral=event_table_neutral,
                                          number_of_brackets=10,
                                          simulation=simulation)
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




