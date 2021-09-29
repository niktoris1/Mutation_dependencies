import random

from Simulator.VGsim._BirthDeath import BirthDeathModel
import matplotlib.pyplot as plt
from simulation import iterations, bRate, dRate, sRate, mRate, popModel, \
    susceptible, lockdownModel, suscepTransition, samplingMultiplier, sampleSize

from Simulator.VGsim.interface import Simulator

from likelyhood_estimation import LikelyhoodEstimationDismembered
from random import randrange

import sys

def Simulate(iterations, bRate, dRate, sRate, mRate, popModel,
                                 susceptible, lockdownModel, suscepTransition, samplingMultiplier, rndseed):

    print("Random seed is: ", rndseed)

    simulator = Simulator(infection=30, uninfection=15, Sr=1, Sp=1/3, num_sites=1, mut_rate=1, mut_target_rate=[0, 0, 0],
                    num_pop=1, size_pop=10000000, contact_density=1.0, total_mig_rate=0.0, lockdown=[1, 100, 100],
                    sampling_multiplier=1.0, susc_type=None, susceptible=None, susc_trans=None)

    simulator.set_Infection(1, 0)
    simulator.set_Infection(2, 0)
    simulator.set_Uninfection(1, 0)
    simulator.set_Uninfection(2, 0)
    simulator.set_S(1, 0)
    simulator.set_S(2, 0)

    simulator.set_MutRate(0, 0, 1, 1)
    simulator.set_MutRate(1, 0, 0, 0)
    simulator.set_MutRate(2, 0, 0, 0)
    simulator.set_MutRate(3, 0, 3, 1)

    simulator.print_Rates()
    simulator.create_class(rndseed)


    #simulation = BirthDeathModel(bRate, dRate, sRate, mRate, populationModel=popModel,
     #                                susceptible=susceptible, lockdownModel=lockdownModel,
     #                            suscepTransition=suscepTransition, samplingMultiplier=samplingMultiplier, rndseed=rndseed)
    simulator.simulation.SimulatePopulation(iterations, sampleSize, time=10)
    simulator.simulation.GetGenealogy()
    tdm = simulator.simulation.gettdm() #get tdm object
    trees_funct, trees_neutral = tdm.Dismember() #перед получением таблиц, нужно разчленить дерево
                #получение таблиц
    event_table_funct, event_table_neutral = tdm.getEventTable() #[{time: [n_samples, n_coals]}]



    LED = LikelyhoodEstimationDismembered(event_table_funct=event_table_funct,
                                          event_table_neutral=event_table_neutral,
                                          number_of_brackets=10,
                                          simulation=simulator.simulation)
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
                                 susceptible, lockdownModel, suscepTransition, samplingMultiplier, randrange(sys.maxsize))




