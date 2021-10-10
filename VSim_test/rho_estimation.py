import random

from Simulator.VGsim._BirthDeath import BirthDeathModel
import matplotlib.pyplot as plt
from simulation import iterations, bRate, dRate, sRate, mRate, popModel, \
    susceptible, lockdownModel, suscepTransition, samplingMultiplier, sampleSize

from Simulator.VGsim.interface import Simulator

from likelyhood_estimation import LikelyhoodEstimationDismembered
from random import randrange

import sys

def Simulate(iterations, rndseed, frequency):

    print("Random seed is: ", rndseed)

    #simulator = Simulator(infection=30, uninfection=15, Sr=1, Sp=1, num_sites=1, mut_rate=1, mut_target_rate=[0, 0, 0],
    #                num_pop=1, size_pop=1000000, contact_density=1.0, total_mig_rate=0.0, lockdown=[1, 100, 100],
    #                sampling_multiplier=1.0, susc_type=None, susceptible=None, susc_trans=None)


    simulator = Simulator(infectious_rate=30, uninfectious_rate=0, sampling_rate=15, sampling_probability=None,
                 sites_number=1, mutation_rate=1, mutation_probabilities=[0, 0, 0], populations_number=1,
                 population_size=1000000, contact_density=1.0, total_migration_probability=0.0, lockdown=[1, 100, 100],
                 sampling_multiplier=1.0, immunity_type=None, susceptibility=None, total_immunity_transition=None)

    #simulator.set_infectious_rate(1, 0)
    #simulator.set_infectious_rate(2, 0)
    simulator.set_uninfectious_rate(1, 0)
    simulator.set_uninfectious_rate(2, 0)
    simulator.set_sampling_rate(1, 0)
    simulator.set_sampling_rate(2, 0)

    simulator.set_mutation_rate(0, 0, rate=1, probabilities=[0, 0, 1])
    simulator.set_mutation_rate(1, 0, rate=0, probabilities=[0, 0, 1])
    simulator.set_mutation_rate(2, 0, rate=0, probabilities=[0, 0, 1])
    simulator.set_mutation_rate(3, 0, rate=1, probabilities=[1, 0, 0])

    simulator.set_infectious_rate(3, 30)

    simulator.print_basic_rates()


    simulator.initialize(_seed=rndseed)
    simulator.simulate(_sampleSize=sampleSize, _iterations=iterations)
    simulator.epidemiology_timeline()
    tdm = simulator.simulation.gettdm() #get tdm object
    trees_funct, trees_neutral = tdm.Dismember() #перед получением таблиц, нужно разчленить дерево
                #получение таблиц
    event_table_funct, event_table_neutral = tdm.getEventTable() #[{time: [n_samples, n_coals]}]

    #log_dynamics = simulator.log_dynamics(step=100, output_file=False)

    LED = LikelyhoodEstimationDismembered(event_table_funct=event_table_funct,
                                          event_table_neutral=event_table_neutral,
                                          number_of_brackets=frequency,
                                          simulation=simulator)
    #LED = LikelyhoodEstimationDismembered(et1, et2, 1, None)
    optimum = LED.OptimiseLLH()
    rho = optimum.x
    LLH_observed = optimum.fun
    #LED.PlotLLH()
    print("Rho equals:", rho)

    return rho, LLH_observed

#for randomiztion use randrange(sys.maxsize)
#793948375341945111 and 0.01 -
rho, LLH_observed = Simulate(1000000, 5956145092605054515, frequency=20)




