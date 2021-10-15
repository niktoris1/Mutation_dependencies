import os

from simulation import sampleSize

from Simulator.VGsim.interface import Simulator

simulator = Simulator(infectious_rate=30, uninfectious_rate=0, sampling_rate=15, sampling_probability=None,
                      sites_number=1, mutation_rate=1, mutation_probabilities=[0, 0, 0], populations_number=1,
                      population_size=100000, contact_density=1.0, total_migration_probability=0.0,
                      lockdown=[1, 100, 100],
                      sampling_multiplier=1.0, immunity_type=None, susceptibility=None, total_immunity_transition=None)

simulator.set_uninfectious_rate(1, 0)
simulator.set_uninfectious_rate(2, 0)
simulator.set_sampling_rate(1, 0)
simulator.set_sampling_rate(2, 0)

simulator.set_mutation_rate(0, 0, rate=1, probabilities=[0, 0, 0.1])
simulator.set_mutation_rate(1, 0, rate=0, probabilities=[0, 0, 1])
simulator.set_mutation_rate(2, 0, rate=0, probabilities=[0, 0, 1])
simulator.set_mutation_rate(3, 0, rate=1, probabilities=[0.1, 0, 0])

simulator.set_infectious_rate(3, 45)

simulator.print_basic_rates()

simulator.initialize(_seed=4810512211332756943)
simulator.simulate(_sampleSize=sampleSize, _iterations=200000)
simulator.epidemiology_timeline()
os.system('say "your program has finished"')

