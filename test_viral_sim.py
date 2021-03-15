#!/usr/bin/env python3

import argparse
import time
from BirthDeath import BirthDeathModel, PopulationModel, Population
from IO import ReadRates, ReadPopulations, ReadMigrationRates

from tree_functions import ArrayTreeToTreeClass
from likelyhood_estimation import LikelyhoodEstimation


frates_file = 'test/test.rt'

bRate, dRate, sRate, mRate = ReadRates(frates_file)
populationModel_args = None
debug_mode = False
iterations = 100

if populationModel_args == None:
    populationModel = PopulationModel([Population()], [[]])
else:
    populations = ReadPopulations(populationModel_args[0])
    migrationRates = ReadMigrationRates(populationModel_args[1])
    populationModel = PopulationModel(populations, migrationRates)

simulation = BirthDeathModel(bRate, dRate, sRate, mRate, debug = debug_mode, populationModel = populationModel)
t1 = time.time()
simulation.SimulatePopulation(iterations)
t2 = time.time()
simulation.GetGenealogy()
t3 = time.time()
simulation.Report()
print("Time to process the simulation - ", t2 - t1)
print("Time to process retrieve the genealogy - ", t3 - t2)


newtree = ArrayTreeToTreeClass(simulation.Tree, simulation.times)
newtree.show()
#print(simulation.times)

#newtree[0].data = MutationOnNode(mutation_name="1C", old_nucleotyde="G", new_nucleotyde="T", time_of_birth=0)


ls = LikelyhoodEstimation(newtree)

#print(ls.events_sequence)


estimation = ls.GetEstimation()
print(estimation)





# print(tree1.newTree)
# print(tree1.nodeSampling)
# print(tree1.times)
# print(tree1.newTimes)
