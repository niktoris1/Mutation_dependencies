#!/usr/bin/env python3

import argparse
import sys
import time
from Simulator.VGsim._BirthDeath import BirthDeathModel
from Simulator.VGsim.IO import ReadRates, ReadPopulations, ReadMigrationRates, ReadSusceptibility
from random import randrange
from get_subtree import SubtreeCreation

from tree_functions import EventsFromSimulation, TreeSequenceToTreeClass, TreeEventsFromSimulation, GetStartAndFinishtTimeFromTrees
from likelyhood_estimation import LikelyhoodEstimation

parser = argparse.ArgumentParser(description='Migration inference from PSMC.')

#parser.add_argument('frate',
#                    help='file with rates')

iterations = 100000
susceptibility = None
seed = 4733791303853627493 #5750693369156385614
populationModel = ['test/test.pp', 'test/test.mg']
frate = 'test/test.rt'


#parser.add_argument('--iterations', '-it', nargs=1, type=int, default=1000,
#                    help='number of iterations (default is 1000)')
#parser.add_argument('--populationModel', '-pm', nargs=2, default=None,
#                    help='population model: a file with population sizes etc, and a file with migration rate matrix')
#parser.add_argument('--susceptibility', '-su', nargs=1, default=None,
#                    help='susceptibility file')
# parser.add_argument('--lockdownModel', '-ld', nargs=1, default=None,
#                     help='lockdown model: a file with parameters for lockdowns')
#parser.add_argument('--seed', '-seed', nargs=1, type=int, default=None,
#                    help='random seed')
parser.add_argument("--createNewick", '-nwk',
                    help="Create a newick file of tree *.nwk ",
                    action="store_true")
parser.add_argument("--writeMutations", '-tsv',
                    help="Create a mutation file *.tsv ",
                    action="store_true")
parser.add_argument("--createTables", '-ctb',
                    help="Create event count and sample fraction tables",
                    action="store_true")

clargs = parser.parse_args()

if isinstance(frate, list):
    frate = frate[0]
if isinstance(iterations, list):
    iterations = iterations[0]
if isinstance(susceptibility, list):
    susceptibility = susceptibility[0]
# if isinstance(clargs.lockdownModel, list):
#     clargs.lockdownModel = clargs.lockdownModel[0]
if isinstance(seed, list):
    seed = seed[0]

bRate, dRate, sRate, mRate = ReadRates(frate)

if populationModel == None:
    popModel = None
    lockdownModel = None
else:
    populations, lockdownModel = ReadPopulations(populationModel[0])
    migrationRates = ReadMigrationRates(populationModel[1])
    popModel = [populations, migrationRates]

if susceptibility == None:
    susceptible = None
else:
    susceptible = ReadSusceptibility(susceptibility)


if seed == None:
    rndseed = randrange(sys.maxsize)
else:
    rndseed = seed
print("Seed: ", rndseed)

#lockdownModel = None # artificially switching off lockdown

#Beware,
simulation = BirthDeathModel(iterations, bRate, dRate, sRate, mRate, populationModel=popModel, susceptible=susceptible, lockdownModel=lockdownModel, rndseed=rndseed)
# simulation.Debug()
t1 = time.time()
simulation.SimulatePopulation(iterations)
#simulation.Debug()
t2 = time.time()
simulation.GetGenealogy()
# simulation.Debug()
t3 = time.time()
simulation.Report()
print(t2 - t1)
print(t3 - t2)
print("_________________________________")

treeevents = TreeEventsFromSimulation(simulation = simulation)
print('got events')
treeclasstree = TreeSequenceToTreeClass(simulation=simulation, tree_event_sequence=treeevents, is_AA_mutation_in_root_node = True)
print('converted tree')