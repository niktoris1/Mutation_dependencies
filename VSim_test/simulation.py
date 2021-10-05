#!/usr/bin/env python3

import argparse
import sys
import time
from VGsim import BirthDeathModel, PopulationModel, Population, Lockdown
from VGsim.IO import ReadRates, ReadPopulations, ReadMigrationRates, ReadSusceptibility, ReadSusceptibilityTransition, writeGenomeNewick, writeMutations
from random import randrange



iterations = 1000000
susceptibility = None
seed = None
sampleSize = None
suscepTransition = None
populationModel = ['test/test.pp', 'test/test.mg']
frate = 'test/test.rt'

if isinstance(frate, list):
    frate = frate[0]
if isinstance(iterations, list):
    iterations = iterations[0]
if isinstance(sampleSize, list):
    sampleSize = sampleSize[0]
if isinstance(susceptibility, list):
    susceptibility = susceptibility[0]
if isinstance(suscepTransition, list):
    suscepTransition = suscepTransition[0]
if isinstance(seed, list):
    seed = seed[0]

bRate, dRate, sRate, mRate = ReadRates(frate)

if sampleSize == None:
    sampleSize = iterations

if populationModel == None:
    popModel = None
    lockdownModel = None
else:
    populations, lockdownModel, samplingMultiplier = ReadPopulations(populationModel[0])
    migrationRates = ReadMigrationRates(populationModel[1])
    popModel = [populations, migrationRates]

if susceptibility == None:
    susceptible = None
else:
    susceptible = ReadSusceptibility(susceptibility)

if suscepTransition == None:
    suscepTransition = None
else:
    suscepTransition = ReadSusceptibilityTransition(suscepTransition)

if seed == None:
    rndseed = randrange(sys.maxsize)
else:
    rndseed = seed
#print("Seed: ", rndseed)

#simulation = BirthDeathModel(bRate, dRate, sRate, mRate, populationModel=popModel, susceptible=susceptible,
#                             suscepTransition=suscepTransition, lockdownModel=lockdownModel, samplingMultiplier=samplingMultiplier, rndseed=rndseed)
# simulation.Debug()
# t1 = time.time()
#simulation.SimulatePopulation(iterations, sampleSize, time=1)
# simulation.Debug()
# t2 = time.time()
#simulation.GetGenealogy()
# simulation.Debug()
# t3 = time.time()
# simulation.Report()
# print(t2 - t1)
# print(t3 - t2)
#print("_________________________________")

#pruferSeq, times, mut, populations = simulation.Output_tree_mutations()