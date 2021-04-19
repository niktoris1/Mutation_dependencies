#!/usr/bin/env python3

import argparse
import sys
import time
from BirthDeathCython import BirthDeathModel, PopulationModel, Population, Lockdown
from IO import ReadRates, ReadPopulations, ReadMigrationRates, ReadSusceptibility
from random import randrange
import numpy as np
from get_subtree import SubtreeCreation, SubtreeCreation2

from tree_functions import ArrayTreeToTreeClass
from likelyhood_estimation import LikelyhoodEstimation

parser = argparse.ArgumentParser(description='Migration inference from PSMC.')

#parser.add_argument('frate',
#                    help='file with rates')

iterations = 50000
susceptibility = None
seed = None
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

# if clargs.lockdownModel == None:
#     lockdownModel = None
# else:
#     lockdownModel = ReadLockdown(clargs.lockdownModel)

if seed == None:
    rndseed = randrange(sys.maxsize)
else:
    rndseed = seed
print("Seed: ", rndseed)

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


tree = simulation.GetTree()
times = simulation.GetTimes()
mut = simulation.GetMut()
currentTime = simulation.GetCurrentTime()

newtree = ArrayTreeToTreeClass(tree, times, mut) # need to get self.tree, self.times and self.muts from BirthDeathClass - тут должно быть дерево, времена создания каждой из нод и мутации на нодах

#newtree.show()
#print("Time is", currentTime)

sc = SubtreeCreation2(A_nucleotyde = 'A', A_cite = 0, B_nucleotyde = 'A', B_cite = 1, tree = newtree)

subtree_AA = SubtreeCreation2.GetABsubtrees(sc)

#for tree in subtree_AA:
#    tree.show()

if len(subtree_AA) > 0:
    print('Subtree AA is not empty')
    ls_AA = LikelyhoodEstimation(subtree_AA)

    tree_size = 0
    for tree in subtree_AA:
        tree_size = tree_size + tree.size()

    print("Tree size is ", tree_size)

    t1 = time.time()
    es_ls_AA = ls_AA.GetEstimation()
    t2 = time.time()
    print('es_ls_AA =', es_ls_AA)
    print('Time spent on estimation: ', t2 - t1)

    time_start = 999

    for timestamp in [x / 10000 for x in range(1, 100000)]:
        if ls_AA.DistinctLineages(timestamp) == 5:
            print('On time', timestamp, 'there was', ls_AA.DistinctLineages(timestamp), 'distinct lineages')
            time_start = timestamp
            break

    print('Current time ', currentTime)
    print('Time start: ', time_start)
    if currentTime - time_start < 0:
        print('Never was 5 simulatanious linages')
    else:
        print('Time passed: ', currentTime - time_start)

else:
    print('Subtree AA is empty')




def writeMutations(mut):
    #digits replacement
    alleles = ["A","T","C","G"]
    for i in [1,3]:
        for j in range(len(mut[i])):
            mut[i][j] = alleles[mut[i][j]]

    mutations_dict = {}
    for nodeId in mut[0]:
        if nodeId in mutations_dict: #adding mutation for existing node
            mutations_dict[nodeId] += str(mut[1][mut[0].index(nodeId)]) \
                                      + str(mut[2][mut[0].index(nodeId)]) \
                                      + str(mut[3][mut[0].index(nodeId)])+','
        else:
            mutations_dict[nodeId] = str(mut[1][mut[0].index(nodeId)]) \
                                     + str(mut[2][mut[0].index(nodeId)]) \
                                     + str(mut[3][mut[0].index(nodeId)])+','
    #removing extra comma
    for nodeId in mutations_dict:
        mutations_dict[nodeId] = mutations_dict[nodeId][:-1]

    f_mut = open('mutation_output.tsv', 'w')
    for i in range(len(tree)):
        if i in mutations_dict:
            f_mut.write(str(i)+'\t'+str(mutations_dict[i])+'\n')
        else:
            f_mut.write(str(i)+'\n')
    f_mut.close()

# count of childrens
def frequentCart(nodes, sequence):
    result = dict()
    #print("nodes=", nodes)
    for node in nodes:
        result[node] = 0
    for parent in sequence:
        if parent in nodes:
            result[parent] = result[parent] + 1
    return result

# place of childrens
def allChildrens(nodes, sequence):
    result = dict()
    for node in nodes:
        result[node] = []
    for index in range(len(sequence)):
        if sequence[index] in nodes:
            result[sequence[index]].append(index)
    return result

def getOutputDict(nodes, times):
    result = dict()
    for node in nodes:
        result[node] = "{0}:{1}".format(node, times[node])
    return result

def phase3_LookForParents(resultOutput, listOfLeefs, pruferSeq, allChildren):
    alreadyFinishedParent = []
    #print('phase 3')
    parentFutureLeeves = []
    futureLeeves = []
    alreadyFinishedLeeves = []
    for leef in listOfLeefs:
        if leef in alreadyFinishedLeeves:
            continue
        parent = int(pruferSeq[leef])
        # root is found
        if parent == -1:
            resultOutput[leef] = "(" + resultOutput[leef] + ")"
            continue
        if parent in alreadyFinishedParent:
            continue
        alreadyFinishedParent.append(parent)
        listOfChildren = allChildren[parent]
        isAllChildrenLeefs = True
        for child in listOfChildren:
            if not child in listOfLeefs:
                isAllChildrenLeefs = False
                break
        if isAllChildrenLeefs:
            parentFutureLeeves.append(parent)
            message = ""
            for child in listOfChildren:
                childSplit = resultOutput[child].split(':')
                absTimeChild = float(childSplit[-1])
                absTimeParent = float(resultOutput[parent].split(':')[-1])
                time = absTimeChild - absTimeParent
                childSplit[-1] = str(time)
                resultOutput[child] = ":".join(childSplit)

                message += resultOutput[child] + "," ##?????
                alreadyFinishedLeeves.append(child)
                resultOutput.pop(child)
            resultOutput[parent] = "(" + message[:-1] + ")" + resultOutput[parent]#!!!!!
        else:
            futureLeeves.append(leef)
    return parentFutureLeeves, futureLeeves


def writeGenomeNewick(pruferSeq, times):
    #pruferSeq = pruferSeq.astype(int)
    for i in range(len(pruferSeq)):
        if pruferSeq[i] == i:
            pruferSeq[i] = -1
    #number of nodes
    numberOfNodes = len(pruferSeq)
    listOfNodes = [i for i in range(numberOfNodes)]
    frequencyCart = frequentCart(listOfNodes, pruferSeq)
    allChildren = allChildrens(listOfNodes, pruferSeq)
    resultOutput = getOutputDict(listOfNodes, times)

    #phase 2: look for normal leefs and parents
    listOfLeefs = []
    for key in frequencyCart:
        if frequencyCart[key] == 0:
            # find leefs
            resultOutput[key] = "{0}:{1}".format(key, times[key])
            listOfNodes.remove(key)
            listOfLeefs.append(key)

    #phase 3: look for parents
    parentFutureLeeves, futureLeeves = phase3_LookForParents(resultOutput, listOfLeefs, pruferSeq, allChildren)

    #phase 4: union lists
    while(True):
        listsOfNextLeeves = parentFutureLeeves + futureLeeves
        noParents = True
        actualLeafList = []
        for leef in listsOfNextLeeves:
            if pruferSeq[int(leef)] != -1.0:
                actualLeafList.append(leef)
                noParents = False
        if noParents:
             break
        else:
            parentFutureLeeves, futureLeeves = phase3_LookForParents(resultOutput, actualLeafList, pruferSeq, allChildren)

    f_nwk = open('newick_output.nwk', 'w')
    for key in resultOutput:
        f_nwk.write(resultOutput[key])
    f_nwk.write(';')
    f_nwk.close()
    #print(len(times))

if clargs.createNewick:
    writeGenomeNewick(tree, times)
if clargs.writeMutations:
    writeMutations(mut)
