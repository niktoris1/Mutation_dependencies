#!/usr/bin/env python3

import time
from BirthDeath.BirthDeath import BirthDeathModel, PopulationModel, Population
from BirthDeath.IO import ReadRates, ReadPopulations, ReadMigrationRates
from get_subtree import SubtreeCreation

from tree_functions import ArrayTreeToTreeClass
from likelyhood_estimation import LikelyhoodEstimation


frates_file = 'test/test.rt'

bRate, dRate, sRate, mRate = ReadRates(frates_file)
populationModel_args = None
debug_mode = False
iterations = 50000

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


newtree = ArrayTreeToTreeClass(simulation.genealogy, simulation.genealogyTimes, simulation.mutations_g)
#newtree.show()
#print("Time is", simulation.currentTime)
currentTime = simulation.currentTime

sc = SubtreeCreation()

subtree_AA = SubtreeCreation.GetABsubtrees(sc, A_nucleotyde = 'A', A_cite = 0, B_nucleotyde = 'A', B_cite = 1, base_tree = newtree)


#subtree_AT = SubtreeCreation.GetABsubtrees(sc, A_nucleotyde = 'A', A_cite = 0, B_nucleotyde = 'T', B_cite = 1, base_tree = newtree)
#print("GetATsubtrees")
#subtree_AC = SubtreeCreation.GetABsubtrees(sc, A_nucleotyde = 'A', A_cite = 0, B_nucleotyde = 'C', B_cite = 1, base_tree = newtree)
#print("GetACsubtrees")
#subtree_AG = SubtreeCreation.GetABsubtrees(sc, A_nucleotyde = 'A', A_cite = 0, B_nucleotyde = 'G', B_cite = 1, base_tree = newtree)
#print("GetAGsubtrees")
#subtree_TA = SubtreeCreation.GetABsubtrees(sc, A_nucleotyde = 'T', A_cite = 0, B_nucleotyde = 'A', B_cite = 1, base_tree = newtree)
#print("GetTAsubtrees")
#subtree_TT = SubtreeCreation.GetABsubtrees(sc, A_nucleotyde = 'T', A_cite = 0, B_nucleotyde = 'T', B_cite = 1, base_tree = newtree)
#print("GetTTsubtrees")
#subtree_TC = SubtreeCreation.GetABsubtrees(sc, A_nucleotyde = 'T', A_cite = 0, B_nucleotyde = 'C', B_cite = 1, base_tree = newtree)
#print("GetTCsubtrees")
#subtree_TG = SubtreeCreation.GetABsubtrees(sc, A_nucleotyde = 'T', A_cite = 0, B_nucleotyde = 'G', B_cite = 1, base_tree = newtree)
#print("GetTGsubtrees")
#subtree_CA = SubtreeCreation.GetABsubtrees(sc, A_nucleotyde = 'C', A_cite = 0, B_nucleotyde = 'A', B_cite = 1, base_tree = newtree)
#print("GetCAsubtrees")
#subtree_CT = SubtreeCreation.GetABsubtrees(sc, A_nucleotyde = 'C', A_cite = 0, B_nucleotyde = 'T', B_cite = 1, base_tree = newtree)
#print("GetCTsubtrees")
#print("GetCCsubtrees")
#subtree_CG = SubtreeCreation.GetABsubtrees(sc, A_nucleotyde = 'C', A_cite = 0, B_nucleotyde = 'G', B_cite = 1, base_tree = newtree)
#print("GetCGsubtrees")
#subtree_GA = SubtreeCreation.GetABsubtrees(sc, A_nucleotyde = 'G', A_cite = 0, B_nucleotyde = 'A', B_cite = 1, base_tree = newtree)
#print("GetGAsubtrees")
#subtree_GT = SubtreeCreation.GetABsubtrees(sc, A_nucleotyde = 'G', A_cite = 0, B_nucleotyde = 'T', B_cite = 1, base_tree = newtree)
#print("GetGTsubtrees")
#subtree_GC = SubtreeCreation.GetABsubtrees(sc, A_nucleotyde = 'G', A_cite = 0, B_nucleotyde = 'C', B_cite = 1, base_tree = newtree)
#print("GetGCsubtrees")
#subtree_GG = SubtreeCreation.GetABsubtrees(sc, A_nucleotyde = 'G', A_cite = 0, B_nucleotyde = 'G', B_cite = 1, base_tree = newtree)
#print("GetGGsubtrees")


if len(subtree_AA) > 0:
    ls_AA = LikelyhoodEstimation(subtree_AA)

    t1 = time.time()
    es_ls_AA = ls_AA.GetEstimation()
    t2 = time.time()
    print('es_ls_AA =', es_ls_AA)
    print('Time spent on estimation: ', t2 - t1)

    time_start = 999

    for timestamp in [x / 100 for x in range(1, 1000)]:
        if ls_AA.DistinctLineages(timestamp) == 5:
            print(timestamp, ls_AA.DistinctLineages(timestamp))
            time_start = timestamp
            break

    print('Current time ', currentTime)
    print('Time start: ', time_start)
    if currentTime - time_start < 0:
        print('Never was 5 similatanious linages')
    else:
        print('Time passed: ', currentTime - time_start)


# if len(subtree_AT) > 0:
#     ls_AT = LikelyhoodEstimation(subtree_AT)
#     t1 = time.time()
#     es_ls_AT = ls_AT.GetEstimation()
#     t2 = time.time()
#     print('es_ls_AT =', es_ls_AT)
#     print('Time spent on estimation: ', t2 - t1)

# if len(subtree_AC) > 0:
#     ls_AC = LikelyhoodEstimation(subtree_AC)
#     t1 = time.time()
#     es_ls_AC = ls_AC.GetEstimation()
#     t2 = time.time()
#     print('es_ls_AC =', es_ls_AC)
#     print('Time spent on estimation: ', t2 - t1)
#
# if len(subtree_AG) > 0:
#     ls_AG = LikelyhoodEstimation(subtree_AG)
#     t1 = time.time()
#     es_ls_AG = ls_AG.GetEstimation()
#     t2 = time.time()
#     print('es_ls_AG =', es_ls_AG)
#     print('Time spent on estimation: ', t2 - t1)
#
# if len(subtree_TA) > 0:
#     ls_TA = LikelyhoodEstimation(subtree_TA)
#     t1 = time.time()
#     es_ls_TA = ls_TA.GetEstimation()
#     t2 = time.time()
#     print('es_ls_TA =', es_ls_TA)
#     print('Time spent on estimation: ', t2 - t1)
#
# if len(subtree_TT) > 0:
#     ls_TT = LikelyhoodEstimation(subtree_TT)
#     t1 = time.time()
#     es_ls_TT = ls_TT.GetEstimation()
#     t2 = time.time()
#     print('es_ls_TT =', es_ls_TT)
#     print('Time spent on estimation: ', t2 - t1)
#
# if len(subtree_TC) > 0:
#     ls_TC = LikelyhoodEstimation(subtree_TC)
#     t1 = time.time()
#     es_ls_TC = ls_TC.GetEstimation()
#     t2 = time.time()
#     print('es_ls_TC =', es_ls_TC)
#     print('Time spent on estimation: ', t2 - t1)
#
# if len(subtree_TG) > 0:
#     ls_TG = LikelyhoodEstimation(subtree_TG)
#     t1 = time.time()
#     es_ls_TG = ls_TG.GetEstimation()
#     t2 = time.time()
#     print('es_ls_TG =', es_ls_TG)
#     print('Time spent on estimation: ', t2 - t1)
#
# if len(subtree_CA) > 0:
#     ls_CA = LikelyhoodEstimation(subtree_CA)
#     t1 = time.time()
#     es_ls_CA = ls_CA.GetEstimation()
#     t2 = time.time()
#     print('es_ls_CA =', es_ls_CA)
#     print('Time spent on estimation: ', t2 - t1)
#
# if len(subtree_CT) > 0:
#     ls_CT = LikelyhoodEstimation(subtree_CT)
#     t1 = time.time()
#     es_ls_CT = ls_CT.GetEstimation()
#     t2 = time.time()
#     print('es_ls_CT =', es_ls_CT)
#     print('Time spent on estimation: ', t2 - t1)
#
# if len(subtree_CC) > 0:
#     ls_CC = LikelyhoodEstimation(subtree_CC)
#     t1 = time.time()
#     es_ls_CC = ls_CC.GetEstimation()
#     t2 = time.time()
#     print('es_ls_CC =', es_ls_CC)
#     print('Time spent on estimation: ', t2 - t1)
#
# if len(subtree_CG) > 0:
#     ls_CG = LikelyhoodEstimation(subtree_CG)
#     t1 = time.time()
#     es_ls_CG = ls_CG.GetEstimation()
#     t2 = time.time()
#     print('es_ls_CG =', es_ls_CG)
#     print('Time spent on estimation: ', t2 - t1)
#
# if len(subtree_GA) > 0:
#     ls_GA = LikelyhoodEstimation(subtree_GA)
#     t1 = time.time()
#     es_ls_GA = ls_GA.GetEstimation()
#     t2 = time.time()
#     print('es_ls_GA =', es_ls_GA)
#     print('Time spent on estimation: ', t2 - t1)
#
# if len(subtree_GT) > 0:
#     ls_GT = LikelyhoodEstimation(subtree_GT)
#     t1 = time.time()
#     es_ls_GT = ls_GT.GetEstimation()
#     t2 = time.time()
#     print('es_ls_GT =', es_ls_GT)
#     print('Time spent on estimation: ', t2 - t1)
#
# if len(subtree_GC) > 0:
#     ls_GC = LikelyhoodEstimation(subtree_GC)
#     t1 = time.time()
#     es_ls_GC = ls_GC.GetEstimation()
#     t2 = time.time()
#     print('es_ls_GC =', es_ls_GC)
#     print('Time spent on estimation: ', t2 - t1)
#
# if len(subtree_GG) > 0:
#     ls_GG = LikelyhoodEstimation(subtree_GG)
#     t1 = time.time()
#     es_ls_GG = ls_GG.GetEstimation()
#     t2 = time.time()
#     print('es_ls_GG =', es_ls_GG)
#     print('Time spent on estimation: ', t2 - t1)

