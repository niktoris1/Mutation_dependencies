import time
from likelyhood_estimation import LikelyhoodEstimation
from get_subtree import SubtreeCreation
from simulation import simulation


#newtree.show()
#print("Time is", currentTime)
A_nucleotyde = 'A'
B_nucleotyde = 'T'

sc = SubtreeCreation(A_nucleotyde = A_nucleotyde, A_cite = 0, B_nucleotyde = B_nucleotyde, B_cite = 1, tree = treeclasstree)

result = SubtreeCreation.GetABsubtrees(sc) # tree node data is empty


if len(result) > 0:
    print('Subtree ', A_nucleotyde, B_nucleotyde, ' is not empty')
    ls = LikelyhoodEstimation(result)

    tree_size = 0
    for tree in result:
        tree_size = tree_size + tree.size()

    print("Tree size is ", tree_size)

    t1 = time.time()
    es_ls = ls.GetEstimation()
    t2 = time.time()
    print('es_ls =', es_ls)
    print('Time spent on estimation: ', t2 - t1)

    time_start = 999

    for event in ls.es.tree_sequence:
        if ls.DistinctLineages(event.tree_time) == 10 and event.tree_time < time_start:
            time_start = event.tree_time

    print('Current time ', simulation.GetCurrentTime())
    print('Time start: ', time_start)
    time_passed = simulation.GetCurrentTime() - time_start
    if time_passed < 0:
        print('Never was 10 simulatanious linages')
    else:
        print('Time passed: ', time_passed)

    tree_size, coal_rate, program_time, time_passed = tree_size, es_ls[1], t2 - t1, time_passed

    #sucs = simulation.GetSucseptibles()
    #infs = simulation.GetInfectious()

    #plt.plot(sucs)
    #plt.plot(infs)

    #plt.show()
    total_sucs = 0
    total_infs = 0
    start_time, finish_time = GetStartAndFinishtTimeFromTrees(result)


    total_sucs = total_sucs + simulation.GetAverageSucseptiblesOnTimeframe(start_time, finish_time)
    total_infs = total_infs + simulation.GetAverageInfectiousOnTimeframe(start_time, finish_time)

    coal_rate_change = total_sucs / total_infs
    coal_rate = coal_rate * coal_rate_change

    print("Coal rate change is", coal_rate_change)
    print("Adjusted coal rate is", coal_rate)

else:
    print('Subtree is empty')

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
if clargs.createTables:
    print('tdm')
    tdm = simulation.gettdm() #get tdm object
    t5 = time.time()
    trees_funct, trees_neutral = tdm.Dismember() #перед получением таблиц, нужно разчленить дерево
    #получение таблиц
    event_table_funct, event_table_neutral = tdm.getEventTable() #[{time: [n_samples, n_coals]}]
    sample_fraction_table = tdm.getSampleFracTable([0.3, 0.6, 0.9, 1.2, 1.5, 1.8, 2.1, 2.4, 2.7]) # {time_bin: fraction}; fraction = -1 if I1 / 0
    t6 = time.time()
    #print(event_table_funct[0])
    tdm.debug()
    print('tdm done for', t6-t5)
