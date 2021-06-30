from treelib import Tree
import re
import time

class Event:
    def __init__(self, event_type = None, event_time = None, iteration = None,
                 haplotype = None, current_sucseptible = None,
                 current_infectious = None):
        self.event_type = event_type
        self.event_time = event_time
        self.iteration = iteration
        self.haplotype = haplotype
        self.current_sucseptible = current_sucseptible
        self.current_infectious = current_infectious


class TreeEvent:
    def __init__(self, node_id = None, tree_type = None, tree_time = None,
                 is_a_mutation = None, number_of_children = None, old_nucleotyde = None,
                 new_nucleotyde = None, mutation_cite = None):

        self.node_id = node_id
        self.tree_type = tree_type
        self.tree_time = tree_time
        self.number_of_children = number_of_children
        self.is_a_mutation = is_a_mutation
        if is_a_mutation == True:
            self.old_nucleotyde = old_nucleotyde
            self.new_nucleotyde = new_nucleotyde
            self.mutation_cite = mutation_cite
        else:
            self.old_nucleotyde = None
            self.new_nucleotyde = None
            self.mutation_cite = None


class EventSequence:
    def __init__(self, event_sequence):
        self.event_sequence = event_sequence

    def TimeFromIteration(self, iteration):
        if iteration < 0:
            return 0
        return self.event_sequence[iteration].event_time

    def GetSlice(self, time_start, time_finish):
        start = None
        finish = None
        for i in range(len(self.event_sequence)):
            if self.event_sequence[i].event_time >= time_start:
                start = i

        for i in range(len(self.event_sequence), 0, -1):
            if self.event_sequence[i].event_time <= time_finish:
                finish = i

        if start >= finish:
            raise ValueError("Slicing is incorrent")

        es = TreeEventSequence(self.event_sequence[start : finish + 1])

        return es

    def GetAverageSucseptible(self):
        sucs = 0
        for i in range(len(self.event_sequence)):
            sucs = sucs + self.event_sequence[i].current_sucseptible

        return sucs / len(self.event_sequence)

    def GetCurrentInfectious(self):
        inf = 0
        for i in range(len(self.event_sequence)):
            inf = inf + self.event_sequence[i].current_infectious

        return inf / len(self.event_sequence)


class TreeEventSequence:
    def __init__(self, tree_sequence):
        self.tree_sequence = tree_sequence

    def TimeFromIteration(self, iteration):
        if iteration < 0:
            return 0
        return self.tree_sequence[iteration].tree_time

    def GetSlice(self, time_start, time_finish):
        start = None
        finish = None
        for i in range(len(self.tree_sequence)):
            if self.tree_sequence[i].tree_time >= time_start:
                start = i

        for i in range(len(self.tree_sequence), 0, -1):
            if self.tree_sequence[i].tree_time <= time_finish:
                finish = i

        if start >= finish:
            raise ValueError("Slicing is incorrent")

        ts = TreeEventSequence(self.tree_sequence[start : finish + 1])

        return ts



def EventsFromSimulation(simulation, is_AA_mutation_in_root_node = False):
    es = EventSequence(event_sequence = [])

    for i in range(simulation.GetNumberOfEvents()):
        event_type = simulation.GetEventTypes()[len(simulation.GetEventTypes()) - 1 - i] # here we have a 5 types of events
        event_time = simulation.GetAllTimes()[len(simulation.GetAllTimes()) - 1 - i] # might be incorrenct, since times are backwards
        iteration = i # the number of event in a sequence
        haplotype = simulation.GetHaplotypes()[len(simulation.GetHaplotypes()) - 1 - i] # which haplotype was in place, when event occured
        current_sucseptible = simulation.GetSucseptibles()[i] # it's worrying, that this array is forward-time, while others are backward-time
        current_infectious = simulation.GetInfectious()[i] # same goes here
        event = Event(event_type= event_type, event_time= event_time, iteration = iteration,
                      haplotype = haplotype, current_sucseptible = current_sucseptible,
                      current_infectious = current_infectious)
        es.event_sequence.append(event)
    return es


def TreeEventsFromSimulation(simulation, is_AA_mutation_in_root_node = False):

    tes = TreeEventSequence(tree_sequence = [])

    def GetNodeIdByEventIteration(iteration):
        t1 = time.time()
        curtime = simulation.GetAllTimes()[iteration]
        if curtime in simulation.GetTreeTimes():
            t2 = time.time()
            i = simulation.GetTreeTimes().index()
            return i
        else:
            t2 = time.time()
            return None

    def IsThereAMutationOnNodeId(node_id):
        if node_id in simulation.GetTreeMutsNodeIds():
            return True
        else:
            return False

    def NumberOfChildrenFromNodeId(node_id):
        if node_id in simulation.GetTree():
            number_of_children = simulation.GetTree().count(node_id)
            return number_of_children
        else:
            return None

    def GetSite(mutation_name):  # gets a cite name from the name of mutation
        return int(re.findall('\d+', mutation_name)[0])

    def NumberToLetter(number):
        if number == 0:
            return 'A'
        if number == 1:
            return 'T'
        if number == 2:
            return 'C'
        if number == 3:
            return 'G'

    a = 0
    print(simulation.GetI())

    len_tree = len(simulation.GetTree())
    number_of_children_array = [0 for _ in range(len_tree)]
    tree_muts_ASs = simulation.GetTreeMutsASs()
    tree_muts_DSs = simulation.GetTreeMutsDSs()
    tree_muts_sites = simulation.GetTreeMutsSites()




    for i in range(len_tree):
        whose_child = simulation.GetTree()[i]
        number_of_children_array[whose_child] = number_of_children_array[whose_child] + 1

    for i in range(len_tree):
        node_id = i
        tree_time = simulation.GetTreeTimes()[i]
        number_of_children = number_of_children_array[i]
        if number_of_children == 0:
            tree_type = 'adding lineage'
        else:
            tree_type = 'coalescence'

        if node_id in simulation.GetTreeMutsNodeIds():
            is_a_mutation = True
            old_nucleotyde = tree_muts_ASs[node_id]
            new_nucleotyde = tree_muts_DSs[node_id]
            mutation_cite = tree_muts_sites[node_id]
        else:
            is_a_mutation = False
            old_nucleotyde = None
            new_nucleotyde = None
            mutation_cite = None

        tree_event = TreeEvent(node_id = node_id, tree_type = tree_type, tree_time = tree_time,
                 is_a_mutation = is_a_mutation, number_of_children = number_of_children, old_nucleotyde = old_nucleotyde,
                 new_nucleotyde = new_nucleotyde, mutation_cite = mutation_cite) # we do not initialize children here
        a = a+1
        print(a)
        tes.tree_sequence.append(tree_event)

    return tes


def TreeSequenceToTreeClass(simulation, tree_event_sequence):

    tree_size = len(tree_event_sequence.tree_sequence)

    tree = Tree()
    root_id = 'Unknown'
    for i in range(tree_size):
        if simulation.GetTree()[i] == -1:
            root_id = i
            tree.create_node(root_id, root_id, data=None) # placeholder on root
            break  # there can be only one root

    if root_id == 'Unknown':
        raise ValueError("There is no root in this tree")

    for i in range(tree_size):
        if i != root_id:
            tree.create_node(i, i, parent=root_id, data=None) # placeholder on other places

    for i in range(tree_size):
        if i != root_id:
            tree.move_node(i, simulation.GetTree()[i])

    for i in range(tree_size):
        noc = len(tree.get_node(i).fpointer) # number of children
        if i == root_id:
            tree.update_node(i, data=TreeEvent(is_a_mutation = True, number_of_children = noc, old_nucleotyde = 'A',
                                                  new_nucleotyde = 'A', mutation_cite = 0))
        else:
            tree.update_node(i, data=TreeEvent(is_a_mutation = None, number_of_children = noc, old_nucleotyde = None,
                                                  new_nucleotyde = None, mutation_cite = None))
    return tree


def IterationFromTime(time, es):
    def IterationFromTimeStartFinish(time, start, finish):
        middle = (start + finish) // 2

        if time <= es.sequence[start].event_time:
            return start
        if time >= es.sequence[finish].event_time:
            return finish

        if time == es.sequence[middle].event_time:
            return middle
        if time > es.sequence[middle].event_time:
            return IterationFromTimeStartFinish(time, middle + 1, finish)
        else:
            return IterationFromTimeStartFinish(time, start, middle - 1)

    return IterationFromTimeStartFinish(time, 0, len(es.sequence) - 1)


def GetEventsFromTree(tree_list):

    nodes_array = []

    for tree in tree_list:
        nodes_array = nodes_array + tree.all_nodes()

    es = TreeEventSequence(sequence=[])

    for node in nodes_array:
        es.sequence.append(node.data)

    def takeBirth(elem):
        return elem.time

    es.sequence.sort(key=takeBirth)

    return es

