from treelib import Tree
import re

class Event:
    def __init__(self, type = None, time = None, iteration = None,
                 haplotype = None, node_id = None, number_of_children = None, current_sucseptible = None,
                 current_infectious = None, is_a_mutation = None, old_nucleotyde = None,
                 new_nucleotyde = None, mutation_cite = None, mutation_name = None):
        self.type = type
        self.time = time
        self.iteration = iteration
        self.haplotype = haplotype
        self.node_id = node_id
        self.number_of_children = number_of_children
        self.current_sucseptible = current_sucseptible
        self.current_infectious = current_infectious
        self.is_a_mutation = is_a_mutation
        if is_a_mutation == True:
            self.old_nucleotyde = old_nucleotyde
            self.new_nucleotyde = new_nucleotyde
            self.mutation_cite = mutation_cite
            self.mutation_name = mutation_name



class EventSequence:
    def __init__(self, sequence):
        self.sequence = sequence

    def TimeFromIteration(self, iteration):
        if iteration < 0:
            return 0
        return self.sequence[iteration].event_time

    def GetSlice(self, time_start, time_finish):
        start = None
        finish = None
        for i in range(len(self.sequence)):
            if self.sequence[i].event_sequence >= time_start:
                start = i

        for i in range(len(self.sequence), 0, -1):
            if self.sequence[i] <= time_finish:
                finish = i

        if start >= finish:
            raise ValueError("Slicing is incorrent")

        es = EventSequence(self.sequence[start : finish + 1])

        return es

    def GetAverageSucseptible(self):
        sucs = 0
        for i in range(len(self.s)):
            sucs = sucs + self.sequence[i].current_sucseptible

        return sucs / len(self.sequence)

    def GetCurrentInfectious(self):
        inf = 0
        for i in range(len(self.sequence)):
            inf = inf + self.sequence[i].current_infectious

        return inf / len(self.sequence)

def EventsFromSimulation(simulation, is_AA_mutation_in_root_node = False):
    es = EventSequence(sequence = [])

    def GetNodeIdByEventIteration(iteration):
        time = simulation.GetAllTimes()[len(simulation.GetAllTimes()) - 1 - iteration]
        tree_times= simulation.GetTreeTimes()
        if time in tree_times:
            node_id = tree_times.index(time)
            return node_id
        else:
            return None

    def IsThereAMutationOnNodeId(node_id):
        if node_id in simulation.GetTreeMuts[0]:
            return True
        else:
            return False

    def NumberOfChildrenFromNodeId(node_id):
        if node_id in simulation.GetTree():
            number_of_children = simulation.GetTree().count(node_id)
            return number_of_children
        else:
            return None

    def GetCite(mutation_name):  # gets a cite name from the name of mutation
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

    for i in range(simulation.GetNumberOfEvents()):
        type = simulation.GetEventTypes()[len(simulation.GetEventTypes()) - 1 - i] # here we have a 5 types of events
        time = simulation.GetAllTimes()[len(simulation.GetAllTimes()) - 1 - i] # might be incorrenct, since times are backwards
        iteration = i # the number of event in a sequence
        haplotype = simulation.GetHaplotypes()[len(simulation.GetHaplotypes()) - 1 - i] # which haplotype was in place, when event occured
        node_id = GetNodeIdByEventIteration(iteration)
        number_of_children = NumberOfChildrenFromNodeId(node_id)
        current_sucseptible = events.newSucseptibles[i] # it's worrying, that this array is forward-time, while others are backward-time
        current_infectious = events.newInfectious[i] # same goes here
        is_a_mutation = IsThereAMutationOnNodeId(node_id)
        old_nucleotyde = NumberToLetter(simulation.GetTreeMuts[1][node_id])
        new_nucleotyde = NumberToLetter(simulation.GetTreeMuts[3][node_id])
        mutation_cite = simulation.GetTreeMuts[2][node_id]
        mutation_name = GetCite(mutation_cite)

        event = Event(type = type, time = time, iteration = iteration,
                 haplotype = haplotype, node_id = node_id, number_of_children = number_of_children, current_sucseptible = current_sucseptible,
                 current_infectious = current_infectious, is_a_mutation = is_a_mutation, old_nucleotyde = old_nucleotyde,
                 new_nucleotyde = new_nucleotyde, mutation_cite = mutation_cite, mutation_name = mutation_name)

        es.sequence.append(event)

    if is_AA_mutation_in_root_node == True:
        es.sequence[0].is_a_mutation = True
        es.sequence[0].old_nucleotyde = 'A'
        es.sequence[0].new_nucleotyde = 'A'
        es.sequence[0].mutation_cite = 0

    return es

def EventSequenceToTreeClass(simulation, event_sequence):
    tree = Tree()

    root_id = 'Unknown'

    for i in range(len(event_sequence.sequence)):
        if simulation.GetTree[i] == -1:
            root_id = i
            tree.create_node(root_id, root_id, data=event_sequence[i])
            break  # there can be only one root

    if root_id == 'Unknown':
        raise ValueError("There is no root in this tree")

    for i in range(len(event_sequence.sequence)):
        if i != root_id:
            tree.create_node(i, i, parent=root_id, data=event_sequence[i])

    for i in range(len(event_sequence.sequence)):
        if i != root_id:
            tree.move_node(i, simulation.GetTree[i])

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

    es = EventSequence(sequence=[])

    for node in nodes_array:
        es.sequence.append(node.data)

    def takeBirth(elem):
        return elem.time

    es.sequence.sort(key=takeBirth)

    return es

