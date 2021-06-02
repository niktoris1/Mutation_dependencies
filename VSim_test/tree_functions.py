from treelib import Tree
import re

def number_to_letter(number):
    if number == 0:
        return 'A'
    if number == 1:
        return 'T'
    if number == 2:
        return 'C'
    if number == 3:
        return 'G'

class Event:
    def __init__(self, vertex_tag, event_type, event_time,
                 number_of_children, vertex_id, iteration = None, haplotype=None):
        self.event_type = event_type
        self.vertex_tag = vertex_tag
        self.event_time = event_time
        self.number_of_children = number_of_children
        self.vertex_id = vertex_id
        self.iteration = iteration
        self.haplotype = haplotype

class EventSequence:
    def __init__(self, events_sequence):
        self.events_sequence = events_sequence

    def TimeFromIteration(self, iteration):
        if iteration < 0:
            return 0
        return self.events_sequence[iteration].event_time

class DataOnNode: # defines a mutation on specific node
    def __init__(self, is_mutation, mutation_name = None, old_nucleotyde = None,
                 new_nucleotyde = None, time_of_birth = None, mutation_cite = None,
                 haplotype=None, current_sucseptible = None, current_infectious = None):
        self.is_mutation = is_mutation

        if is_mutation == True:
            self.mutation_name = mutation_name
            self.mutation_cite = mutation_cite
            self.old_nucleotyde = old_nucleotyde
            self.new_nucleotyde = new_nucleotyde
        else:
            self.mutation_name = None
            self.mutation_cite = None
            self.old_nucleotyde = None
            self.new_nucleotyde = None

        self.time_of_birth = time_of_birth
        self.haplotype = haplotype
        self.current_sucseptible = current_sucseptible
        self.current_infectious = current_infectious


    def getcite(self, mutation_name):  # gets a cite name from the name of mutation
        return int(re.findall('\d+', mutation_name)[0])

def ArraysToTreeClass(array_tree, array_times, array_mutations, sucseptible_array, infectious_array, is_AA_mutation_in_root_node = False):

        # sets everythong with a placeholder mutations

        # if is_AA_mutation_in_root_node == True, we assume that we have a predefined AA nucleotydes in the root node

    tree = Tree()

    root_id = 'Unknown'
    for i in range(len(array_tree)):
        if array_tree[i] == -1:
            root_id = i
            tree.create_node(root_id, root_id, data=DataOnNode(is_mutation=False, time_of_birth=array_times[root_id],
                                                               current_sucseptible = sucseptible_array[i], current_infectious=infectious_array[i]))
            break # there can be only one root

    if root_id == 'Unknown':
        raise ValueError("There is no root in this tree")

    for i in range(len(array_tree)):
        if i != root_id:
            tree.create_node(i, i, parent=root_id, data=DataOnNode(is_mutation=False, time_of_birth=array_times[i],
                                                                   current_sucseptible=sucseptible_array[i],
                                                                   current_infectious=infectious_array[i]
                                                                   ))

    for i in range(len(array_tree)):
        if i != root_id:
            tree.move_node(i, array_tree[i])

    if tree.root not in array_mutations[0] and is_AA_mutation_in_root_node == True: # adding a mutation on tree root
        array_mutations[0].append(tree.root)
        array_mutations[1].append(None) # old nucleotyde is zero
        array_mutations[2].append(None) # we consider a cite, there we had a mutation as zero
        array_mutations[3].append(None) # new nucleotyde is zero
        array_times.append(0)

    for i in range(len(array_mutations[0])): # i - это номер, представленный в списке мутаций. Совпадает с номером в списке эвентов
        tree.update_node(array_mutations[0][i], data = DataOnNode(is_mutation=True, mutation_name=str(number_to_letter(array_mutations[1][i])) +
            "_to_" + str(number_to_letter(array_mutations[3][i])) +"_on_time_" + str(array_times[i]), old_nucleotyde=number_to_letter(array_mutations[1][i]),
                                                                  new_nucleotyde=number_to_letter(array_mutations[3][i]),
                                                                  time_of_birth=array_times[array_mutations[0][i]], mutation_cite = array_mutations[2][i],
                                                                  current_sucseptible = sucseptible_array[array_mutations[0][i]], current_infectious = infectious_array[array_mutations[0][i]])
                         )

    return tree


def IterationFromTime(time, es):
    def IterationFromTimeStartFinish(time, start, finish):
        middle = (start + finish) // 2

        if time <= es.events_sequence[start].event_time:
            return start
        if time >= es.events_sequence[finish].event_time:
            return finish

        if time == es.events_sequence[middle].event_time:
            return middle
        if time > es.events_sequence[middle].event_time:
            return IterationFromTimeStartFinish(time, middle + 1, finish)
        else:
            return IterationFromTimeStartFinish(time, start, middle - 1)

    return IterationFromTimeStartFinish(time, 0, len(es.events_sequence) - 1)



def GetTime(node):
    return node.data.time_of_birth

def GetEventsFromTree(tree_list):

    nodes_array = []

    for tree in tree_list:
        nodes_array = nodes_array + tree.all_nodes()

    def EventTypeFromNode(node):
        if len(node.fpointer) == 0:
            return "adding_lineage"
        else:
            return "coalescence"


    es = EventSequence(events_sequence=[])

    for event_number in range(len(nodes_array)):
        es.events_sequence.append(Event(vertex_tag=nodes_array[event_number].tag,
                                           event_time=nodes_array[event_number].data.time_of_birth,
                                           event_type=EventTypeFromNode(nodes_array[event_number]),
                                           number_of_children=len(nodes_array[event_number].fpointer),
                                           vertex_id = nodes_array[event_number].identifier,
                                           haplotype=None))

    def takeBirth(elem):
        return elem.event_time

    es.events_sequence.sort(key=takeBirth)

    return es


def GetTimesFromEvents(events_array):
    times_array = [0] * len(events_array)

    for time_number in range(len(times_array)):
        times_array[time_number] = events_array[time_number].event_time

    return times_array

def GetTotalSucseptibleByTree(tree):
    ts = 0
    for node in tree.all_nodes():
        ts = ts + node.data.current_sucseptible
    return ts

def GetTotalInfectiousByTree(tree):
    ti = 0
    for node in tree.all_nodes():
        ti = ti + node.data.current_infectious
    return ti
