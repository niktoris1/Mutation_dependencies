from treelib import Node, Tree
import re

def getcite(mutation_name): # gets a cite name from the name of mutation
    return int(re.findall('\d+', mutation_name)[0])

class MutationOnNode: # defines a mutation on specific node
    def __init__(self, mutation_name, old_nucleotyde, new_nucleotyde, time_of_birth):
        self.mutation_name = mutation_name
        self.mutation_cite = getcite(mutation_name)
        self.old_nucleotyde = old_nucleotyde
        self.new_nucleotyde = new_nucleotyde
        self.time_of_birth = time_of_birth

def ArrayTreeToTreeClass(array_tree, array_times): # sets everythong with a placeholder mutations

    tree = Tree()

    tree.create_node(0, 0, data=MutationOnNode(mutation_name="0", old_nucleotyde="None", new_nucleotyde="None", time_of_birth=array_times[0]))

    for i in range(1, len(array_tree)):
        tree.create_node(i, i, parent=array_tree[i], data=MutationOnNode(mutation_name="0", old_nucleotyde="None", new_nucleotyde="None", time_of_birth=array_times[i]))

    return tree


class Event:
    def __init__(self, vertex_tag, event_type, event_time):
        self.event_type = event_type
        self.vertex_tag = vertex_tag
        self.event_time = event_time


# class EventSequence: #list of events, where the index is the order of iterations - might be deprecated
#    def __init__(self, event_sequence):
#        self.event_sequence = event_sequence

def GetTime(node):
    return node.data.time_of_birth


def GetEventsFromTree(tree):
    nodes_array = tree.all_nodes()
    events_array = [0] * len(nodes_array)

    for event_number in range(0, len(events_array)):
        events_array[event_number] = Event(vertex_tag=nodes_array[event_number].tag, \
                                           event_time=nodes_array[event_number].data.time_of_birth,
                                           event_type="Unknown")

    for event_number in range(0, len(events_array)):
        if len(tree.children(nodes_array[event_number].tag)) == 0:
            events_array[event_number].event_type = "adding_lineage"
        else:
            events_array[event_number].event_type = "coalescence"

    def takeBirth(elem):
        return elem.event_time

    events_array.sort(key=takeBirth)

    return events_array


def GetTimesFromEvents(events_array):
    times_array = [0] * len(events_array)

    for time_number in range(0, len(times_array)):
        times_array[time_number] = events_array[time_number].event_time

    return times_array