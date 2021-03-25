from treelib import Node, Tree
import re

def getcite(mutation_name): # gets a cite name from the name of mutation
    return int(re.findall('\d+', mutation_name)[0])

def number_to_letter(number):
    if number == 0:
        return 'A'
    elif number == 1:
        return 'T'
    elif number == 2:
        return 'C'
    elif number == 3:
        return 'G'
    else:
        raise('Error, nucleotyde should be from 0 to 3')



class MutationOnNode: # defines a mutation on specific node
    def __init__(self, mutation_name, old_nucleotyde, new_nucleotyde, time_of_birth, mutation_cite = 9999999):
        self.mutation_name = mutation_name
        self.mutation_cite = mutation_cite if mutation_cite != 9999999 else getcite(mutation_name)
        self.old_nucleotyde = old_nucleotyde
        self.new_nucleotyde = new_nucleotyde
        self.time_of_birth = time_of_birth

def ArrayTreeToTreeClass(array_tree, array_times, array_mutations): # sets everythong with a placeholder mutations

    tree = Tree()

    for i in range(0, len(array_tree)):
        if array_tree[i] == -1:
            root_id = i
            tree.create_node(root_id, root_id, data=MutationOnNode(mutation_name="999999999", old_nucleotyde="None", new_nucleotyde="None", time_of_birth=array_times[root_id]))
            break # there can be only one root

    for i in range(0, len(array_tree)):
        if i != root_id:
            tree.create_node(i, i, parent=root_id, data=MutationOnNode(mutation_name="999999999", old_nucleotyde="None", new_nucleotyde="None", time_of_birth=array_times[i]))

    for i in range(0, len(array_tree)):
        if i != root_id:
            tree.move_node(i, array_tree[i])

    for mutation in array_mutations:
        tree.update_node(mutation.nodeId, data = MutationOnNode(mutation_name=str(number_to_letter(mutation.AS))+ \
            "to"+str(number_to_letter(mutation.DS))+"on"+str(mutation.time), old_nucleotyde=number_to_letter(mutation.AS), new_nucleotyde=number_to_letter(mutation.DS), \
                                                                   time_of_birth=array_times[mutation.nodeId], mutation_cite = mutation.position))

    return tree



class Event:
    def __init__(self, vertex_tag, event_type, event_time):
        self.event_type = event_type
        self.vertex_tag = vertex_tag
        self.event_time = event_time


def GetTime(node):
    return node.data.time_of_birth


def GetEventsFromTree(tree_list): # returns a list of events in chronological order
    nodes_array = []

    for tree in tree_list:
        for node in tree.all_nodes():
            nodes_array.append(node)

    events_array = [0] * len(nodes_array)

    def event_type_from_children(children):
        if children == 0:
            return "adding_lineage"
        else:
            return "coalescence"

    for event_number in range(0, len(events_array)):
        events_array[event_number] = Event(vertex_tag=nodes_array[event_number].tag, \
                                           event_time=nodes_array[event_number].data.time_of_birth,
                                           event_type=event_type_from_children(nodes_array[event_number].fpointer))

    def takeBirth(elem):
        return elem.event_time

    events_array.sort(key=takeBirth)

    return events_array


def GetTimesFromEvents(events_array):
    times_array = [0] * len(events_array)

    for time_number in range(0, len(times_array)):
        times_array[time_number] = events_array[time_number].event_time

    return times_array