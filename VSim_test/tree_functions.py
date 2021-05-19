from treelib import Tree
import re

def getcite(mutation_name): # gets a cite name from the name of mutation
    return int(re.findall('\d+', mutation_name)[0])

def number_to_letter(number):
    if number == 0:
        return 'A'
    if number == 1:
        return 'T'
    if number == 2:
        return 'C'
    if number == 3:
        return 'G'


class MutationOnNode: # defines a mutation on specific node
    def __init__(self, mutation_name, old_nucleotyde, new_nucleotyde, time_of_birth, mutation_cite = 9999999):
        self.mutation_name = mutation_name
        self.mutation_cite = mutation_cite if mutation_cite != 9999999 else getcite(mutation_name)
        self.old_nucleotyde = old_nucleotyde
        self.new_nucleotyde = new_nucleotyde
        self.time_of_birth = time_of_birth

def ArrayTreeToTreeClass(array_tree, array_times, array_mutations, is_AA_mutation_in_root_node = False):

    # sets everythong with a placeholder mutations

    # if is_AA_mutation_in_root_node == True, we assume that we have a predefined AA nucleotydes in the root node

    tree = Tree()

    root_id = 'Unknown'
    for i in range(len(array_tree)):
        if array_tree[i] == -1:
            root_id = i
            tree.create_node(root_id, root_id, data=MutationOnNode(mutation_name="999999999", old_nucleotyde="None", new_nucleotyde="None", time_of_birth=array_times[root_id]))
            break # there can be only one root

    if root_id == 'Unknown':
        raise ValueError("There is no root in this tree")

    for i in range(len(array_tree)):
        if i != root_id:
            tree.create_node(i, i, parent=root_id, data=MutationOnNode(mutation_name="999999999", old_nucleotyde="None", new_nucleotyde="None", time_of_birth=array_times[i]))

    for i in range(len(array_tree)):
        if i != root_id:
            tree.move_node(i, array_tree[i])

    if tree.root not in array_mutations[0] and is_AA_mutation_in_root_node == True: # adding a mutation on tree root
        array_mutations[0].append(tree.root)
        array_mutations[1].append(0) # old nucleotyde is zero
        array_mutations[2].append(0) # we consider a cite, there we had a mutation as zero
        array_mutations[3].append(0) # new nucleotyde is zero
        array_times.append(0)

    for i in range(len(array_mutations[0])):
        tree.update_node(array_mutations[0][i], data = MutationOnNode(mutation_name=str(number_to_letter(array_mutations[1][i]))+
            "_to_"+str(number_to_letter(array_mutations[3][i]))+"_on_time_"+str(array_times[i]), old_nucleotyde=number_to_letter(array_mutations[1][i]), new_nucleotyde=number_to_letter(array_mutations[3][i]),
                                                                   time_of_birth=array_times[array_mutations[0][i]], mutation_cite = array_mutations[2][i]))

    return tree


class Event:
    def __init__(self, vertex_tag, event_type, event_time, number_of_children, vertex_id):
        self.event_type = event_type
        self.vertex_tag = vertex_tag
        self.event_time = event_time
        self.number_of_children = number_of_children
        self.vertex_id = vertex_id

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


    events_array = []

    for event_number in range(len(nodes_array)):
        events_array.append(Event(vertex_tag=nodes_array[event_number].tag,
                                           event_time=nodes_array[event_number].data.time_of_birth,
                                           event_type=EventTypeFromNode(nodes_array[event_number]),
                                           number_of_children=len(nodes_array[event_number].fpointer),
                                           vertex_id = nodes_array[event_number].identifier))

    def takeBirth(elem):
        return elem.event_time

    events_array.sort(key=takeBirth)

    return events_array


def GetTimesFromEvents(events_array):
    times_array = [0] * len(events_array)

    for time_number in range(len(times_array)):
        times_array[time_number] = events_array[time_number].event_time

    return times_array