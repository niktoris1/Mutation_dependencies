import math

from get_subtree import test_tree


class Event:
    def __init__(self, vertex_tag, event_type, event_time):
        self.event_type = event_type
        self.vertex_tag = vertex_tag
        self.event_time = event_time



#class EventSequence: #list of events, where the index is the order of iterations - might be deprecated
#    def __init__(self, event_sequence):
#        self.event_sequence = event_sequence

def GetTime(node):
    return node.data.time_of_birth

def GetEventsFromTree(tree):

    nodes_array = tree.all_nodes()
    events_array = [0] * len(nodes_array)

    for event_number in range(0, len(events_array)):
        events_array[event_number] = Event(vertex_tag=nodes_array[event_number].tag, \
                                     event_time = nodes_array[event_number].data.time_of_birth, event_type="Unknown")

    for event_number in range(0, len(events_array)):
        if len(test_tree.children(nodes_array[event_number].tag)) == 0:
            events_array[event_number].event_type = "adding_lineage"
        else:
            events_array[event_number].event_type = "coalescence"

    def takeBirth(elem):
        return elem.event_time

    events_array.sort(key=takeBirth)

    return events_array

test_events_sequence = GetEventsFromTree(test_tree)

