import math

from get_subtree import test_tree


class Event:
    def __init__(self, vertex_tag, event_type, event_time):
        self.event_type = event_type
        self.vertex_tag = vertex_tag
        self.event_time = event_time



class EventSequence: #list of events, where the index is the order of iterations
    def __init__(self, event_sequence):
        self.event_sequence = event_sequence

def GetTime(node):
    return node.data.time_of_birth

def EventsFromTree(estimated_tree):
    events = estimated_tree.all_nodes()
    return events

def GetEventsFromTree(tree):


    events_array = EventsFromTree(tree)

    for event_number in range(0, len(events_array)):
        if events_array[event_number].data is None:
            events_array[event_number] = [events_array[event_number].tag, 0]
        else:
            events_array[event_number] = [events_array[event_number].tag, events_array[event_number].data.time_of_birth]

    def takeSecond(elem):
        return elem[1]

    events_array.sort(key=takeSecond)

    return events_array

events_array = GetEventsFromTree(test_tree)

for event_number in range(0, len(events_array)):
    events_array[event_number] = Event(vertex_tag = events_array[event_number][0], event_time = events_array[event_number][1], event_type='Unknown')
    if len(test_tree.children(events_array[event_number].vertex_tag)) == 0:
        events_array[event_number].event_type = "adding_lineage"
    else:
        events_array[event_number].event_type = "coalescence"


events_sequence = EventSequence(events_array)

