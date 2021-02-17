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


class LikelyhoodEstimation:

    def __init__(self, estimated_tree, A_nucleotyde, B_nucleotyde):
        self.estimated_tree = estimated_tree
        self.A_nucleotyde = A_nucleotyde
        self.B_nucleotyde = B_nucleotyde

    def EventsFromTree(self, estimated_tree):
        events = estimated_tree.all_nodes()
        return events


    def LLH_function(self, iteration, coal_rate, coal_iteration, number_of_lineages, event_probability, event_type):
        return LikelyhoodEstimation.LLH_function(self, iteration - 1, coal_rate, coal_iteration, number_of_lineages, event_probability, event_type) * \
               math.exp(coal_rate * (LikelyhoodEstimation.TimeFromIteration(coal_iteration) - LikelyhoodEstimation.TimeFromIteration(coal_iteration - 1))) \
               * math.comb(LikelyhoodEstimation.DistinctLineages(LikelyhoodEstimation.TimeFromIteration(coal_iteration - 1)), 2) * \
               LikelyhoodEstimation.EventProbability(event_type, coal_rate, LikelyhoodEstimation.DistinctLineages(LikelyhoodEstimation.TimeFromIteration(iteration)))


    def TimeFromIteration(self, iteration):
        time = iteration ## to do
        return time

    def DistinctLineages(self, iteration):
        number_of_lineages = iteration
        return number_of_lineages

    def EventFromIteration(self, iteration):
        return 0



    def EventProbability(self, event_type, coal_rate, distinct_lineages):
        if event_type == 1:
            return 1
        else:
            probability = 1
            for i in range(0, distinct_lineages):
                probability = probability * coal_rate * math.comb(distinct_lineages - i + 1, 2)

def GetTime(node):
    return node.data.time_of_birth

def GetEventsFromTree(tree):

    le = LikelyhoodEstimation

    events_array = LikelyhoodEstimation.EventsFromTree(le, tree)

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


print(events_array)

