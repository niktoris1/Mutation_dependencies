import math





class Event:
    def __init__(self, event_type, number_of_lineages, lineages):
        self.event_type = event_type
        self.number_of_lineages = number_of_lineages
        self.lineages = lineages

        if event_type is not 'adding_lineage' or 'coalescence':
            raise('Error - wrong event type')
        if event_type is 'adding_lineage' and len(lineages) is not 1:
            raise('Error - trying to add more or less than 1 lineage')

    def LineagesChange(self):
        if self.event_type == 'adding_lineage':
            return 1
        else:
            return - (self.number_of_lineages - 1)



class EventSequence: #list of events, where the index is the order of iterations
    def __init__(self, event_sequence):
        self.event_sequence = event_sequence

    def DistinctLineages(self, iteration):
        number_of_lineages = EventSequence.DistinctLineages(self, iteration - 1) + self.event_sequence[iteration].LineagesChange()
        return number_of_lineages


class LikelyhoodEstimation:

    def __init__(self, estimated_tree, A_nucleotyde, B_nucleotyde):
        self.estimated_tree = estimated_tree
        self.A_nucleotyde = A_nucleotyde
        self.B_nucleotyde = B_nucleotyde



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



    def EventProbability(event_type, coal_rate, distinct_lineages):
        if event_type == 1:
            return 1
        else:
            probability = 1
            for i in range(0, distinct_lineages):
                probability = probability * coal_rate * math.comb(distinct_lineages - i + 1, 2)