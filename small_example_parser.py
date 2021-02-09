import treelib
from treelib import Tree

def parse_parents_plus_times():
    with open('parents_times_events.txt', 'r') as file:
        data = file.read().split('\n')

    parents_str = data[1][1:-1].split(", ")
    parents = [int(item) for item in parents_str]


    times_str = data[4][1:-1].split(", ")
    times = [float(item) for item in times_str]

    events_str = data[7][1:-1].split(", ")
    events = [int(item) for item in events_str]


    return [parents, times, events]

def tree_from_array(parents_array, times_array, events_array):

    tree = Tree()
    tree.create_node(0, 0, data=times_array[0]) # id and tag are the same here


    for i in range(len(parents_array) - 1):
        tree.create_node(i+1, i+1, parent=parents_array[i+1], data=times_array[i+1] - times_array[parents_array[i+1]]) #data is a distance to the parent

    return tree

def event_sequence_from_tree():
    return 0


[parents, times, events] = parse_parents_plus_times()
test_tree = tree_from_array(parents, times)
#test_tree.show()