import treelib
from treelib import Tree

def parse_parents_plus_times():
    with open('parents_plus_times.txt', 'r') as file:
        data = file.read().split('\n')

    parents_str = data[1][1:-1].split(", ")

    parents = [int(item) for item in parents_str]


    times_str = data[4][1:-1].split(", ")

    times = [float(item) for item in times_str]


    return [parents, times]

def tree_from_array(parents_array, times_array):

    tree = Tree()
    tree.create_node(0, 0, data=times_array[0]) # id and tag are the same here


    for i in range(len(parents_array) - 1):
        tree.create_node(i+1, i+1, parent=parents_array[i+1], data=times_array[i+1])

    return tree


[parents, times] = parse_parents_plus_times()
test_tree = tree_from_array(parents, times)
test_tree.show()