import treelib
from treelib import Tree


def parse_parents_plus_times():
    with open('parents_times_events.txt', 'r') as file:
        data = file.read().split('\n')

    parents_str = data[1][1:-1].split(", ")
    parents = [int(item) for item in parents_str]

    times_str = data[4][1:-1].split(", ")
    times = [float(item) for item in times_str]

    return [parents, times]

def tree_from_array(parents_array, times_array):

    tree = Tree()

    tree.create_node(0, 0, data=times_array[0]) # id and tag are the same here

    for node_number in range(0, len(parents_array)):
        for potential_child_number in range(node_number + 1, len(parents_array)): # generating children
            if parents_array[potential_child_number] == node_number:
                tree.create_node(potential_child_number, potential_child_number, parent=node_number, \
                    data=times_array[potential_child_number] - times_array[node_number]) #data is a distance to the parent - from a new birth or from death

    return tree



[parents, times] = parse_parents_plus_times()
test_tree = tree_from_array(parents, times)
test_tree.show()

print(test_tree[3].data)

