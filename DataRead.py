import re

class Node():
    def __init__(self, name):
        name = self.name
        parent = self.parent
        children = self.children

    def init_nodes(self):
        return 0

    def get_parent(self):
        return 0

    def get_children(self):
        return 0


def data_read():
    file = open("data.txt", "r")
    raw_data = [line.split(',') for line in file]

    name_array = []
    parent_array = []

    for element in raw_data[0]:
        nodes = re.split('\:', element)





        name_array.append(nodes)



    #for element in raw_data[0]:
    #    nodes = element.split(')')
    #    for node in nodes:
    #        name_array.append(node)

    #for name in name_array:
    #    name = name.replace('(', '')




    return name_array

print(data_read())