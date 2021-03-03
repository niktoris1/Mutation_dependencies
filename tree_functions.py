import re

def getcite(mutation_name): # gets a cite name from the name of mutation
    return int(re.findall('\d+', mutation_name)[0])

class MutationOnNode: # defines a mutation on specific node
    def __init__(self, mutation_name, old_nucleotyde, new_nucleotyde, time_of_birth):
        self.mutation_name = mutation_name
        self.mutation_cite = getcite(mutation_name)
        self.old_nucleotyde = old_nucleotyde
        self.new_nucleotyde = new_nucleotyde
        self.time_of_birth = time_of_birth

name_id_dict = {} # dictionary in type id: name

def add_children(some_tree, node): #adds all children to the tree by checking all parent ids
    if node["parentid"] is not None:
        some_tree.create_node(node["name"], node["id"], parent=node["parentid"])
    else:
        some_tree.create_node(node["name"], node["id"])
    name_id_dict[node["name"]] = node["id"]
    for child_node in node["children"]:
        add_children(some_tree, child_node)