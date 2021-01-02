import re
from treelib import Node, Tree


with open('data.txt', 'r') as file:
    data = file.read().replace('\n', '')

def parse(newick):
    tokens = re.findall(r"([^:;,()\s]*)(?:\s*:\s*([\d.]+)\s*)?([,);])|(\S)", newick+";")

    def recurse(nextid = 0, parentid = -1): # one node
        thisid = nextid;
        children = []

        name, length, delim, ch = tokens.pop(0)
        if ch == "(":
            while ch in "(,":
                node, ch, nextid = recurse(nextid+1, thisid)
                children.append(node)
            name, length, delim, ch = tokens.pop(0)
        return {"id": thisid, "name": name, "length": float(length) if length else None,
                "parentid": parentid, "children": children}, delim, nextid

    return recurse()[0]

raw_nodes = parse(data)


tree = Tree()
tree.create_node(raw_nodes["name"], raw_nodes["id"])

def add_children(node):
    tree.create_node(node["name"], node["id"], parent=node["parentid"])
    for child_node in node["children"]:
        add_children(child_node)

for children_node in raw_nodes["children"]:
    add_children(children_node)

tree.save2file('tree.txt')


#tree.show()
