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
tree.create_node(-1, -1)

def add_children(node):
    tree.create_node(node["id"], node["id"], parent=node["parentid"])
    for child_node in node["children"]:
        add_children(child_node)


add_children(raw_nodes)


print('READY')

print(tree.depth())

#tree.show()
