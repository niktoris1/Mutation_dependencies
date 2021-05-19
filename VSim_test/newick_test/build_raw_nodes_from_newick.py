import re

def nodes_from_newick_file(filename):

    with open(filename, 'r') as file:
        data = file.read().replace('\n', '')

    print('STARTED PARSING')

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

    print('PARSING ENDED')

    return raw_nodes




#tree.show()
