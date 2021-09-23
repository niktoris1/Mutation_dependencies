import re
from treelib import Tree, Node

def parse_newick_string(newick):
    tokens = re.finditer(r"([^:;,()\s]*)(?:\s*:\s*([\d.]+)\s*)?([,);])|(\S)", newick+";")

    def recurse(nextid = 0, parentid = -1): # one node
        thisid = nextid;
        children = []

        name, length, delim, ch = next(tokens).groups(0)
        if ch == "(":
            while ch in "(,":
                node, ch, nextid = recurse(nextid+1, thisid)
                children.append(node)
            name, length, delim, ch = next(tokens).groups(0)
        return {"id": thisid, "name": name, "length": float(length) if length else None,
                "parentid": parentid, "children": children}, delim, nextid

    return recurse()[0]

def ParseNodeName(node_name):
    if node_name == '':
        return {'country': 'Unknown', 'sample_name':'Unknown', 'year':'Unknown', 'haplotype':'Unknown', 'date':'Unknown'}
    else:
        result = re.findall(r"[a-zA-Z0-9\-._]+", node_name)
        if len(result) == 7:
            return {'country': result[2], 'sample_name': result[0]+result[1]+result[2], 'year': result[4], 'haplotype': result[5],
                    'date': result[6]}
        if len(result) == 6:
            return {'country': result[3], 'sample_name': result[0]+result[1], 'year': result[2], 'haplotype': result[4],
                    'date': result[5]}
        if len(result) == 5:
            return {'country': result[0], 'sample_name': result[1], 'year': result[2], 'haplotype': result[3], 'date': result[4]}
        if len(result) == 4:
            return {'country': result[0], 'sample_name': result[1], 'year': result[2], 'haplotype': 'Unknown', 'date': result[3]}
        # in assumption, that all node_names of length 4 are the same
        elif len(result) == 3:
            return {'country': 'Unknown', 'sample_name': result[0], 'year': 'Unknown', 'haplotype': result[1],
                    'date': result[2]} # in assumption, that all node_names of length 3 are the same
        elif len(result) == 2:
            return {'country': 'Unknown', 'sample_name': 'Unknown', 'year': 'Unknown', 'haplotype': result[0],
                    'date': result[1]}
        elif len(result) == 1:
            return {'country': 'Unknown', 'sample_name': result[0], 'year': 'Unknown', 'haplotype': 'Unknown',
                    'date': 'Unknown'} # in assumption, that all node_names of length 1 are the same
        else:
            print(result)
           # raise ValueError


def TreeArrayToTreeClass(tree_array):
    tree_class_tree = Tree()

    def AddNode(node_as_dict, tree_class_tree_funct_arg):
        parsed_data = ParseNodeName(node_as_dict['name'])
        node_as_node_class = Node(tag=parsed_data['sample_name'], identifier=node_as_dict['id'],
                                  data={'country': parsed_data['country'], 'haplotype': parsed_data['haplotype'], 'date': parsed_data['date']}) # no data yet
        if node_as_dict['parentid'] == -1:
            tree_class_tree_funct_arg.add_node(node_as_node_class, parent=None)
        else:
            tree_class_tree_funct_arg.add_node(node_as_node_class, parent=node_as_dict['parentid'])
        for child_node in node_as_dict['children']:
            AddNode(child_node, tree_class_tree_funct_arg)

        return tree_class_tree_funct_arg

    result = AddNode(tree_array, tree_class_tree)

    return result

class DataForNode:
    def __init__(self, haplotype_name, date):
        self.haplotype_name = haplotype_name
        self.date = date



