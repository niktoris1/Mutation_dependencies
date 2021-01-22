with open('mutations.txt', 'r') as file:
    data = file.read().split('\n')

class MutationNodes:
    def __init__(self, old_nucleotyde, new_nucleotyde, nodes_names, clade_sizes):
        self.old_nucleotyde = old_nucleotyde
        self.new_nucleotyde = new_nucleotyde
        self.nodes_names = nodes_names
        self.clade_sizes = clade_sizes

class AltAlleles:
    def __init__(self, nucleotydes, numbers):
        self.nucleotydes = nucleotydes
        self.numbers = numbers


class Mutation:
    def __init__(self, mutation_name, alt_alleles, parsimony_score, mutation_nodes, flagged_leaves):
        # alt alleles define the nucleotyde variants on this cite
        self.mutation_name = mutation_name
        self.alt_alleles = alt_alleles
        self.parsimony_score = parsimony_score
        self.mutation_nodes = mutation_nodes
        self.flagged_leaves = flagged_leaves


mutations = []

## basically every mutation is defined by its cite in genome. We have an original genome and some changes in this cite, which happen at some point of a tree

for unparsed_mutation in data:
    unparsed_mutation_splited = unparsed_mutation.split('\t')
    alt_alleles = []
    mutation_nodes = []
    flagged_leaves = []
    for parse_element in unparsed_mutation_splited:
        if parse_element is unparsed_mutation_splited[0]:
            mutation_name = parse_element
        if "alt_alleles" in parse_element:
            alt_alleles.append(AltAlleles(parse_element[0], parse_element[-1]))
        if "parsimony_score" in parse_element:
            parsimony_score = parse_element[-1]
        if "mutation_nodes" in parse_element:
            mutation_nodes.append(MutationNodes(parse_element[0], parse_element[2], parse_element.split("=")[1].split(',') , 'Unknown')) # old nycleotype, new nycleotyde, node name, clade size
        if "clade_size" in parse_element:
            for mutation_node in mutation_nodes:
                if mutation_node.old_nucleotyde == parse_element[0] and mutation_node.new_nucleotyde == parse_element[2] and mutation_node.clade_sizes == 'Unknown':
                    mutation_node.clade_sizes = parse_element[-1]
        if "flagged_leaves" in parse_element:
            parse_element = parse_element.split("=")[1]
            for leaves in parse_element.split(','):
                flagged_leaves.append(leaves)

    mutation = Mutation(mutation_name, alt_alleles, parsimony_score, mutation_nodes, flagged_leaves)
    mutations.append(mutation)

















