import csv

from datetime import date
from parse_data_classes import PhyloDataSimpleImport

tsv_file = open("/Users/LAB-SCG-125/Desktop/SarsData/UShERtree_3thSeptTree_PlusRUS.FullTree.DEEPMOSTSimpleImportsNEW.tsv")
read_tsv = csv.reader(tsv_file, delimiter="\t")

phylo_data = PhyloDataSimpleImport(dictByCountries={})

for node_as_list in read_tsv:
    node_as_node_class = PhyloDataSimpleImport.Node(nodeName=node_as_list[0], country=node_as_list[1],
                              numberOfChildrenInCountry=int(node_as_list[2]), samplingDate=date(int(node_as_list[3][0:4]),
                                                                                                int(node_as_list[3][5:7]),
                                                                                                int(node_as_list[3][8:10])))
    if node_as_node_class.country in phylo_data.dictByCountries.keys():
        phylo_data.dictByCountries[node_as_node_class.country].append(node_as_node_class)
    else:
        phylo_data.dictByCountries[node_as_node_class.country] = [node_as_node_class]

tsv_file.close()