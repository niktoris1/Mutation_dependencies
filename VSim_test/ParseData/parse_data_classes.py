class PhyloDataSimpleImport:
    def __init__(self, dictByCountries):
        self.dictByCountries = dictByCountries

    class Node:
        def __init__(self, nodeName, country, numberOfChildrenInCountry, samplingDate):
            self.nodeName = nodeName
            self.country = country
            self.numberOfChildrenInCountry = numberOfChildrenInCountry
            self.samplingDate = samplingDate