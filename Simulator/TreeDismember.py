import numpy as np

class TreeDismemberIO:
    def gettdm(self):
        return TreeDismember(self)

    def getrtopo(self, topo):
        rtopo = -np.ones((len(topo),2), dtype=int)
        for i in range(len(topo)-1):
            if rtopo[topo[i]][0] == -1:
                rtopo[topo[i]][0] = i
            else:
                rtopo[topo[i]][1] = i
        return rtopo

    def __init__(self, genealogy, times, mutations, **kwargs):
        self.topo = genealogy
        self.rtopo = self.getrtopo(genealogy)
        self.mutations_raw = mutations
        self.times = times
        self.Tm = []
        self.T = []

class TreeDismember:
    def __init__(self, TreeIO):
        self.topo = TreeIO.topo
        self.mutations_raw = TreeIO.mutations_raw
        self.Mtor = np.zeros(len(self.topo))
        self.rtopo = TreeIO.rtopo
        self.times = TreeIO.times

        #debug
        self.nodes_copies_funct = [0]*len(self.topo)
        self.nodes_copies_neutral = [0]*len(self.topo)
        self.time_events_funct = {x: 0 for x in self.times}
        self.time_events_neutral = {x: 0 for x in self.times}

    def getmut(self, topo, mut, allele): #no site
        # returns array of 1 and -1
        # 1 if DS = 3, -1 - if AS = 3
        al = {'A' : 0, 'T' : 1, 'C' : 2, 'G' : 3}
        nod = mut[0]
        AS = mut[1]
        DS = mut[3]
        #time = mut[4]

        su_AS = -np.ones(len(topo), dtype=int)
        su_DS = -np.ones(len(topo), dtype=int)
        M = np.zeros(len(topo), dtype=int)
        M_count = np.zeros(len(topo), dtype=int)

        for i in range(len(nod)-1, -1, -1):
            if su_AS[nod[i]] == -1:
                su_AS[nod[i]] = AS[i]
            su_DS[nod[i]] = DS[i]
            if su_DS[nod[i]] == su_AS[nod[i]]:
                su_AS[nod[i]] = -1
                su_DS[nod[i]] = -1

        for i in range(len(topo)):
            if su_DS[i] == al[allele]:
                M[i] = 1
            elif su_AS[i] == al[allele]:
                M[i] = -1
            M_count[i] += 1
        return M, M_count

    def debug(self, tprint=True):
        mcount = 0
        for m in self.M:
            if m == 1:
                mcount+=1
        bmcount = 0
        for bm in self.M:
            if bm == -1:
                bmcount+=1

        f_dtf = 0
        f_dtn = 0
        for x, y in self.time_events_funct.items():
            if y > 1:
                f_dtf += 1
        for x, y in self.time_events_neutral.items():
            if y > 1:
                f_dtn += 1

        for x, y in self.time_events_neutral.items():
            if y > 1:
                f_dtn += 1

        f_dnf = 0
        f_dnn = 0
        for x in self.nodes_copies_funct:
            if x > 1:
                f_dnf += 1
        for x in self.nodes_copies_neutral:
            if x > 1:
                f_dnn += 1

        debug_dict = {'Mutated trees (all)' : len(self.trees_funct),
                      'Simple trees (all):' : len(self.trees_neutral),
                      'number of not unique nodes(in funct trees):' : f_dnf,
                      'number of not unique nodes(in neutral trees):' : f_dnn,
                      'Tree len (num of vertexes):' : len(self.topo),
                      'number of start mutations (any>'+self.allele+'):' : mcount,
                      'number of back mutations ('+self.allele+'>any):' : bmcount,
                      'number of mutated tables: ' : len(self.event_table_funct),
                      'number of simple tables: ' : len(self.event_table_neutral),
                      'number of not unique events(in funct tables):' : f_dtf,
                      'number of not unique events(in neutral tables):' : f_dtn}
                      #'sample fraction table:' : self.sample_fraction_table}

        if tprint:
            for name, value in debug_dict.items():
                print(name, value)
        return debug_dict

    def what_type(self, root):
        al = {'A' : 0, 'T' : 1, 'C' : 2, 'G' : 3}
        allele_type = 0
        if self.M[root] == 1:
            allele_type = 1
        elif self.M[root] == -1:
            allele_type = -1
        else:
            #print('root without mut:',root)
            if self.mutations_raw[1] != [] and self.mutations_raw[1][-1] == al[self.allele]:
                allele_type = 1
            else:
                allele_type = -1
        return allele_type


    def Dismember(self, allele='G'):
        rtopo = self.rtopo
        self.M, self.M_count = self.getmut(self.topo, self.mutations_raw, allele)
        self.allele = allele
        mut = self.M
        trees_funct = []
        trees_neutral = []
        node = len(rtopo)-1

        NewRoots = [] # корни поддеревьев
        NewRoots.append(node)

        while NewRoots != []:   #для каждого корня ищем то, что отрезать
            Subtree = [] #для формирования поддерева
            Subtree_is_sample = []
            S = [] # для обхода
            S.append(NewRoots[0])
            mut_type = self.what_type(NewRoots[0])
            while S != []:   #ищем корни для деревьев, которые отрежем
                node = S.pop(-1)
                if mut[node] == 0 or node == NewRoots[0]: #если не нашли мутацию или оказались в рассматриваемом корне
                    Subtree.append(node) #добавляем вершину
                    Subtree_is_sample.append(1) #превентивно - семплирование

                    if mut_type == 1:  # если мутация функциональная, то пометим вершину как мут
                        self.Mtor[node] = 1
                        self.nodes_copies_funct[node] += 1 #debug
                    else:
                        self.nodes_copies_neutral[node] += 1 #debug
                        self.Mtor[node] = 0

                    if (rtopo[node][1] != -1 or rtopo[node][0] != -1):
                        if rtopo[node][1] != -1:
                            S.append(rtopo[node][1])  #идем дальше

                        if rtopo[node][0] != -1:
                            S.append(rtopo[node][0])
                        Subtree_is_sample[-1] = 0 #меняем семплирование на коал. так как есть ветвление
                else:   # иначе добавляем вершину к новым корням
                    NewRoots.append(node)
                    Subtree.append(node) # добавить вершину-обрезок
                    Subtree_is_sample.append(1) #обрезок = семплирование

            NewRoots.pop(0)
            if mut_type == 1:  # если мутация функциональная, то добавляем поддерево к мутировавшим
                Subtree_zip = [Subtree, Subtree_is_sample]
                trees_funct.append([Subtree, Subtree_is_sample])
            else:    # если мутация - откат, или ее нет (главный корень), то добавляем поддерево к обычным
                Subtree_zip = [Subtree, Subtree_is_sample]
                trees_neutral.append([Subtree, Subtree_is_sample])
        self.trees_funct = trees_funct
        self.trees_neutral = trees_neutral
        return trees_funct, trees_neutral # [наборы вершин поддеревьев], [наборы вершин поддеревьев]


    def getEventTable(self, ignore_single_node=True):
        event_table_funct = []
        for tree_i in range(len(self.trees_funct)):
            if len(self.trees_funct[tree_i][0]) != 1 or not ignore_single_node: #костыль, отметаем деревья из 1 вершины
                TreeTable = []
                times_ind = {}
                funct_subtree = self.trees_funct[tree_i][0]
                funct_subtree_is_sample = self.trees_funct[tree_i][1]
                for node, is_sample in zip(funct_subtree, funct_subtree_is_sample):
                    time_i = times_ind.get(self.times[node], -1)
                    if  time_i == -1:
                        times_ind[self.times[node]] = len(TreeTable)
                        TreeTable.append([self.times[node], 0, 0])
                        time_i = times_ind[self.times[node]]
                        #debug
                        self.time_events_funct[self.times[node]] += 1
                        #debug
                    if is_sample:
                        TreeTable[time_i][1] += 1
                    else :
                        TreeTable[time_i][2] += 1
                event_table_funct.append(TreeTable)
        event_table_neutral = []
        for tree_i in range(len(self.trees_neutral)):
            if len(self.trees_neutral[tree_i][0]) != 1  or not ignore_single_node: #костыль, отметаем деревья из 1 вершины
                TreeTable = []
                times_ind = {}
                neutral_subtree = self.trees_neutral[tree_i][0]
                neutral_subtree_is_sample = self.trees_neutral[tree_i][1]
                for node, is_sample in zip(neutral_subtree, neutral_subtree_is_sample):
                    time_i = times_ind.get(self.times[node], -1)
                    if time_i == -1:
                        times_ind[self.times[node]] = len(TreeTable)
                        TreeTable.append([self.times[node], 0, 0])
                        time_i = times_ind[self.times[node]]
                        #debug
                        self.time_events_neutral[self.times[node]] += 1
                        #debug
                    if is_sample:
                        TreeTable[time_i][1] += 1
                    else :
                        TreeTable[time_i][2] += 1
                event_table_neutral.append(TreeTable)
        self.event_table_funct = event_table_funct
        self.event_table_neutral = event_table_neutral
        return event_table_funct, event_table_neutral
        #event_table_funct (массив) - таблицы для каждого поддерева с мутацией
        #event_table_neutral (массив) - таблицы для каждого поддерева без мутации
        #вид таблицы: [[время, кол-во семплов, кол-во коал.], ...]


    def getSampleFracTable(self, tb):   #нужно предоставить массив моментов времени, которыми разбиваем время на интервалы
        sample_fraction_table = [[bracket, 0] for bracket in tb[:-1]]    #tb = ([моменты времени (левая и правая границы - явно), которые разбивают время])
        Mtor_times = np.array([self.times, self.Mtor])
        Mtor_times = np.array([self.times, self.Mtor]).transpose()
        Mtor_times = sorted(Mtor_times, key=lambda x: x[0])


        time_bin = 0
        I1 = 0
        I2 = 0
        for node in range(len(self.topo)):
            while time_bin < len(tb) - 1 and Mtor_times[node][0] > tb[time_bin + 1]: #(classic)
                time_bin += 1
                I1 = 0
                I2 = 0

            if Mtor_times[node][0] < tb[time_bin] or (Mtor_times[node][0] > tb[time_bin] and time_bin == len(tb) - 1):
                pass

            else:
                if self.rtopo[node][0]==self.rtopo[node][1]:
                    if Mtor_times[node][1]:
                        I1 += 1
                    else:
                        I2 += 1
                if I2 == 0:
                    sample_fraction_table[time_bin][1] = -1
                else:
                    sample_fraction_table[time_bin][1] = I1/I2

        self.sample_fraction_table = sample_fraction_table
        return sample_fraction_table
        #sample_fraction_table - таблица с отношениями количеств семплов
        #вид таблицы: [[time_bin, fraction], ...]; fraction = -1 если I1 / 0;
        #time_bin -  левая граница временного интервала
