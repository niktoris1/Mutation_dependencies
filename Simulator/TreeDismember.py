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

    def getmut(self, topo, mut, allele): #no site
        # returns array of 1 and -1
        # 1 if DS = 3, -1 - if AS = 3
        al = {'A' : 0, 'T' : 1, 'C' : 2, 'G' : 3}
        nod = mut[0]
        AS = mut[1]
        DS = mut[3]
        M = np.zeros(len(topo), dtype=int)
        for i in range(len(nod)):
            if DS[i] == al[allele]:
                M[nod[i]] = 1
            elif AS[i] == al[allele]:
                M[nod[i]] = -1
        return M

    def debug(self):
        mcount = 0
        for m in self.M:
            if m == 1:
                mcount+=1
        bmcount = 0
        for bm in self.M:
            if bm == -1:
                bmcount+=1
        print('Mutated trees (all):', len(self.trees_funct))
        print('Simple trees (all):', len(self.trees_neutral))
        print('Tree len (num of vertexes):', len(self.topo))
        print('number of start mutations (any>'+self.allele+'):', mcount)
        print('number of back mutations ('+self.allele+'>any):', bmcount)
        print('number of mutated tables: ', len(self.event_table_funct))
        print('number of simple tables: ', len(self.event_table_neutral))
        print('sample fraction table:', self.sample_fraction_table)

    def Dismember_old(self, allele='G'):
        rtopo = self.rtopo
        self.M = self.getmut(self.topo, self.mutations_raw, allele)
        self.allele = allele
        mut = self.M
        trees_funct = []
        trees_neutral = []
        node = len(rtopo)-1

        NewRoots = [] # корни поддеревьев
        NewRoots.append(node)

        while NewRoots != []:   #для каждого корня ищем то, что отрезать
            Subtree = [] #для формирования поддерева
            InStack = [0]*len(rtopo)
            S = [] # для обхода
            S.append(NewRoots[0])
            while S != []:   #ищем корни для деревьев, которые отрежем
                node = S.pop(-1)
                if mut[node] == 0 or node == NewRoots[0]: #если не нашли мутацию или оказались в рассматриваемом корне
                    f = 0
                    if (rtopo[node][1] != -1 or rtopo[node][0] != -1) and (not InStack[rtopo[node][1]] or not InStack[rtopo[node][0]]):
                        f = 1
                        if rtopo[node][1] != -1 and not InStack[rtopo[node][1]] :
                            S.append(rtopo[node][1])  #идем дальше соблюдая порядок планарного представления
                            InStack[rtopo[node][1]] = 1

                        S.append(node)

                        if rtopo[node][0] != -1 and not InStack[rtopo[node][0]]:
                            S.append(rtopo[node][0])
                            InStack[rtopo[node][0]] = 1
                    if f==0: # добавляем вершину от которой больше не ответвляемся
                        Subtree.append(node)
                        if mut[NewRoots[0]] == 1:  # если мутация функциональная, то пометим вершину как мут
                            self.Mtor[node] = 1
                        else:    # если мутация - откат, или ее нет (главный корень), то добавляем поддерево к обычным
                            self.Mtor[node] = 0
                else:   # иначе добавляем вершину к новым корням
                    NewRoots.append(node)
                    Subtree.append(node) # добавить вершину-обрезок
            if mut[NewRoots.pop(0)] == 1:  # если мутация функциональная, то добавляем поддерево к мутировавшим
                trees_funct.append(Subtree)
            else:    # если мутация - откат, или ее нет (главный корень), то добавляем поддерево к обычным
                trees_neutral.append(Subtree)
        self.trees_funct = trees_funct
        self.trees_neutral = trees_neutral
        return trees_funct, trees_neutral # [наборы вершин поддеревьев], [наборы вершин поддеревьев]


    def Dismember(self, allele='G'):
        rtopo = self.rtopo
        self.M = self.getmut(self.topo, self.mutations_raw, allele)
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
            while S != []:   #ищем корни для деревьев, которые отрежем
                node = S.pop(-1)
                if mut[node] == 0 or node == NewRoots[0]: #если не нашли мутацию или оказались в рассматриваемом корне
                    Subtree.append(node) #добавляем вершину
                    Subtree_is_sample.append(1) #превентивно - семплирование
                    if mut[NewRoots[0]] == 1:  # если мутация функциональная, то пометим вершину как мут
                        self.Mtor[node] = 1
                    else:
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
            if mut[NewRoots.pop(0)] == 1:  # если мутация функциональная, то добавляем поддерево к мутировавшим
                Subtree_zip = [Subtree, Subtree_is_sample]
                trees_funct.append([Subtree, Subtree_is_sample])
            else:    # если мутация - откат, или ее нет (главный корень), то добавляем поддерево к обычным
                Subtree_zip = [Subtree, Subtree_is_sample]
                trees_neutral.append([Subtree, Subtree_is_sample])
        self.trees_funct = trees_funct
        self.trees_neutral = trees_neutral
        return trees_funct, trees_neutral # [наборы вершин поддеревьев], [наборы вершин поддеревьев]


    def getEventTable(self):
        event_table_funct = []
        for tree_i in range(len(self.trees_funct)):
            if len(self.trees_funct[tree_i]) != 1: #костыль, отметаем деревья из 1 вершины
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
                    if is_sample:
                        TreeTable[time_i][1] += 1
                    else :
                        TreeTable[time_i][2] += 1
                event_table_funct.append(TreeTable)
        event_table_neutral = []
        for tree_i in range(len(self.trees_neutral)):
            if len(self.trees_neutral[tree_i]) != 1: #костыль, отметаем деревья из 1 вершины
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
