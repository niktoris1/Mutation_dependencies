cdef class Mutations:
    cdef:
        vector[Py_ssize_t] nodeId, AS, DS, site
        vector[double] time

    def __init__(self):#AS = ancestral state, DS = derived state
        pass

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cdef void AddMutation(self, Py_ssize_t nodeId, Py_ssize_t haplotype, Py_ssize_t newHaplotype, double time):
        cdef:
            Py_ssize_t ASDSdigit4, site, digit4
        ASDSdigit4 = int(abs(newHaplotype - haplotype))
        site = 0
        while ASDSdigit4 >= 4:
            ASDSdigit4 = ASDSdigit4 / 4
            site += 1
        digit4 = int(4**site)
        self.nodeId.push_back(nodeId)
        self.DS.push_back(int(floor(newHaplotype/digit4) % 4))
        self.AS.push_back(int(floor(haplotype/digit4) % 4))
        self.site.push_back(int(site))
        self.time.push_back(time)

        def GetData(self):
            return {'node': self.nodeId, 'AS': self.AS, 'DS': self.DS, 'site': self.site, 'time': self.time}
            # print("MutType, AS, DS: ", site, self.AS[self.AS.size()-1], self.DS[self.DS.size()-1])

cdef class Migrations:
    cdef:
        vector[Py_ssize_t] nodeId, oldPop, newPop
        vector[double] time

    def __init__(self):
        pass

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cdef void AddMigration(self, Py_ssize_t nodeId, double time, Py_ssize_t oldPop, Py_ssize_t newPop):
        self.nodeId.push_back(nodeId)
        self.time.push_back(time)
        self.oldPop.push_back(oldPop)
        self.newPop.push_back(newPop)

class Population:
    def __init__(self, size=1000000, contactDensity=1.0):
        self.size = size
        self.contactDensity = contactDensity

class Lockdown:
    def __init__(self, conDenAfterLD=0.1, startLD=2, endLD=1):
        self.conDenAfterLD = conDenAfterLD
        self.startLD = startLD
        self.endLD = endLD


cdef class PopulationModel:
    cdef:
        Py_ssize_t globalInfectious

        long[::1] sizes, totalSusceptible, totalInfectious, lockdownON
        long[:,::1] susceptible

        double[::1] contactDensity, contactDensityBeforeLockdown, contactDensityAfterLockdown, startLD, endLD, samplingMultiplier

    def __init__(self, populations, susceptible_num, lockdownModel=None, samplingMultiplier=None):
        sizePop = len(populations)

        self.sizes = np.zeros(sizePop, dtype=np.int64)
        for i in range(sizePop):
            self.sizes[i] = populations[i].size

        self.totalSusceptible = np.zeros(sizePop, dtype=np.int64)
        self.totalInfectious = np.zeros(sizePop, dtype=np.int64)
        self.globalInfectious = 0

        self.susceptible = np.zeros((sizePop, susceptible_num), dtype=np.int64)
        for i in range(sizePop):
            self.susceptible[i, 0] = populations[i].size
            self.totalSusceptible[i] = populations[i].size

        self.contactDensity = np.zeros(sizePop, dtype=float)
        for i in range(sizePop):
            self.contactDensity[i] = populations[i].contactDensity

        self.contactDensityBeforeLockdown = np.zeros(sizePop, dtype=float)
        self.contactDensityAfterLockdown = np.zeros(sizePop, dtype=float)
        self.startLD = np.zeros(sizePop, dtype=float)
        self.endLD = np.zeros(sizePop, dtype=float)
        if lockdownModel != None:
            for i in range(sizePop):
                self.contactDensityBeforeLockdown[i] = populations[i].contactDensity
                self.contactDensityAfterLockdown[i] = lockdownModel[i].conDenAfterLD
                self.startLD[i] = lockdownModel[i].startLD*self.sizes[i]/100.0
                self.endLD[i] = lockdownModel[i].endLD*self.sizes[i]/100.0
        else:
            for i in range(sizePop):
                self.contactDensityBeforeLockdown[i] = 0
                self.contactDensityAfterLockdown[i] = 0
                self.startLD[i] = 1.01*self.sizes[i]
                self.endLD[i] = 1.0*self.sizes[i]
        self.lockdownON = np.zeros(sizePop, dtype=np.int64)
        if samplingMultiplier != None:
            self.samplingMultiplier = np.asarray(samplingMultiplier)
        else:
            self.samplingMultiplier = np.ones(sizePop)


    @cython.boundscheck(False)
    @cython.wraparound(False)
    cdef inline void NewInfection(self, Py_ssize_t pi, Py_ssize_t si):
        self.susceptible[pi, si] -= 1
        self.totalSusceptible[pi] -= 1
        self.totalInfectious[pi] += 1
        self.globalInfectious += 1

    @cython.boundscheck(False)
    @cython.wraparound(False)
    cdef inline void NewRecovery(self, Py_ssize_t pi, Py_ssize_t si):
        self.susceptible[pi, si] += 1
        self.totalSusceptible[pi] += 1
        self.totalInfectious[pi] -= 1
        self.globalInfectious -= 1
