import matplotlib.pyplot as plt
from VGsim_test import susceptibleArray, infectiousArray, all_times
#from VGsim._BirthDeath import LogDynamics

#data = simulation.LogDynamics(1000) # the argument is simply number of timepoints
# Data looks like [timestamps, arrayYfh]

#overall_infected = []
#distinct_lineages = []
#for i in range(len(data[0])):
#    overall_infected.append(sum(data[1][i][0]))
#    distinct_lineages.append(ls.DistinctLineages(data[0][i]))
#plt.plot(data[0], overall_infected)
#plt.plot(data[0], distinct_lineages)
#plt.show()

susceptibleArray = list(filter(lambda a: a != 0, susceptibleArray))
infectiousArray = list(filter(lambda a: a != 0, infectiousArray))
all_times = [0] + list(filter(lambda a: a != 0, all_times))

plt.plot(all_times, susceptibleArray) #TODO - for some reason the
plt.plot(all_times, infectiousArray)

plt.show()