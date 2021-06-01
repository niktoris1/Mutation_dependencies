import matplotlib.pyplot as plt
from VGsim_test import tree,  alltree, LikelyhoodEstimation, simulation, event_types, es_ls, ls
#from VGsim._BirthDeath import LogDynamics



data = simulation.LogDynamics(1000) # the argument is simply number of timepoints
# Data looks like [timestamps, arrayYfh]

overall_infected = []
distinct_lineages = []
for i in range(len(data[0])):
    overall_infected.append(sum(data[1][i][0]))
    distinct_lineages.append(ls.DistinctLineages(data[0][i]))
plt.plot(data[0], overall_infected)
plt.plot(data[0], distinct_lineages)
plt.show()