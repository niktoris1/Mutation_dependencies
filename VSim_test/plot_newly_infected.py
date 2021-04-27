import matplotlib.pyplot as plt
from VGsim_test import tree,  newtree, LikelyhoodEstimation, simulation, events
#from VGsim._BirthDeath import LogDynamics

#ls1 = LikelyhoodEstimation(newtree)
#max_time = ls1.events_sequence[-1].event_time

#infected_time = []
#diff_lineages = []
#frequency = 1000
#for i in range(frequency + 1):
#    infected_time.append(i * (max_time / frequency))
#    diff_lineages.append(ls1.DistinctLineages(infected_time[i]))

#plt.plot(infected_time, diff_lineages, 'bo')
#plt.show()

data = simulation.LogDynamics(1000)
#print(data[0])
print(data[1])

for e_id in range(events.ptr):
    pass


y = []
for i in range(len(data[1])):
    y.append(sum(data[1][i][0]))
plt.plot(data[0], y)
plt.show()