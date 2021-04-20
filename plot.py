import matplotlib.pyplot as plt

myFiler = open("report_data", "r")
Colist = str(myFiler.read()).split('\n')

tree_size_list = Colist[0].split(' ')[1:]
coal_rate_list = Colist[1].split(' ')[1:]
program_time_list = Colist[2].split(' ')[1:]
time_passed_list = Colist[3].split(' ')[1:]

tree_size_list = [float(item) for item in tree_size_list]
coal_rate_list = [float(item) for item in coal_rate_list]
program_time_list = [float(item) for item in program_time_list]
time_passed_list = [float(item) for item in time_passed_list]


plt.plot(tree_size_list, coal_rate_list, 'bo')
plt.show()
