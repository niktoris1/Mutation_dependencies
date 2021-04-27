from VGsim_test import tree_size, coal_rate, program_time, time_passed

myFiler = open("report_data", "r")
Colist = str(myFiler.read()).split('\n')

if time_passed > 0:
    myFilew = open("report_data", "w")
    myFilew.write(Colist[0] + ' ' + str(tree_size) + '\n')
    myFilew.write(Colist[1] + ' ' + str(coal_rate) + '\n')
    myFilew.write(Colist[2] + ' ' + str(program_time) + '\n')
    myFilew.write(Colist[3] + ' ' + str(time_passed) + '\n')
