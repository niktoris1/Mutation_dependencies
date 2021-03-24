import matplotlib.pyplot as plt

x = [15, 20, 25, 30]
y = [14.78, 23.32, 24.41, 28.53]

plt.plot(x, y)
plt.suptitle('Example dependence of coalescent estimated rate from birth rate')
plt.show()

x1 = [1.39, 1.74, 1.35, 1.35, 1.70, 1.77, 1.52, 1.60, 1.35, 1.66, 1.53, 1.81, 1.80, 2.02, 1.82, 1.62, 1.63]
y1 = [20.30, 15.36, 24.83, 23.84, 17.54, 16.23, 18.54, 19.70, 15.80, 17.62, 20.98, 19.43, 17.38, 13.94, 19.80, 13.45, 19.68]

#plt.plot(x1, y1, 'bo')
#plt.suptitle('Example dependence of coalescent estimated rate from total time passed')
#plt.show()