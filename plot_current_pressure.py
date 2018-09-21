import sys
import matplotlib.pyplot as plt
import numpy as np

case_name = sys.argv[1]
filename = 'current_pressure.dat'
p = np.loadtxt('./Results/'+case_name+filename)

plt.figure()
plt.plot(p)
plt.draw()
plt.pause(1)
input('enter to close')
plt.close()
