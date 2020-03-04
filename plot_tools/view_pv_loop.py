import matplotlib
import matplotlib.pyplot as plt
import sys
import numpy as np

pv_file_input = sys.argv[1]
time, pressure, volume, temp = np.loadtxt(pv_file_input, skiprows=3, unpack=True)

#time = data[:,0]
#pressure = data[:,1]
#volume = data[:2]

plt.plot(volume,pressure)
plt.show()


