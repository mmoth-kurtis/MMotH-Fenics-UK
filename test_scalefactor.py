import numpy as np
import matplotlib.pyplot as plt


tstep = np.arange(0.0,100.0,0.5)
step_size =0.5
print tstep[0]
scalefactor = np.zeros(200)
for i in np.arange(200):
    #print i
    scalefactor[i] = 0.93*np.sin((3.14*(np.exp(-((((tstep[i]+16*step_size)/850.)-.01)*23.5)**2.35)))*2*3.14/360.)
    print scalefactor[i]

fig = plt.figure()
plt.plot(tstep,scalefactor)
plt.show()
