import numpy as np
import matplotlib.pyplot as plt


tstep = np.arange(0.0,100.0,0.5)
step_size =0.5
print tstep[0]
scalefactor = np.zeros(200)
scalefactor_2 = np.zeros(200)
pops = np.zeros(200)
for i in np.arange(200):
    #print i
    #scalefactor[i] = 1.55*np.sin((3.14*(np.exp(-((((tstep[i]+20.00*step_size)/850.)-.01)*23.5)**2.35)))*2*3.14/360.)
    #scalefactor_2[i] = .93*np.sin((3.14*(np.exp(-((((tstep[i]+20.00*step_size)/850.)-.01)*23.5)**2.35)))*2*3.14/360.)
    pops[i] = 0.01*(1+np.sin((1.0/16.0)*tstep[i]+80.2))

print scalefactor[0:2]
fig = plt.figure()
#plt.plot(tstep,scalefactor)
#plt.plot(tstep,scalefactor_2)
plt.plot(tstep,pops)
plt.show()
