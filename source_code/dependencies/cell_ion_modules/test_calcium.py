import numpy as np
import matplotlib.pyplot as plt
import three_state_calcium as ca

t = np.arange(0.0,170,0.5)

params = {"dummy_param": [0]}

ca_class = ca.three_state_ca(params)

fdataCa = open( "calcium_.txt", "w", 0)
ca_value = np.zeros(340)
for i in np.arange(np.shape(t)[0]):
    ca_value[i] = ca_class.calculate_concentrations(1,t[i],fdataCa)

plt.plot(t,ca_value)
plt.show()
fdataCa.close()
