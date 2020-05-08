import numpy as np

## Just want a sine curve active stress profile with
# magnitude of 30 kPa

# Simulation is going for 110 ms
def calculate_force(time):

    #Time is in ms
    print "time is " + str(time)
    active_force = 15000*(1+np.sin((1.0/16.0)*time + 80.2))
    print "calculating force" + str(active_force)
    return active_force
