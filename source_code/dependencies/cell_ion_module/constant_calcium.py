import numpy as np


## Uses a constant calcium
def display_success():
    print "model successfully imported"

def calculate_concentrations(cycle,time):

    """# Time is passed in as ms, not seconds
    t = time/1000
    # Don't plan on using this transient much, hard coding some stuff
    t_act = 0.08
    cardiac_period = .17
    t_p = cardiac_period*cycle + t_act+0.01
    fCa = 25
    fCa_2 = 2.5
    calcium_value = 0.0
    pCa = 7

    if t <= t_act:
        pCa = 7
        calcium_value = 10**(-1*pCa)

    elif ((cardiac_period*cycle+t_act) < t) and (t < (t_p)):
        pCa = (t - (cardiac_period*cycle+t_act))/0.02

    elif (t >= t_p):
        pCa = 0.5*np.exp(-np.power((t - t_p)*fCa, fCa_2))

    calcium_value = (0.1 + 1000 * np.sin(3.14*pCa)) * 1E-9"""

    t = time/1000
    t_act = 0.01
    calcium_value = 0.0
    if t > t_act:
        calcium_value = np.power(10.0,-5)

    return calcium_value
