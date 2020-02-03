import numpy as np


## Uses the calcium transient from the three state paper
def display_success():
    print "model successfully imported"

def calculate_concentrations(cycle,time):

    # Don't plan on using this transient much, hard coding some stuff
    t_act = 0.08
    cardiac_period = .17
    t_p = cardiac_period*cycle + t_act+0.01
    fCa = 25
    fCa_2 = 2.5
    calcium_value = 0.0
    pCa = 7

    if time <= t_act:
        pCa = 7
        calcium_value = 10**(-1*pCa)

    elif ((cardiac_period*cycle+t_act) < time) and (time < (t_p)):
        pCa = (time - (cardiac_period*cycle+t_act))/0.02

    elif (time >= t_p):
        pCa = 0.5*np.exp(-np.power((time - t_p)*fCa, fCa_2))

    calcium_value = (0.1 + 1000 * np.sin(3.14*pCa)) * 1E-9

    return calcium_value
