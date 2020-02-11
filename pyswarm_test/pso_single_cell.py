import numpy as np
from pyswarm import pso


# Define objective to be minimized
# Starting with minimizing difference between FE force and target

def muscle_force_difference(x, *args):
    target_force = 30000

    # Need a general way to modify the input file
    # For now, hard coding k3 as the input parameters
    
