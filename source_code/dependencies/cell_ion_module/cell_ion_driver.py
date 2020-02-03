# So far these are only needed to test standalone
import sys
import numpy as np



## Class for cell ion models
#
#calculate calcium concentration (and other concentrations depending on the model) and ion voltage
class cell_ion_driver():

    def __init__(self,params):

        # Specify model to be ran
        temp = params["model"][0]

        base_dir = "cell_ion_module."
        model_name = base_dir + temp
        #model_name = "three_state_calcium"

        # Import the model
        self.model = __import__(model_name)
        #help(self.model)
