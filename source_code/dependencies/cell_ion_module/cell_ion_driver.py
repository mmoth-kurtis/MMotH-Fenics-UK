# So far these are only needed to test standalone
import sys
import numpy as np
sys.path.append("/home/fenics/shared/source_code/dependencies/cell_ion_module")
#sys.path.append("/Users/charlesmann/Academic/UK/fenics/source_code/dependencies/cell_ion_module")

## Class for cell ion models
#
#calculate calcium concentration (and other concentrations depending on the model) and ion voltage
class cell_ion_driver():

    def __init__(self,params):

        # Specify model to be ran
        model_name = params["model"][0]

        #base_dir = "cell_ion_module"
        #model_name = base_dir + temp
        #model_name = "three_state_calcium"

        #print model_name

        # Import the model
        self.model = __import__(model_name)
        #help(self.model)
