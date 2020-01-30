


## Class for cell ion models
#
#calculate calcium concentration (and other concentrations depending on the model) and ion voltage
class cell_ion_module():

    def __init__(self,params):

        # Specify model to be ran
        model = params["model"]
