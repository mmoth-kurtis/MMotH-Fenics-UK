import numpy as np

def init(params):
    print "inside init function"
    model_class = calcium_file_input(params)
    return model_class


class calcium_file_input():

    def __init__(self,params):

        #params is a dictionary defined in the json input files
        # we can be specific here because this module is customizable to user needs
        self.file_path = params["path_to_calcium"][0]
        self.ca = np.load(self.file_path)
        print "inside class init"


    def calculate_concentrations(self,cycle,time,file):

        if time < 0.5:
            calcium_value = 1e-7
        else:
            calcium_value = self.ca[int(2*time),0]
        return calcium_value
