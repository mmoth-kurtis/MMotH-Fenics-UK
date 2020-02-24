import numpy as np


## Make this a class to store info so it's not calculated for each particle?

class positionGenerator:

    def init(self, num_particles, gen_string, init_dict):

        self.num_particles = num_particles
        self.gen_string = gen_string
        self.num_vars = len(init_dict.keys())

        # initialize things needed for certain generator algorithms
        if self.gen_string[0] == "uniform":
            self.refinement = int(float(num_particles)**(1.0/num_vars))


    def uniform_generator(j, init_dict, num_particles):
        print "in uniform generator now \n"
        print init_dict
        # init_dict values are lists with entries bounds, position, velocity


        # Need to work out what to do with the remainder from refinement. Should
        # be less than n, so place these points on the corners of the hypercube?


        # Need an index, and a counter for each variable
        self.index_counter_dict = {}
        for var in init_dict.keys():
            self.index_counter_dict[var] = [0, 0]

        # Create arrays of possible values for each variable
        self.possible_value_dictionary = {}
        for var in init_dict.keys():

            # Each variable has its own bounds
            lb = float(init_dict[var][0][0])
            ub = float(init_dict[var][0][1])

            # and its own interval to separate points
            interval = (ub - lb)/(self.refinement+1)

            # temporary storage list for the possible values for this variable
            temp = []

            for i in range(0,self.refinement):
                temp[i] = lb + i*interval

            # Assign this list to the full tensor of possible values
            self.possible_value_dictionary[var] = temp

        # Now we should have a "num_vars" dimensional dictionary, each value has
        # self.refinement number of elements

        # Assign the appropriate value to the j-th particle
        # j should be the sum of refinement^index, for index ranging from 0 to dimensionality

        # check variable counter to see if we've iterated over it

        # for vars in init_dict.keys()
        #   if checksum less than j
        #       

            # set variable value to one of the possible values based on the index
            init_dict[var][1] = self.possible_value_dictionary[var][[self.index_counter_dict][var][0]]

            # update the index and counter dictionary


        return init_dict


    def generate_initial_position(i, initial_dictionary):
        ## function to parse which generating function to use
        #
        # for now, just going to have the uniformly distributed case
        if self.gen_string[0] == "uniform":

            return_dictionary = self.uniform_generator(i, initial_dictionary)

        return return_dictionary
