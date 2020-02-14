import numpy as np
import random



## Class for particles in the particle swarm optimization
class fenicsParticle:

## Initialization _____________________________________________________________
    def __init__(self, x0, dimensionality, bounds, target_force):

        # x0 is the initial position for this particle
        # randomly selected within the bounds.

        #bounds consists of {[xi_min, xi_max]} pairs for i = 1,..,dimensionality

        # initialize data holders
        self.position = []
        self.velocity = []
        self.best_ind_position = []
        self.best_ind_error = -1
        self.current_error = -1

        # initialize optimization bounds and dimensionality
        self.dimensionality = dimensionality
        self.bounds = bounds

        # For now, just have a single cell target force
        self.target = target_force

        # initialize position and random velocity
        for i in range(0,self.dimensionality):

            # Generalize for more than one variable
            dim_ub = self.bounds[1]
            dim_lb = self.bounds[0]
            dim_range = dim_ub - dim_lb
            self.velocity.append(random.uniform(-dim_range,dim_range))

            # Need to randomly/uniformly assign positions. Hard coding for now
            self.position.append(np.random.uniform(dim_lb,dim_ub))

# Methods ______________________________________________________________________

    ## Evaluate objective function at current position and check against previous errors
    def evaluate_objective(self, objFunction, params):


        # Unpack parameters to pass to fenics script
        sim_params = params[0]
        file_inputs = params[1]
        output_params = params[2]
        passive_params = params[3]
        hs_params = params[4]
        cell_ion_params = params[5]
        monodomain_params = params[6]
        windkessel_params = params[7]

        # Passing base parameters from input file, will update with new particle values
        # Need to map from "position" to param values
        # For now, hard coding
        hs_params["myofilament_parameters"]["k_3"] = [self.position[0], "text"]

        # Just to test, objective function is this
        temp_output = objFunction.fenics(sim_params,file_inputs,output_params,passive_params,hs_params,cell_ion_params,monodomain_params,windkessel_params)
        predicted_force = temp_output["strarray"][-1]
        self.current_error = self.target-predicted_force

        # Need to map self.position to keyword dictionary for inputs
        # Will need something like
        # self.current_error = objFunction(sim_params,file_inputs,output_params,passive_params,hs_params,cell_ion_params,monodomain_params,windkessel_params)
        # where each dictionary of params is updated with the particle's "position"

        # Check if current position is the best for this particle
        if self.current_error < self.best_ind_error or self.best_ind_error == -1:
            # We have a new best error (or we are initializing)
            self.best_ind_position = self.position
            self.best_ind_error = self.current_error


    #Update particle velocity
    #
    # combination of best particle position, best swarm position, and inertia from previous velocity
    def update_particle_velocity(self, w, c1, c2, best_swarm_position):

        for i in range(0,self.dimensionality):

            r1 = random.random()
            r2 = random.random()

            print "This particle's best position " + str(self.best_ind_position)
            print "This particle's current posittion " + str(self.position)
            print "This particle's best error " + str(self.best_ind_error)
            print "This particle's current error " + str(self.current_error)
            vel_cognitive=c1*r1*(self.best_ind_position[i]-self.position[i])
            vel_social = c2*r2*(best_swarm_position[i]-self.position[i])
            self.velocity[i] = w*self.velocity[i]+vel_cognitive+vel_social


    # Update particle position from updated velocity
    def update_particle_position(self):

        for i in range(0,self.dimensionality):

            self.position[i] = self.position[i] + self.velocity[i]

            # Keep the particle within bounds
            # UPDATE FOR GENERAL NUMBER OF BOUNDS
            if self.position[i] > self.bounds[1]:
                self.position[i] = self.bounds[1]

            if self.position[i] < self.bounds[0]:
                self.position[i] = self.bounds[0]
