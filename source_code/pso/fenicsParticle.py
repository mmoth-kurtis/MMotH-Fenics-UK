import numpy as np
import random



## Class for particles in the particle swarm optimization
class fenicsParticle:

## Initialization _____________________________________________________________
    def __init__(self, x0, dimensionality, bounds):

        # x0 is the initial position for this particle
        # randomly selected within the bounds.

        #bounds consists of {[xi_min, xi_max]} pairs for i = 1,..,dimensionality

        # initialize data holders
        self.position = []
        self.velocity = []
        self.best_ind_position = []
        self.best_ind_error = []
        self.current_error = -1

        # initialize optimization bounds and dimensionality
        self.dimensionality = dimensionality
        self.bounds = bounds

        # initialize position and random velocity
        for i in range(0,self.dimensionality):

            dim_ub = self.bounds[i][1]
            dim_lb = self.bounds[i][0]
            dim_range = dim_ub - dim_lb
            self.velocity.append(random.uniform(-dim_range,dim_range))

# Methods ______________________________________________________________________

    ## Evaluate objective function at current position and check against previous errors
    def evaluate_objective(self, obj_function):

        self.current_error = obj_function(self.position)

        # Check if current position is the best for this particle
        if self.current_error < self.best_ind_error or self.best_ind_error == -1:
            # We have a new best error (or we are initializing)
            self.best_ind_position = self.position
            self.best_ind_error = self.current_error


    #Update particle velocity
    #
    # combination of best particle position, best swarm position, and inertia from previous velocity
    def update_particle_velocity(w, c1, c2, best_swarm_position):

        for i in range(0,dimensionality):

            r1 = random.random()
            r2 = random.random()

            vel_cognitive=c1*r1*(self.best_ind_position[i]-self.position[i])
            vel_social = c2*r2*(best_swarm_position[i]-self.position[i])
            self.velocity[i] = w*self.velocity[i]+vel_cognitive+vel_social


    # Update particle position from updated velocity
    def update_particle_position(self):

        for i in range(0,self.dimensionality):

            self.position[i] = self.position[i] + self.velocity[i]

            # Keep the particle within bounds
            if self.position[i] > self.bounds[i][1]:
                self.position[i] = bounds[i][1]

            if self.position[i] < self.bounds[i][0]:
                self.position[i] = bounds[i][0]
