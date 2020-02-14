import sys
import json
sys.path.append("/home/fenics/shared/source_code/pso")
sys.path.append("/home/fenics/shared/source_code/fenics_cases")
import fenicsParticle
import numpy as np

## Call this from fenics_driver? Pass in parameters?
def particle_swarm_optimization(pso_params,sim_params,file_inputs,output_params,passive_params,hs_params,cell_ion_params,monodomain_params,windkessel_params,objFunction):

    params = [sim_params, file_inputs, output_params, passive_params,hs_params,cell_ion_params,monodomain_params,windkessel_params]
    # Assign optimization parameters
    max_iterations = pso_params["max_iterations"][0]
    num_particles = pso_params["num_particles"][0]
    dimensionality = pso_params["num_variables"][0]
    w = pso_params["w"][0] # Velocity inertial parameter
    c1 = pso_params["c1"][0] # Weight for particle position
    c2 = pso_params["c2"][0] # Weight for best swarm position

    # hard coding this for now
    target_force = pso_params["singlecell_target_force"][0]

    # Set up bounds for each input parameters (not general right now)
    bounds = pso_params["variable_and_bounds"][1]
    # Need a mapping between parameters and dimensions of "x"


    # Initialize best global error and best global position
    best_global_error = -1
    best_global_position = []


    # Initialize swarm
    swarm = []
    for i in range(0,num_particles):
        x0 = []

        # Initialize position for particle
        for j in range(0,dimensionality):
            # need special case for one variable?
            x0.append(np.random.uniform([bounds[0],bounds[1]]))
            # Need to just work with dictionary of parameters?
            # For key in dictionary, if it's in "variable and bounds"
            # modify them in relevant dictionary of input parameters,
            # and pass these to the fenicsParticle class?

        temp_particle = fenicsParticle.fenicsParticle(x0, dimensionality, bounds, target_force)
        swarm.append(temp_particle)


    ## Begin optimization
    iter = 0

    while iter < max_iterations:

        print "iteration # " + str(iter)
        print "current minimum objective" + str(best_global_error)

        ## Iterate through swarm
        for j in range(0,num_particles):
            swarm[j].evaluate_objective(objFunction, params)

            # Determine if this is the best particle
            if swarm[j].current_error < best_global_error or best_global_error == -1:
                best_global_position = swarm[j].position
                best_global_error = float(swarm[j].current_error)

        ## Update velocities and position (start new loop to make sure
        # we are using the best global position after objective is evaluated for
        # entire swarm)
        for j in range(0,num_particles):
            swarm[j].update_particle_velocity(w,c1,c2,best_global_position)
            swarm[j].update_particle_position()

        iter += 1

    ## Map components of position back to input keys, get "final_inputs"
    # The best location will have been "found", don't want to re-run.

    return(best_global_position,output_dictionary,best_global_error)
