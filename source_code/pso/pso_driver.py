import sys
import json
sys.path.append("/home/fenics/shared/source_code/pso")
sys.path.append("/home/fenics/shared/source_code/fenics_cases")
import fenicsParticle
import matplotlib.pyplot as plt
import numpy as np

## Call this from fenics_driver? Pass in parameters?
def particle_swarm_optimization(pso_params,sim_params,file_inputs,output_params,passive_params,hs_params,cell_ion_params,monodomain_params,windkessel_params,objFunction):

    params = [sim_params, file_inputs, output_params, passive_params,hs_params,cell_ion_params,monodomain_params,windkessel_params]
    # Assign optimization parameters
    max_iterations = pso_params["max_iterations"][0]
    num_particles = pso_params["num_particles"][0]
    dimensionality = pso_params["num_variables"][0]
    variables_and_bounds = pso_params["variables_and_bounds"] # dictionary of variables and their bounds
    w = pso_params["w"][0] # Velocity inertial parameter
    c1 = pso_params["c1"][0] # Weight for particle position
    c2 = pso_params["c2"][0] # Weight for best swarm position

    # initalize individual particle error holder
    particle_errors = np.zeros((max_iterations,num_particles))

    # hard coding this for now
    target_force = pso_params["singlecell_target_force"][0]

    # defining a dictionary to store optimization history information
    opt_history = { \
    "global_error_history": [] \
    }


    # Initialize best global error and best global position
    best_global_error = -1
    #best_global_position = []
    best_global_position = variables_and_bounds


    # Initialize swarm
    swarm = []
    for i in range(0,num_particles):
        x0 = []

        # Initialize position for particle
        #for j in range(0,dimensionality):

            #x0.append(np.random.uniform([bounds[j][0],bounds[j][1]]))
            # Need to just work with dictionary of parameters?
            # For key in dictionary, if it's in "variable and bounds"
            # modify them in relevant dictionary of input parameters,
            # and pass these to the fenicsParticle class?

        # Right now, x0 initialization is within the particle class
        temp_particle = fenicsParticle.fenicsParticle(dimensionality, target_force, variables_and_bounds, params)
        swarm.append(temp_particle)


    ## Begin optimization
    iter = 0

    fig, ax = plt.subplots()

    xdata = np.arange(max_iterations)

    while iter < max_iterations:

        print "iteration # " + str(iter)
        print "current minimum objective " + str(best_global_error)
        print best_global_position

        ## Iterate through swarm
        for j in range(0,num_particles):
            swarm[j].evaluate_objective(objFunction)

            # Update the error for this particle for visualizing on the fly
            particle_errors[iter,j] = swarm[j].current_error

            #plot particle errors?
            ax.plot(xdata[0:iter],particle_errors[0:iter])
            ax.set(xlabel='Iteration', ylabel='Particle Error', title='Error for each iteration')
            ax.grid()
            fig.savefig(output_params["output_path"][0]+"test.png")

            #plt.show()
            #plt.get_current_fig_manager().full_screen_toggle()
            #plt.show()
            #fig.canvas.flush_events()
            #time.sleep(1)


            # Determine if this is the best particle
            if swarm[j].current_error < best_global_error or best_global_error == -1:
                best_global_position = swarm[j].working_dict
                best_global_error = float(swarm[j].current_error)
                output_dictionary = swarm[j].output_dict

        opt_history["global_error_history"].append(best_global_error)

        ## Update velocities and position (start new loop to make sure
        # we are using the best global position after objective is evaluated for
        # entire swarm)
        for j in range(0,num_particles):
            swarm[j].update_particle_velocity(w,c1,c2,best_global_position)
            swarm[j].update_particle_position()

        iter += 1

    ## Map components of position back to input keys, get "final_inputs"
    # The best location will have been "found", don't want to re-run.

    return(best_global_position,output_dictionary,opt_history)
