from __future__ import division
import sys
sys.path.append("/home/fenics/shared/source_code/dependencies/")
sys.path.append("/home/fenics/shared/revised_structure_attempt")
import os as os
from dolfin import *
import numpy as np
from forms import Forms
from nsolver import NSolver as NSolver
import math
import Python_MyoSim.half_sarcomere.half_sarcomere as half_sarcomere
import Python_MyoSim.half_sarcomere.implement as implement
from cell_ion_module import cell_ion_driver
from edgetypebc import *
import pandas as pd
import copy
from methods import mesh_import
from methods.mesh_import import mesh_import as mesh_import
from methods.assign_initial_hsl import assign_initial_hsl as assign_hsl
from methods.assign_local_coordinate_system import assign_local_coordinate_system as lcs
from methods.assign_heterogeneous_params import assign_heterogeneous_params as assign_params
from methods.set_boundary_conditions import set_bcs as set_bcs
import recode_dictionary
import json


# For now, sticking to hieracrchy that this is called by fenics_driver.py
#def fenics(sim_params,file_inputs,output_params,passive_params,hs_params,cell_ion_params,monodomain_params,windkessel_params):
def fenics(sim_params):

    # declare these as indices for things like as_tensor()
    i,j = indices(2)
    m,k = indices(2)


    ## Assign simulation parameters
    sim_geometry = sim_params["simulation_geometry"][0] #unit_cube, cylinder, ventricle
    if sim_geometry == "unit_cube":
        geo_options = {}
    else:
        geo_options = sim_params["geometry_options"]

    sim_protocol = sim_params["protocol"] # contains simulation dependent options

    sim_timestep = sim_protocol["simulation_timestep"][0]
    sim_duration = sim_protocol["simulation_duration"][0] # can be overwritten in ventricle protocol once number of cycles and heartrate specified
    no_of_time_steps = int(sim_duration/sim_timestep)
    save_cell_output = sim_params["save_cell_output"][0] # myosim output
    save_visual_output = sim_params["save_visual_output"][0] # paraview files for visualization
    output_path = sim_params["output_path"][0]

    # assign amount of random variation in f0 (cube and cylinder simulations, 0 means normal alignment)
    gaussian_width = sim_params["fiber_randomness"][0]

#------------------------------------------------------------------------------
#           Mesh Information
#------------------------------------------------------------------------------

    ## Create/load mesh geometry
    mesh,lv_options = mesh_import.import_mesh(sim_geometry, geo_options)

    File(output_path + '/test_mesh_import.pvd') << mesh

    # define communicator, for running with multiple cores in parallel
    comm = mesh.mpi_comm()

    facetboundaries = MeshFunction('size_t', mesh, mesh.topology().dim()-1)
    edgeboundaries = MeshFunction('size_t', mesh, mesh.topology().dim()-2)


    # from the mesh, define some things
    if sim_geometry == "cylinder" or sim_geometry == "unit_cube" :
        no_of_int_points = 4 * np.shape(mesh.cells())[0]
        deg = 2
        ds = dolfin.ds(subdomain_data = facetboundaries)
    else:
        #ventricle modeling
        deg = 4
        no_of_int_points = 14 * np.shape(mesh.cells())[0]
        #set surface id numbers
        topid = 4
        LVendoid = 2
        epiid = 1

    parameters["form_compiler"]["quadrature_degree"]=deg
    parameters["form_compiler"]["representation"] = "quadrature"
    edgeboundaries = MeshFunction('size_t', mesh, mesh.topology().dim()-2)
    N = FacetNormal (mesh)
    X = SpatialCoordinate (mesh)
    dx = dolfin.dx(mesh,metadata = {"integration_order":2})
    isincomp = True


#------------------------------------------------------------------------------
#           Initialize myosim info for active stress calculation
#------------------------------------------------------------------------------

    filament_compliance_factor = hs_params["myofilament_parameters"]["filament_compliance_factor"][0]
    no_of_states = hs_params["myofilament_parameters"]["num_states"][0]
    no_of_attached_states = hs_params["myofilament_parameters"]["num_attached_states"][0]
    no_of_detached_states = no_of_states-no_of_attached_states
    no_of_transitions = hs_params["myofilament_parameters"]["num_transitions"][0]
    state_attached = hs_params["myofilament_parameters"]["state_attached"][0]
    cb_extensions = hs_params["myofilament_parameters"]["cb_extensions"][0]
    k_cb_multiplier = hs_params["myofilament_parameters"]["k_cb_multiplier"][0]
    k_cb_pos = hs_params["myofilament_parameters"]["k_cb_pos"][0]
    k_cb_neg = hs_params["myofilament_parameters"]["k_cb_neg"][0]
    cb_number_density = hs_params["cb_number_density"][0]
    alpha_value = hs_params["myofilament_parameters"]["alpha"][0]
    x_bin_min = hs_params["myofilament_parameters"]["bin_min"][0]
    x_bin_max = hs_params["myofilament_parameters"]["bin_max"][0]
    x_bin_increment = hs_params["myofilament_parameters"]["bin_width"][0]
    xfiber_fraction = hs_params["myofilament_parameters"]["xfiber_fraction"][0]
    k_myo_damp = hs_params["myofilament_parameters"]["k_myo_damp"][0]

    # Create x interval for cross-bridges
    xx = np.arange(x_bin_min, x_bin_max + x_bin_increment, x_bin_increment)
    # Define number of intervals cross-bridges are defined over
    no_of_x_bins = np.shape(xx)[0]
    # Define the length of the populations vector
    n_array_length = no_of_attached_states * no_of_x_bins + no_of_detached_states + 2
    n_vector_indices = [[0,0], [1,1], [2,2+no_of_x_bins-1]]
    #hsl0 = hs_params["initial_hs_length"][0]

#------------------------------------------------------------------------------
#           Storage data frames
#------------------------------------------------------------------------------

    ## Create files for saving information if needed.
    # cell level info is saved through a pandas command later
    if save_visual_output:
        # Can visualize pretty much anything. For now, just looking at deformation
        # and the active stress magnitude
        displacement_file = File(output_path + "u_disp.pvd")
        active_stress_file = File(output_path + "active_stress_magnitude.pvd")
        hsl_file = File(output_path + "hsl_mesh.pvd")
        #alpha_file = File(output_path + "alpha_mesh.pvd")

        if (sim_geometry == "ventricle") or (sim_geometry == "ellipsoid"):
            # initialize a file for pressures and volumes in windkessel
            if (MPI.rank(comm) == 0):
                fdataPV = open(output_path + "PV_.txt", "w", 0)

    if save_cell_output:

        # storage arrays to print to file
        # active stress
        active_stress_ds = pd.DataFrame(np.zeros(no_of_int_points),index=None)
        active_stress_ds = active_stress_ds.transpose()

        # myosim populations
        dumped_populations_ds = pd.DataFrame(np.zeros((no_of_int_points,n_array_length)))

        # passive stressed in myofiber
        p_f_array_ds = pd.DataFrame(np.zeros(no_of_int_points),index=None)
        p_f_array_ds = p_f_array_ds.transpose()

        # guccione bulk passive stress in fiber direction
        pgf_array_ds = pd.DataFrame(np.zeros(no_of_int_points),index=None)
        pgf_array_ds = pgf_array_ds.transpose()

        # guccione bulk passive stress in transverse directions
        pgt_array_ds = pd.DataFrame(np.zeros(no_of_int_points),index=None)
        pgt_array_ds = pgt_array_ds.transpose()

        # guccione bulk passive shear stress
        pgs_array_ds =pd.DataFrame(np.zeros(no_of_int_points),index=None)
        pgs_array_ds = pgs_array_ds.transpose()

        # thick and thin filament overlaps
        temp_overlap_ds = pd.DataFrame(np.zeros(no_of_int_points),index=None)
        temp_overlap_ds = temp_overlap_ds.transpose()

        # stretch
        alpha_array_ds = pd.DataFrame(np.zeros(no_of_int_points),index=None)
        alpha_array_ds = alpha_array_ds.transpose()

        # half-sarcomere length
        hsl_array_ds =pd.DataFrame(np.zeros(no_of_int_points),index=None)
        hsl_array_ds = hsl_array_ds.transpose()

        # change in half-sarcomere length
        delta_hsl_array_ds = pd.DataFrame(np.zeros(no_of_int_points),index=None)
        delta_hsl_array_ds = delta_hsl_array_ds.transpose()

    # Initialize some data holders that are necessary
    temp_overlap = np.zeros((no_of_int_points))
    y_vec_array_new = np.zeros(((no_of_int_points)*n_array_length))
    j3_fluxes = np.zeros((no_of_int_points,no_of_time_steps))
    j4_fluxes = np.zeros((no_of_int_points,no_of_time_steps))
    y_interp = np.zeros((no_of_int_points)*n_array_length)
    calcium = np.zeros(no_of_time_steps)

    # Put any needed expressions here
    #--------------------------------
    # initialize LV cavity volume, only updated if windkessel is called
    LVCavityvol = Expression(("vol"), vol=0.0, degree=2)

    # displacement boundary expression for end of cell or fiber sims
    u_D = Expression(("u_D"), u_D = 0.0, degree = 2)

    # traction boundary condition for end of cell/fiber, could use this to apply
    # a traction to epicardium or something
    Press = Expression(("P"), P=0.0, degree=0)

    # Sometimes define an expression for active stress
    #cb_force2 = Expression(("f"), f=0.0, degree=0)


#-------------------------------------------------------------------------------
#           Initialize finite elements and function spaces
#-------------------------------------------------------------------------------

    # Vector element at gauss points (for fibers)
    VQuadelem = VectorElement("Quadrature", mesh.ufl_cell(), degree=deg, quad_scheme="default")
    VQuadelem._quad_scheme = 'default'

    # General quadrature element whose points we will evaluate myosim at
    Quadelem = FiniteElement("Quadrature", tetrahedron, degree=deg, quad_scheme="default")
    Quadelem._quad_scheme = 'default'

    # Vector element for displacement
    Velem = VectorElement("CG", mesh.ufl_cell(), 2, quad_scheme="default")
    Velem._quad_scheme = 'default'

    # Quadrature element for pressure
    Qelem = FiniteElement("CG", mesh.ufl_cell(), 1, quad_scheme="default")
    Qelem._quad_scheme = 'default'

    # Real element for rigid body motion boundary condition
    Relem = FiniteElement("Real", mesh.ufl_cell(), 0, quad_scheme="default")
    Relem._quad_scheme = 'default'

    # Mixed element for rigid body motion. One each for x, y displacement. One each for
    # x, y, z rotation
    VRelem = MixedElement([Relem, Relem, Relem, Relem, Relem])

    # ------- Define function spaces on mesh using above elements --------------

    # Quadrature space for information needed at gauss points, such as
    # hsl, cb_force, passive forces, etc.
    Quad = FunctionSpace(mesh, Quadelem)

    # Function space for myosim populations
    Quad_vectorized_Fspace = FunctionSpace(mesh, MixedElement(n_array_length*[Quadelem]))

    # Function space for local coordinate system (fiber, sheet, sheet-normal)
    fiberFS = FunctionSpace(mesh, VQuadelem)

    if sim_geometry == "cylinder" or sim_geometry == "unit_cube":
        W = FunctionSpace(mesh, MixedElement([Velem,Qelem]))
        x_dofs = W.sub(0).sub(0).dofmap().dofs() # will use this for x rxn forces later
    else:
        W = FunctionSpace(mesh, MixedElement([Velem,Qelem,Relem,VRelem]))

    # Function space to differentiate fibrous tissue from contractile for fiber simulations
    marker_space = FunctionSpace(mesh, 'CG', 1)

#-------------------------------------------------------------------------------
#           Initialize functions on the above spaces
#-------------------------------------------------------------------------------

    # fiber, sheet, and sheet-normal functions
    f0 = Function(fiberFS)
    s0 = Function(fiberFS)
    n0 = Function(fiberFS)

    # put these in a dictionary to pass to function for assignment
    coord_params = {
        "f0":f0,
        "s0":s0,
        "n0":n0,
        "fiberFS":fiberFS,
        "marker_space":marker_space,
        "sim_geometry":sim_geometry,
        "mesh":mesh,
        "Quad":Quad,
        "no_of_int_points":no_of_int_points,
        "geo_options":geo_options,
        "facetboundaries":facetboundaries,
        "edgeboundaries":edgeboundaries
    }

    # create lists of dictionaries that hold parameters for each gauss point
    hs_params_list = [{}]*no_of_int_points
    passive_params_list = [{}]*no_of_int_points

    # Must make a deep copy so each item in the list is independent, and not linked
    # to the original paramter dictionary
    for jj in np.arange(np.shape(hs_params_list)[0]):
        hs_params_list[jj] = copy.deepcopy(hs_params)
        passive_params_list[jj] = copy.deepcopy(passive_params)

    # parameters that are heterogeneous declared here as functions
    # Do these need to come from the input file? As part of declaration, "heterogenous = true"?
    # for key in input dictionary:
    #    if dict[key] has value "heterogeneous =  true":
    #        add this param to the list, but how to declare as a function? Hard code for now...
    c_param = Function(Quad)
    c2_param = Function(Quad)
    c3_param = Function(Quad)

    # create heterogeneous function list to be passed in to method "assign_heterogeneous_params"
    # Then
    heterogeneous_fcn_list = [c_param, c2_param, c3_param]

    # functions for the weak form
    w     = Function(W)
    dw    = TrialFunction(W)
    wtest = TestFunction(W)

    if (sim_geometry == "ventricle") or (sim_geometry == "ellipsoid"):
        # Need the pressure function
        du,dp,dpendo,dc11 = TrialFunctions(W)
        (u,p,pendo,c11)   = split(w)
        (v,q,qendo,v11)   = TestFunctions(W)
        ventricle_params  = {
            "lv_volconst_variable": pendo,
            "lv_constrained_vol":LVCavityvol,
            "LVendoid": LVendoid,
            "LVendo_comp": 2,
        }

    else:
        du,dp = TrialFunctions(W)
        (u,p) = split(w)
        (v,q) = TestFunctions(W)

    # Initial and previous timestep half-sarcomere length functions
    hsl0    = Function(Quad)
    hsl_old = Function(Quad)

    # Defining functions to calculate how far hsl is from reference, to be used
    # to get back to reference length over time
    pseudo_alpha = Function(Quad)
    pseudo_old = Function(Quad)
    pseudo_old.vector()[:] = 1.0
    hsl_diff_from_reference = Function(Quad)
    hsl_diff_from_reference.vector()[:] = 0.0

    # Population solution holder for myosim, allows active stress to be calculated
    # and adjusted as the Newton Solver tries new displacements
    y_vec = Function(Quad_vectorized_Fspace)
    y_vec_array = y_vec.vector().get_local()[:]

#-------------------------------------------------------------------------------
#           Assign function values
#-------------------------------------------------------------------------------

    hsl0 = assign_hsl.assign_initial_hsl(lv_options,hs_params,sim_geometry,hsl0)
    f0,s0,n0,geo_options = lcs.assign_local_coordinate_system(lv_options,coord_params,sim_params)

    # Assign the heterogeneous parameters
    heterogeneous_fcn_list,hs_params_list,passive_params_list = assign_params.assign_heterogeneous_params(sim_params,hs_params_list,passive_params_list,geo_options,heterogeneous_fcn_list,no_of_int_points)

#-------------------------------------------------------------------------------
#           Save initial values
#-------------------------------------------------------------------------------

    #save initial f0, s0, n0, hsl0
    hsl_temp = project(hsl0,FunctionSpace(mesh,'DG',1))
    hsl_temp.rename("hsl_temp","half-sarcomere length")
    hsl_file << hsl_temp

    File(output_path + "fiber.pvd") << project(f0, VectorFunctionSpace(mesh, "CG", 1))
    File(output_path + "sheet.pvd") << project(s0, VectorFunctionSpace(mesh, "CG", 1))
    File(output_path + "sheet-normal.pvd") << project(n0, VectorFunctionSpace(mesh, "CG", 1))

#-------------------------------------------------------------------------------
#           Initialize the solver and forms parameters, continuum tensors
#-------------------------------------------------------------------------------

    # parameters for forms file
    params= {"mesh": mesh,
             "facetboundaries": facetboundaries,
             "facet_normal": N,
             "mixedfunctionspace": W,
             "mixedfunction": w,
             "displacement_variable": u,
             "pressure_variable": p,
             "fiber": f0,
             "sheet": s0,
             "sheet-normal": n0,
             "incompressible": isincomp,
             "hsl0": hsl0,
             "Kappa":Constant(1e5)}

    # update passive params because now they are heterogeneous functions
    # need to generalize this?
    params.update(passive_params)
    if (sim_geometry == "ventricle") or (sim_geometry == "ellipsoid"):
        params.update(ventricle_params)
    params["c"] = c_param
    params["c2"] = c2_param
    params["c3"] = c3_param

    # Initialize the forms module
    uflforms = Forms(params)

    #LVCavityvol.vol = uflforms.LVcavityvol() may initialize this later to avoid another if statement here

    # Get deformation gradient
    Fmat = uflforms.Fmat()

    # Get right cauchy stretch tensor
    Cmat = (Fmat.T*Fmat)

    # Get Green strain tensor
    Emat = uflforms.Emat()

    # Polar decomposition stretch tensor, used for Kroon (for now)
    Umat = uflforms.Umat()

    # jacobian of deformation gradient
    J = uflforms.J()

    # facet normal in current config
    n = J*inv(Fmat.T)*N

#-------------------------------------------------------------------------------
#           Initialize boundary conditions
#-------------------------------------------------------------------------------
    # returns a dictionary of bcs and potentially a test_marker_fcn for work loops
    bc_output = set_bcs.set_bcs(sim_geometry,sim_protocol,mesh,W,facetboundaries,u_D)
    bcs = bc_output["bcs"]

#-------------------------------------------------------------------------------
#           Active stress calculation
#-------------------------------------------------------------------------------
    # Start with active stress calculation here to validate mmoth_vent, then try
    # to move it to the forms file

    # Calculate a pseudo stretch, not based on deformation gradient
    hsl_old.vector()[:] = hsl0.vector()[:]
    hsl_diff_from_reference = hsl_old - hsl0
    pseudo_alpha = pseudo_old*(1.-k_myo_damp*hsl_diff_from_reference)
    alpha_f = sqrt(dot(f0, Cmat*f0)) # actual stretch based on deformation gradient
    hsl = pseudo_alpha*alpha_f*hsl0
    delta_hsl = hsl - hsl_old

    cb_force = Constant(0.0)

    y_vec_split = split(y_vec)
    Wp = uflforms.PassiveMatSEF(hsl)

    for jj in range(no_of_states):
        f_holder = Constant(0.0)
        temp_holder = Constant(0.0)

        if state_attached[jj] == 1:
            cb_ext = cb_extensions[jj]

            for kk in range(no_of_x_bins):
                dxx = xx[kk] + delta_hsl * filament_compliance_factor
                n_pop = y_vec_split[n_vector_indices[jj][0] + kk]
                temp_holder = n_pop * k_cb_multiplier[jj] * (dxx + cb_ext) * conditional(gt(dxx + cb_ext,0.0), k_cb_pos, k_cb_neg)
                f_holder = f_holder + temp_holder

            f_holder = f_holder * cb_number_density * 1e-9
            f_holder = f_holder * alpha_value

        cb_force = cb_force + f_holder

    Pactive = cb_force * as_tensor(f0[m]*f0[k], (m,k))+ xfiber_fraction*cb_force * as_tensor(s0[m]*s0[k], (m,k))+ xfiber_fraction*cb_force * as_tensor(n0[m]*n0[k], (m,k))

#-------------------------------------------------------------------------------
#           Now hsl function is initiated, make sure all arrays are initialized
#-------------------------------------------------------------------------------

    #create hsl_array from projection
    hsl_array = project(hsl, Quad).vector().get_local()[:]

    #initialize y_vec_array to put all hads in SRX and all binding sites to off
    for init_counter in range(0,n_array_length * no_of_int_points,n_array_length):
        # Initializing myosin heads in the off state
        y_vec_array[init_counter] = 1
        # Initialize all binding sites to off state
        y_vec_array[init_counter-2] = 1

    #Get passive stress tensors from forms
    Pg, Pff, alpha = uflforms.stress(hsl)

    # need p_f_array for myosim
    temp_DG = project(Pff, FunctionSpace(mesh, "DG", 1), form_compiler_parameters={"representation":"uflacs"})
    p_f = interpolate(temp_DG, Quad)
    p_f_array = p_f.vector().get_local()[:]

    # calculate magnitude of passive stress for guccione_fiber, guccione_transverse, and guccione_shear
    if save_cell_output:

        # Magnitude of bulk passive stress in fiber direction
        Pg_fiber = inner(f0,Pg*f0)
        Pg_transverse = inner(n0,Pg*n0)
        Pg_shear = inner(n0,Pg*f0)

        temp_DG_1 = project(alpha, FunctionSpace(mesh, "DG", 1), form_compiler_parameters={"representation":"uflacs"})
        alphas = interpolate(temp_DG_1, Quad)
        alpha_array = alphas.vector().get_local()[:]

        temp_DG_2 = project(Pg_fiber, FunctionSpace(mesh, "DG", 1), form_compiler_parameters={"representation":"uflacs"})
        pgf = interpolate(temp_DG_2, Quad)
        pgf_array = pgf.vector().get_local()[:]
        temp_DG_3 = project(Pg_transverse, FunctionSpace(mesh, "DG", 1), form_compiler_parameters={"representation":"uflacs"})
        pgt = interpolate(temp_DG_3, Quad)
        pgt_array = pgt.vector().get_local()[:]
        temp_DG_4 = project(Pg_shear, FunctionSpace(mesh, "DG", 1), form_compiler_parameters={"representation":"uflacs"})
        pgs = interpolate(temp_DG_4, Quad)
        pgs_array = pgs.vector().get_local()[:]

    # define cb force array to save
    cb_f_array = project(cb_force, Quad).vector().get_local()[:]

# for LV, initialize LV cavity and LV pressure arrays, and get first values from forms

#-------------------------------------------------------------------------------
#           Weak Form
#-------------------------------------------------------------------------------
    # passive material contribution
    F1 = derivative(Wp, w, wtest)*dx

    # active stress contribution (Pactive is PK1, transform to PK2)
    F2 = inner(Fmat*Pactive, grad(v))*dx

    if (sim_geometry == "ventricle") or (sim_geometry == "ellipsoid"):
        # LV volume increase
        Wvol = uflforms.LVV0constrainedE()
        F3 = derivative(Wvol, w, wtest)

        # constrain rigid body motion
        L4 = inner(as_vector([c11[0], c11[1], 0.0]), u)*dx + \
    	 inner(as_vector([0.0, 0.0, c11[2]]), cross(X, u))*dx + \
    	 inner(as_vector([c11[3], 0.0, 0.0]), cross(X, u))*dx + \
    	 inner(as_vector([0.0, c11[4], 0.0]), cross(X, u))*dx

        F4 = derivative(L4, w, wtest)

        Ftotal = F1 + F2 + F3 + F4

        Jac1 = derivative(F1, w, dw)
        Jac2 = derivative(F2, w, dw)
        Jac3 = derivative(F3, w, dw)
        Jac4 = derivative(F4, w, dw)

        Jac = Jac1 + Jac2 + Jac3 + Jac4

    else:
        F3 = inner(Press*N, v)*ds(2, domain=mesh)

        Ftotal = F1 + F2 - F3

        Jac1 = derivative(F1, w, dw)
        Jac2 = derivative(F2, w, dw)
        Jac3 = derivative(F3, w, dw)
        Jac = Jac1 + Jac2 - Jac3

    # Can use Dr. Lee's Nsolver if solver needs debugging
    solverparams = {"Jacobian": Jac,
                    "F": Ftotal,
                    "w": w,
                    "boundary_conditions": bcs,
                    "Type": 0,
                    "mesh": mesh,
                    "mode": 0
                    }

    solver= NSolver(solverparams)

#-------------------------------------------------------------------------------
#           Initial loading (for ventricle)
#-------------------------------------------------------------------------------

    if (sim_geometry == "ventricle") or (sim_geometry == "ellipsoid"):

        # Load the ventricle from the reference configuration
        LVcav_array = np.zeros(no_of_time_steps)
        LVcav_array[0] = uflforms.LVcavityvol()
        Pcav_array = np.zeros(no_of_time_steps)
        Pcav_array[0] = uflforms.LVcavitypressure()*0.0075

        LVCavityvol.vol = uflforms.LVcavityvol()

        # Calculate the increment to LV volume
        end_diastolic_volume = sim_protocol["initial_end_diastolic_volume"][0]
        total_vol_loading = end_diastolic_volume - LVCavityvol.vol
        volume_increment = total_vol_loading/sim_protocol["reference_loading_steps"][0]

        for lmbda_value in range(0, sim_protocol["reference_loading_steps"][0]):

            print "Diastolic loading step " + str(lmbda_value)

            LVCavityvol.vol += volume_increment

            p_cav = uflforms.LVcavitypressure()
            V_cav = uflforms.LVcavityvol()

            hsl_array_old = hsl_array

            #solver.solvenonlinear()
            solve(Ftotal == 0, w, bcs, J = Jac, form_compiler_parameters={"representation":"uflacs"})

            hsl_array = project(hsl, Quad).vector().get_local()[:]           # for Myosim


            temp_DG = project(Pff, FunctionSpace(mesh, "DG", 1), form_compiler_parameters={"representation":"uflacs"})
            p_f = interpolate(temp_DG, Quad)
            p_f_array = p_f.vector().get_local()[:]

            for ii in range(np.shape(hsl_array)[0]):
                if p_f_array[ii] < 0.0:
                    p_f_array[ii] = 0.0

            delta_hsl_array = hsl_array - hsl_array_old

            if save_cell_output:
                temp_DG_1 = project(alpha, FunctionSpace(mesh, "DG", 1), form_compiler_parameters={"representation":"uflacs"})
                alphas = interpolate(temp_DG_1, Quad)
                alpha_array = alphas.vector().get_local()[:]

                temp_DG_2 = project(Pg_fiber, FunctionSpace(mesh, "DG", 1), form_compiler_parameters={"representation":"uflacs"})
                pgf = interpolate(temp_DG_2, Quad)
                pgf_array = pgf.vector().get_local()[:]
                temp_DG_3 = project(Pg_transverse, FunctionSpace(mesh, "DG", 1), form_compiler_parameters={"representation":"uflacs"})
                pgt = interpolate(temp_DG_3, Quad)
                pgt_array = pgt.vector().get_local()[:]
                temp_DG_4 = project(Pg_shear, FunctionSpace(mesh, "DG", 1), form_compiler_parameters={"representation":"uflacs"})
                pgs = interpolate(temp_DG_4, Quad)
                pgs_array = pgs.vector().get_local()[:]


            if(MPI.rank(comm) == 0):

                print >>fdataPV, 0.0, p_cav*0.0075 , 0.0, 0.0, V_cav, 0.0, 0.0, 0.0

                if save_visual_output:
                    displacement_file << w.sub(0)
                    pk1temp = project(inner(f0,Pactive*f0),FunctionSpace(mesh,'DG',1))
                    pk1temp.rename("pk1temp","active_stress")
                    active_stress_file << pk1temp
                    hsl_temp = project(hsl,FunctionSpace(mesh,'DG',1))
                    hsl_temp.rename("hsl_temp","hsl")
                    hsl_file << hsl_temp

            print("cavity-vol = ", LVCavityvol.vol)
            print("p_cav = ", uflforms.LVcavitypressure())

#-------------------------------------------------------------------------------
#           Time Loop
#-------------------------------------------------------------------------------



# initialize solver and forms params
# weak form
# initialize hs and cell_ion classes










#-------------------------------------------------------------------------------
# for stand-alone testing
input_file = sys.argv[1]
# Load in JSON dictionary
with open(input_file, 'r') as json_input:
  input_parameters = json.load(json_input)

# Convert any unicode values to python strings so they work with some cpp libraries.
recode_dictionary.recode(input_parameters)

# Parse out the different types of parameters.
sim_params = input_parameters["simulation_parameters"]
passive_params = input_parameters["forms_parameters"]["passive_law_parameters"]
hs_params = input_parameters["myosim_parameters"]
cell_ion_params = input_parameters["electrophys_parameters"]["cell_ion_parameters"]
#monodomain_params = input_parameters["electrophys_parameters"]["monodomain_parameters"]
#windkessel_params = input_parameters["windkessel_parameters"]
#optimization_params = input_parameters["optimization_parameters"]
fenics(sim_params)
