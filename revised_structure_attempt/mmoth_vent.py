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
from methods.circulatory_module import circulatory_module as cm
from methods.update_boundary_conditions import update_boundary_conditions
import recode_dictionary
import json
import timeit


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
    t = np.linspace(0,sim_duration,no_of_time_steps)
    save_cell_output = sim_params["save_cell_output"][0] # myosim output
    save_visual_output = sim_params["save_visual_output"][0] # paraview files for visualization
    output_path = sim_params["output_path"][0]

    # assign amount of random variation in f0 (cube and cylinder simulations, 0 means normal alignment)
    gaussian_width = sim_params["fiber_randomness"][0]

    # growth parameters
    kroon_time_constant = growth_params["kroon_time_constant"][0]
    ecc_growth_rate = growth_params["ecc_growth_rate"][0]
    set_point = growth_params["passive_set_point"][0]
    k_myo_damp = Constant(growth_params["k_myo_damp"][0])

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
    if sim_geometry == "cylinder" or sim_geometry == "unit_cube" or sim_geometry == "box_mesh" or sim_geometry == "gmesh_cylinder":
        no_of_int_points = 4 * np.shape(mesh.cells())[0]
        deg = 2
        ds = dolfin.ds(subdomain_data = facetboundaries)
        fx_rxn = np.zeros((no_of_time_steps))

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
        # Want to visualize fiber directions through simulation
        fiber_file = File(output_path + "f0_vectors.pvd")
        mesh_file = File(output_path + "mesh_growth.pvd")
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

        calcium_ds = pd.DataFrame(np.zeros(no_of_int_points),index=None)
        calcium_ds = calcium_ds.transpose()

        # myosim populations
        dumped_populations_ds = pd.DataFrame(np.zeros((no_of_int_points,n_array_length)))
        dumped_populations = np.zeros((no_of_int_points,n_array_length))

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
    rxn_force = np.zeros(no_of_time_steps)
    delta_hsl_array = np.zeros(no_of_int_points)
    traction_switch_flag = 0


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

    expressions = {
        "u_D":u_D,
        "Press":Press
    }

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

    Telem2 = TensorElement("Quadrature", mesh.ufl_cell(), degree=deg, shape=2*(3,), quad_scheme='default')
    Telem2._quad_scheme = 'default'
    for e in Telem2.sub_elements():
    	e._quad_scheme = 'default'

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

    # Tensor function space
    TF = TensorFunctionSpace(mesh, 'CG', 2)
    TFQuad = FunctionSpace(mesh, Telem2)


    if sim_geometry == "cylinder" or sim_geometry == "unit_cube" or sim_geometry == "box_mesh" or sim_geometry == "gmesh_cylinder":
        W = FunctionSpace(mesh, MixedElement([Velem,Qelem]))
        x_dofs = W.sub(0).sub(0).dofmap().dofs() # will use this for x rxn forces later
        print "assigned W"
    else:
        W = FunctionSpace(mesh, MixedElement([Velem,Qelem,Relem,VRelem]))

    print str(W)
    # Function space to differentiate fibrous tissue from contractile for fiber simulations
    # Can use thise to mark gauss points according to a user defined law in "assign_heterogeneous_params"
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

    d = u.ufl_domain().geometric_dimension()
    I = Identity(d)
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
    File(output_path + "c param.pvd") << project(c_param,FunctionSpace(mesh,"DG",0))
#-------------------------------------------------------------------------------
#           Save initial values
#-------------------------------------------------------------------------------

    #save initial f0, s0, n0, hsl0
    hsl_temp = project(hsl0,FunctionSpace(mesh,'DG',1))
    hsl_temp.rename("hsl_temp","half-sarcomere length")
    hsl_file << hsl_temp

    File(output_path + "fiber.pvd") << project(f0, VectorFunctionSpace(mesh, "DG", 0))
    File(output_path + "sheet.pvd") << project(s0, VectorFunctionSpace(mesh, "DG", 0))
    File(output_path + "sheet-normal.pvd") << project(n0, VectorFunctionSpace(mesh, "DG", 0))

#-------------------------------------------------------------------------------
#           Initialize the solver and forms parameters, continuum tensors
#-------------------------------------------------------------------------------

    # Create growth tensor. Initialized as identity
    M1ij = project(as_tensor(f0[m]*f0[k], (m,k)),TFQuad)
    M2ij = project(as_tensor(s0[m]*s0[k], (m,k)),TFQuad)
    M3ij = project(as_tensor(n0[m]*n0[k], (m,k)),TFQuad)
    #M1ij = Function(TFQuad)
    #print "m1ij shape = " + str(np.shape(M1ij.vector()))
    #M2ij = Function(TFQuad)
    #M3ij = Function(TFQuad)

    Theta1 = Function(Quad)
    Theta1.vector()[:] = 1.0

    Theta2 = Function(Quad)
    Theta2.vector()[:] = 1.0

    # Based on the material coordinates, we can define different Growth Tensor Construct

    Fg = project(Theta1*(M1ij) +  Theta2*(M2ij + M3ij),TFQuad)
    #print "fg"
    #print str(Fg.vector().get_local())
    #Fg = as_tensor([[1.,0.,0.],[0.,1.,0.],[0.,0.,1.]])
    print "Fg"
    #print Fg[0][0]


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
             "Kappa":Constant(1e5),
             "growth_tensor": Fg}

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
    test_FS = FunctionSpace(mesh,Velem)
    Fmat2 = Function(TF)
    #d = u.ufl_domain().geometric_dimension()
    #Fmat2= Identity(d) + grad(u)
    print "shape of Fmat:"
    print np.shape(Fmat)
    print "Shape of u"
    print np.shape(u)

    # Get right cauchy stretch tensor
    #Cmat = (Fmat.T*Fmat)
    Cmat = uflforms.Cmat()

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
    hsl_diff_from_reference = (hsl_old - hsl0)/hsl0
    print "hsl diff from ref:"
    #print hsl_diff_from_reference.vector().get_local()[:]
    pseudo_alpha = pseudo_old*(1.-(k_myo_damp*(hsl_diff_from_reference)))
    #
    #print "pseudo_alpha:"
    #print project(pseudo_alpha,Quad).vector().get_local()[:]
    alpha_f = sqrt(dot(f0, Cmat*f0)) # actual stretch based on deformation gradient
    #print "alpha_f:"
    #print project(alpha_f,Quad).vector().get_local()[:]
    hsl = pseudo_alpha*alpha_f*hsl0
    #print "hsl:"
    #print project(hsl,Quad).vector().get_local()[:]
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

    #Get passive stress tensor, magnitude of myofiber passive stress
    PK2_passive,Sff = uflforms.stress(hsl)

    # need p_f_array for myosim
    temp_DG = project(Sff, FunctionSpace(mesh, "DG", 1), form_compiler_parameters={"representation":"uflacs"})
    p_f = interpolate(temp_DG, Quad)
    p_f_array = p_f.vector().get_local()[:]

    # calculate magnitude of passive stress for guccione_fiber, guccione_transverse, and guccione_shear
    """if save_cell_output:

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
        pgs_array = pgs.vector().get_local()[:]"""

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
        print "using F4"
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
        print "not using F4"
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


            temp_DG = project(Sff, FunctionSpace(mesh, "DG", 1), form_compiler_parameters={"representation":"uflacs"})
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

                """temp_DG_2 = project(Pg_fiber, FunctionSpace(mesh, "DG", 1), form_compiler_parameters={"representation":"uflacs"})
                pgf = interpolate(temp_DG_2, Quad)
                pgf_array = pgf.vector().get_local()[:]
                temp_DG_3 = project(Pg_transverse, FunctionSpace(mesh, "DG", 1), form_compiler_parameters={"representation":"uflacs"})
                pgt = interpolate(temp_DG_3, Quad)
                pgt_array = pgt.vector().get_local()[:]
                temp_DG_4 = project(Pg_shear, FunctionSpace(mesh, "DG", 1), form_compiler_parameters={"representation":"uflacs"})
                pgs = interpolate(temp_DG_4, Quad)
                pgs_array = pgs.vector().get_local()[:]"""


            if(MPI.rank(comm) == 0):

                print >>fdataPV, 0.0, p_cav*0.0075 , 0.0, 0.0, V_cav, 0.0, 0.0, 0.0

                if save_visual_output:
                    displacement_file << w.sub(0)
                    pk1temp = project(inner(f0,Pactive*f0),FunctionSpace(mesh,'DG',0))
                    pk1temp.rename("pk1temp","active_stress")
                    active_stress_file << pk1temp
                    hsl_temp = project(hsl,FunctionSpace(mesh,'DG',1))
                    hsl_temp.rename("hsl_temp","half-sarcomere length")
                    hsl_file << hsl_temp

            print("cavity-vol = ", LVCavityvol.vol)
            print("p_cav = ", uflforms.LVcavitypressure())

#-------------------------------------------------------------------------------
#           Time Loop
#-------------------------------------------------------------------------------
    # Initialize half-sarcomere class. Methods used to calculate cross-bridges
    # at gauss points
    hs = half_sarcomere.half_sarcomere(hs_params,1)

    # Initialize cell ion module
    cell_ion = cell_ion_driver.cell_ion_driver(cell_ion_params)

    # Initialize calcium concentration from cell_ion module
    calcium[0] = cell_ion.calculate_concentrations(0,0)

    # Load in circulatory module
    if (sim_geometry == "ventricle") or (sim_geometry == "ellipsoid"):
        circ_model = cm.circ_module(windkessel_params)


    for l in np.arange(no_of_time_steps):
        tic = timeit.default_timer()

        print "Time step number " + str(l)
        if (sim_geometry == "ventricle") or (sim_geometry == "ellipsoid"):

            # Update circulatory model
            p_cav = uflforms.LVcavitypressure()
            V_cav = uflforms.LVcavityvol()
            circ_dict = circ_model.update_compartments(p_cav,V_cav,sim_timestep)
            LVCavityvol.vol = V_cav
            LVcav_array[counter] = V_cav
            Pcav_array[counter] = p_cav*0.0075

            # Now print out volumes, pressures, calcium
            if(MPI.rank(comm) == 0):
                print >>fdataPV, tstep, p_cav*0.0075 , Part*.0075, Pven*.0075, V_cav, V_ven, V_art, calcium[counter]

        # update calcium
        calcium[l] = cell_ion.calculate_concentrations(sim_timestep,l)

        # Quick hack
        if l == 0:
            overlap_counter = 1
        else:
            overlap_counter = l

        # At each gauss point, solve for cross-bridge distributions using myosim
        print "calling myosim"
        """for mm in np.arange(no_of_int_points):
            temp_overlap[mm], y_interp[mm*n_array_length:(mm+1)*n_array_length], y_vec_array_new[mm*n_array_length:(mm+1)*n_array_length] = implement.update_simulation(hs, sim_timestep, delta_hsl_array[mm], hsl_array[mm], y_vec_array[mm*n_array_length:(mm+1)*n_array_length], p_f_array[mm], cb_f_array[mm], calcium[l], n_array_length, t,hs_params_list[mm])
            temp_flux_dict, temp_rate_dict = implement.return_rates_fenics(hs)
            j3_fluxes[mm,l] = sum(temp_flux_dict["J3"])
            j4_fluxes[mm,l] = sum(temp_flux_dict["J4"])"""

        if save_cell_output:
            for  i in range(no_of_int_points):
                for j in range(n_array_length):
                    # saving the interpolated populations. These match up with active
                    # stress from previous timestep
                    dumped_populations[i, j] = y_interp[i * n_array_length + j]

        # Update the populations
        y_vec_array = y_vec_array_new # for Myosim

        # Update the population function for fenics
        y_vec.vector()[:] = y_vec_array # for PDE

        # Update the array for myosim
        hsl_array_old = hsl_array

        # Update the hsl_old function for fenics
        hsl_old.vector()[:] = hsl_array_old[:]
        #print "hsl_old before solve"
        #print project(hsl_old,Quad).vector().get_local()[:]

        # including call to nsolver commented out to show it can be used
        #solver.solvenonlinear()

        print "calling Newton Solver"
        # solve for displacement to satisfy balance of linear momentum
        solve(Ftotal == 0, w, bcs, J = Jac, form_compiler_parameters={"representation":"uflacs"},solver_parameters={"newton_solver":{"relative_tolerance":1e-8},"newton_solver":{"maximum_iterations":50},"newton_solver":{"absolute_tolerance":1e-8}})
        f_proj =project(Fmat,TensorFunctionSpace(mesh,"DG",1),form_compiler_parameters={"representation":"uflacs"})
        print "saving def gradient"
        file = File(output_path+'f_proj.pvd')
        file << f_proj
        #Pg_fs_shear = inner(s0,Pg*f0)
        #shear_file = File(output_path+'P_fs_shear.pvd')
        #shear_file << project(Pg_fs_shear,FunctionSpace(mesh,"DG",1), form_compiler_parameters={"representation":"uflacs"})
        #F_matrix = PETScMatrix()
        #f_assembled = assemble(Fmat,tensor=F_matrix)
        #print f_proj.vector().get_local()
        print "guccione passive stress"
        print project(PK2_passive,TensorFunctionSpace(mesh,"DG",1),form_compiler_parameters={"representation":"uflacs"}).vector().get_local()
        #print str(project(p,FunctionSpace(mesh,"CG",1)).vector().get_local())
        #print project(temp,TensorFunctionSpace(mesh,"DG",1),form_compiler_parameters={"representation":"uflacs"}).vector().get_local()
        #print "J"
        #print project(jforms,FunctionSpace(mesh,"DG",0)).vector().get_local()

        #print "f assembly"
        #print f_assembled
        #print "Shape of u"
        #print np.shape(u)

        # Update functions and arrays
        cb_f_array[:] = project(cb_force, Quad).vector().get_local()[:]
        #print "hsl_old after solve"
        #print project(hsl_old,Quad).vector().get_local()[:]
        hsl_old.vector()[:] = project(hsl, Quad).vector().get_local()[:] # for PDE
        pseudo_old.vector()[:] = project(pseudo_alpha, Quad).vector().get_local()[:]
        hsl_array = project(hsl, Quad).vector().get_local()[:]           # for Myosim
        #delta_hsl_array = project(sqrt(dot(f0, Cmat*f0))*hsl0, Quad).vector().get_local()[:] - hsl_array_old # for Myosim
        delta_hsl_array = project(delta_hsl, Quad).vector().get_local()[:]
        #pseudo_old = pseudo_alpha
        #pseudo_old.vector().get_local()[:] = project(pseudo_alpha,Quad).vector().get_local()[:]
        #print "pseudo old"
        #print project(pseudo_old,Quad).vector().get_local()[:]
        #print "pseudo alpha:"
        #print project(pseudo_alpha,Quad).vector().get_local()[:]
        #print "real alpha"
        #print project(alpha_f,Quad).vector().get_local()[:]
        #print "hsl diff from ref "
        #print project(hsl_diff_from_reference,Quad).vector().get_local()[:]
        #print "hsl_old"
        #print project(hsl_old,Quad).vector().get_local()[:]
        #print "hsl0"
        #print project(hsl0,Quad).vector().get_local()[:]
        #print "hsl"
        #print project(hsl,Quad).vector().get_local()[:]

        # For growth
        """temp_DG_2 = project(Pg_fiber, FunctionSpace(mesh, "DG", 1), form_compiler_parameters={"representation":"uflacs"})
        pgf = interpolate(temp_DG_2, Quad)"""

        temp_DG = project(Sff, FunctionSpace(mesh, "DG", 1), form_compiler_parameters={"representation":"uflacs"})
        p_f = interpolate(temp_DG, Quad)
        p_f_array = p_f.vector().get_local()[:]

        for ii in range(np.shape(hsl_array)[0]):
            if p_f_array[ii] < 0.0:
                p_f_array[ii] = 0.0

        # Kroon update fiber orientation?
        if kroon_time_constant != 0.0:
            # update fiber orientations
            print "updating fiber orientation"

        print "updating boundary conditions"
        # Update boundary conditions/expressions (need to include general displacements and tractions)
        bc_update_dict = update_boundary_conditions.update_bcs(bcs,sim_geometry,Ftotal,geo_options,sim_protocol,expressions,t[l],traction_switch_flag,x_dofs)
        bcs = bc_update_dict["bcs"]
        traction_switch_flag = bc_update_dict["traction_switch_flag"]
        rxn_force[l] = bc_update_dict["rxn_force"]
        u_D = bc_update_dict["expr"]["u_D"]
        Press = bc_update_dict["expr"]["Press"]


        # Save visualization info
        if save_visual_output:
            displacement_file << w.sub(0)
            #pk1temp = project(inner(f0,Pactive*f0),FunctionSpace(mesh,'CG',1))
            #pk1temp.rename("pk1temp","active_stress")
            #active_stress_file << pk1temp
            hsl_temp = project(hsl,FunctionSpace(mesh,'DG',1))
            hsl_temp.rename("hsl_temp","half-sarcomere length")
            hsl_file << hsl_temp
            np.save(output_path + 'fx',rxn_force)
            f0_temp = project(f0, VectorFunctionSpace(mesh, "DG", 0))
            f0_temp.rename('f0','f0')
            fiber_file << f0_temp
            #File(output_path + "fiber.pvd") << project(f0, VectorFunctionSpace(mesh, "DG", 0))


        # Save cell info
        tic_save_cell = timeit.default_timer()
        if save_cell_output:

            """temp_DG_2 = project(Pg_fiber, FunctionSpace(mesh, "DG", 1), form_compiler_parameters={"representation":"uflacs"})
            pgf = interpolate(temp_DG_2, Quad)
            pgf_array = pgf.vector().get_local()[:]
            temp_DG_3 = project(Pg_transverse, FunctionSpace(mesh, "DG", 1), form_compiler_parameters={"representation":"uflacs"})
            pgt = interpolate(temp_DG_3, Quad)
            pgt_array = pgt.vector().get_local()[:]
            temp_DG_4 = project(Pg_shear, FunctionSpace(mesh, "DG", 1), form_compiler_parameters={"representation":"uflacs"})
            pgs = interpolate(temp_DG_4, Quad)
            pgs_array = pgs.vector().get_local()[:]"""

            active_stress_ds.iloc[0,:] = cb_f_array[:]
            active_stress_ds.to_csv(output_path + 'active_stress.csv',mode='a',header=False)

            #active_stress_ds = active_stress_ds.transpose()
            hsl_array_ds.iloc[0,:] = hsl_array[:]
            hsl_array_ds.to_csv(output_path + 'half_sarcomere_lengths.csv',mode='a',header=False)

            calcium_ds.iloc[0,:] = calcium[l]
            calcium_ds.to_csv(output_path + 'calcium.csv',mode='a',header=False)

            for i in range(no_of_int_points):
                dumped_populations_ds.iloc[i,:] = dumped_populations[i,:]
            dumped_populations_ds.to_csv(output_path + 'populations.csv',mode='a',header=False)

            #tarray_ds[l] = tarray[l]
            #tarray_ds.to_csv(output_path + 'time.csv',mode='a',header=False)
            np.save(output_path+"time", t) # fix this

            p_f_array_ds.iloc[0,:] = p_f_array[:]
            p_f_array_ds.to_csv(output_path + 'myofiber_passive.csv',mode='a',header=False)

            """pgf_array_ds.iloc[0,:] = pgf_array[:]
            pgf_array_ds.to_csv(output_path + 'gucc_fiber_pstress.csv',mode='a',header=False)

            pgt_array_ds.iloc[0,:] = pgt_array[:]
            pgt_array_ds.to_csv(output_path + 'gucc_trans_pstress.csv',mode='a',header=False)

            pgs_array_ds.iloc[0,:] = pgs_array[:]
            pgs_array_ds.to_csv(output_path + 'gucc_shear_pstress.csv',mode='a',header=False)"""

            temp_overlap_ds.iloc[0,:] = temp_overlap[:]
            temp_overlap_ds.to_csv(output_path + 'overlap.csv',mode='a',header=False)

            #alpha_array_ds.iloc[0,:] = alpha_array[:]
            #alpha_array_ds.to_csv(output_path + 'alpha.csv',mode='a',header=False)

            delta_hsl_array_ds.iloc[0,:] = delta_hsl_array[:]
            delta_hsl_array_ds.to_csv(output_path + 'delta_hsl.csv',mode='a',header=False)

        toc_save_cell = timeit.default_timer() - tic_save_cell
        print "time to save cell info = " + str(toc_save_cell)
        toc = timeit.default_timer() - tic
        print "time loop performance time = " + str(toc)
    if sim_geometry == "ventricle" or sim_geometry == "ellipsoid":
        if(MPI.rank(comm) == 0):
            fdataPV.close()

    # -------------- Attempting growth here --------------------------------

    Fg00 = 1.0
    """for n_grow in np.arange(10):


        #Get passive stress tensors from forms
        PK2_passive= uflforms.stress(hsl)
        Pg_fiber = inner(f0,Pg*f0)


        # For growth
        temp_DG_2 = project(Pg_fiber, FunctionSpace(mesh, "DG", 1), form_compiler_parameters={"representation":"uflacs"})
        pgf = interpolate(temp_DG_2, Quad)
        pgf_array = pgf.vector().get_local()[:]

        temp_DG = project(Sff, FunctionSpace(mesh, "DG", 1), form_compiler_parameters={"representation":"uflacs"})
        p_f = interpolate(temp_DG, Quad)
        p_f_array = p_f.vector().get_local()[:]


        print " deviation from set point = " + str(p_f_array[0]+pgf_array[0]-set_point)
        Fg00 *= 1.+ecc_growth_rate*((p_f_array[0]+pgf_array[0]-set_point)/set_point)
        Theta1.vector()[:] = Fg00
        Theta2.vector()[:] = 1./Fg00
        Fg = project(Theta1*(M1ij) + Theta2*(M2ij + M3ij))
        uflforms.parameters["growth_tensor"] = Fg
        #Fg22 = 1./Fg00
        #M2ij.vector()[:] = Fg11
        #M3ij.vector()[:] = Fg22

        #Fg = as_tensor([[Fg00,0.,0.],[0.,Fg11,0.],[0.,0.,Fg22]])

        print Fg.vector().get_local()[:]
        Fe = uflforms.Fe()
        print Fe
        print "myof passive before solving:"
        print p_f_array[0]
        print "gucc fiber passive before solving:"
        print pgf_array[0]

        #Theta1.vector()[:] = Theta1.vector().array()+ecc_growth_rate*((p_f + pgf-set_point)/set_point)
        #Theta2.vector()[:] = 1./theta1.vector().array()
        solve(Ftotal == 0, w, bcs, J = Jac, form_compiler_parameters={"representation":"uflacs"},solver_parameters={"newton_solver":{"relative_tolerance":1e-8},"newton_solver":{"maximum_iterations":50},"newton_solver":{"absolute_tolerance":1e-8}})
        #Get passive stress tensors from forms
        PK2_passive = uflforms.stress(hsl)
        Pg_fiber = inner(f0,Pg*f0)
        # Update functions and arrays
        cb_f_array[:] = project(cb_force, Quad).vector().get_local()[:]
        hsl_old.vector()[:] = project(hsl, Quad).vector().get_local()[:] # for PDE
        hsl_array = project(hsl, Quad).vector().get_local()[:]           # for Myosim
        delta_hsl_array = project(sqrt(dot(f0, Cmat*f0))*hsl0, Quad).vector().get_local()[:] - hsl_array_old # for Myosim

        # For growth
        temp_DG_2 = project(Pg_fiber, FunctionSpace(mesh, "DG", 1), form_compiler_parameters={"representation":"uflacs"})
        pgf = interpolate(temp_DG_2, Quad)
        pgf_array = pgf.vector().get_local()[:]

        temp_DG = project(Sff, FunctionSpace(mesh, "DG", 1), form_compiler_parameters={"representation":"uflacs"})
        p_f = interpolate(temp_DG, Quad)
        p_f_array = p_f.vector().get_local()[:]
        print "myof passive after solving:"
        print p_f_array[0]
        print "gucc fiber passive after solving:"
        print pgf_array[0]

        for ii in range(np.shape(hsl_array)[0]):
            if p_f_array[ii] < 0.0:
                p_f_array[ii] = 0.0

        if save_cell_output:
            p_f_array_ds.iloc[0,:] = p_f_array[:]
            p_f_array_ds.to_csv(output_path + 'myofiber_passive.csv',mode='a',header=False)

            pgf_array_ds.iloc[0,:] = pgf_array[:]
            pgf_array_ds.to_csv(output_path + 'gucc_fiber_pstress.csv',mode='a',header=False)

            hsl_array_ds.iloc[0,:] = hsl_array[:]
            hsl_array_ds.to_csv(output_path + 'half_sarcomere_lengths.csv',mode='a',header=False)

    ALE.move(mesh, project(u, VectorFunctionSpace(mesh, 'CG', 1)))
    File(output_path + "mesh_grown.pvd") << mesh"""

        #ALE.move(mesh, project(u, VectorFunctionSpace(mesh, 'CG', 1)))
        #File(output_path + "mesh_"+str(n_grow)+".pvd") << mesh
    #print "displacement during growth = " + str(u_D.u_D)
        #displacement_file << w.sub(0)





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
growth_params = input_parameters["growth_and_remodeling"]
#optimization_params = input_parameters["optimization_parameters"]
fenics(sim_params)
