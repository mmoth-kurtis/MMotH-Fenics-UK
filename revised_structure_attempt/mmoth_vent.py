from __future__ import division
import sys
sys.path.append("/home/fenics/shared/source_code/dependencies/")
sys.path.append("/home/fenics/shared/revised_structure_attempt")
import os as os
from dolfin import *
import numpy as np
#from petsc4py import PETSc
#from dependencies import forms
from forms import Forms
from nsolver import NSolver as NSolver
import math
import Python_MyoSim.half_sarcomere.half_sarcomere as half_sarcomere
import Python_MyoSim.half_sarcomere.implement as implement
#import vtk_py
#import mshr
from cell_ion_module import cell_ion_driver
from edgetypebc import *
import pandas as pd
import copy
from mesh_import import import_mesh
import recode_dictionary
import json


# For now, sticking to hieracrchy that this is called by fenics_driver.py
#def fenics(sim_params,file_inputs,output_params,passive_params,hs_params,cell_ion_params,monodomain_params,windkessel_params):
def fenics(sim_params):

    # declare these as indices for things like as_tensor()
    i,j = indices(2)

    ## Assign simulation parameters
    sim_geometry = sim_params["simulation_geometry"][0] #unit_cube, cylinder, ventricle
    geo_options = sim_params["geometry_options"]
    sim_type = sim_params["protocol"]["simulation_type"][0]
    sim_timestep = sim_params["protocol"]["simulation_timestep"][0]
    sim_duration = sim_params["protocol"]["simulation_duration"][0] # can be overwritten in ventricle protocol once number of cycles and heartrate specified
    save_cell_output = sim_params["save_cell_output"][0] # myosim output
    save_visual_output = sim_params["save_visual_output"][0] # paraview files for visualization
    output_path = sim_params["output_path"][0]

#------------------------------------------------------------------------------
#           Mesh Information
#------------------------------------------------------------------------------

    ## Create/load mesh geometry
    mesh = import_mesh(sim_geometry, geo_options)

    File(output_path + '/test_mesh_import.pvd') << mesh

    # define communicator, for running with multiple cores in parallel
    comm = mesh.mpi_comm()

    facetboundaries = MeshFunction('size_t', mesh, mesh.topology().dim()-1)

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
    # Create x interval for cross-bridges
    xx = np.arange(x_bin_min, x_bin_max + x_bin_increment, x_bin_increment)
    # Define number of intervals cross-bridges are defined over
    no_of_x_bins = np.shape(xx)[0]
    # Define the length of the populations vector
    n_array_length = no_of_attached_states * no_of_x_bins + no_of_detached_states + 2
    n_vector_indices = [[0,0], [1,1], [2,2+no_of_x_bins-1]]
    hsl0 = hs_params["initial_hs_length"][0]

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
        #hsl_file = File(output_path + "hsl_mesh.pvd")
        #alpha_file = File(output_path + "alpha_mesh.pvd")

        if sim_geometry == "ventricle":

            # initialize a file for pressures and volumes in windkessel
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
        W = FunctionSpace(mesh, MixedElement([Velem,Qelem,VRelem]))
        x_dofs = W.sub(0).sub(0).dofmap().dofs() # will use this for x rxn forces later
    else:
        W = FunctionSpace(mesh, MixedElement([Velem,Qelem,Relem,VRelem]))

#-------------------------------------------------------------------------------
#           Initialize functions on the above spaces
#-------------------------------------------------------------------------------

    # initialize everything regardless of simulation?

    # fiber, sheet, and sheet-normal functions
    f0 = Function(fiberFS)
    s0 = Function(fiberFS)
    n0 = Function(fiberFS)
    f0_diff = Function(fiberFS)


    # parameters that are heterogeneous here
    c_param = Function(Quad)
    c2_param = Function(Quad)
    c3_param = Function(Quad)

# find the rest of the functions (hsl, etc.)
# assign f0, s0, n0 with outside function?
# assign heterogeneous properties with outside function? e.g. "if cylinder, make ends fiberous,
# if ventricle, could change based on wall location, etc."
# initialize test and trial functions
# initialize solver and forms params
# initialize continuum tensors
# define relationships for pk1, pgt, pgf, etc.
# active stress calculation
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
sim_params = input_parameters["simulation_parameters"]
passive_params = input_parameters["forms_parameters"]["passive_law_parameters"]
hs_params = input_parameters["myosim_parameters"]
cell_ion_params = input_parameters["electrophys_parameters"]["cell_ion_parameters"]
monodomain_params = input_parameters["electrophys_parameters"]["monodomain_parameters"]
windkessel_params = input_parameters["windkessel_parameters"]
optimization_params = input_parameters["optimization_parameters"]
fenics(sim_params)
