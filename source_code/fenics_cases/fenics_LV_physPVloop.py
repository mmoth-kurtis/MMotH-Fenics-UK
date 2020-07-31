from __future__ import division
import sys
sys.path.append("/home/fenics/shared/Research/FEniCS-Myosim/MMotH-Fenics-UK/source_code/dependencies")
import os as os
from dolfin import *
import numpy as np
from matplotlib import pylab as plt
from petsc4py import PETSc
from forms import Forms
from nsolver import NSolver as NSolver
import math
import json
import recode_dictionary
import load_parameters
import Python_MyoSim.half_sarcomere.half_sarcomere as half_sarcomere
import Python_MyoSim.half_sarcomere.implement as implement
import vtk_py
from cell_ion_module import cell_ion_driver
import homogeneous_stress
from edgetypebc import *
import objgraph as obg

## Fenics simulation function
#
# This simulates some integer number of cardiac cycles.
# As of 07/01/2020, the following is implemented:
# - Combination of Guccione strain-energy function (SEF) and myofiber SEF
# - General MyoSim drives active contraction based on [Ca]2+ from file or 3 state paper
# - 3 compartment Windkessel circulatory system
#
# General layout of this script:
#   - Parse input parameters
#   - Read in mesh, assign fiber, sheet, sheet-normal local coordinate system
#   - Define appropriate finite elements and function spaces
#   - Define variational problem and boundary conditions (active stress calculated here so Newton Iteration can vary it)
#       boundary conditions: zero z-displacement at base, mean x,y displacement = 0
#   - Load ventricle by incrementally increasing cavity volume.
#     No contraction occurs here. Solve for new displacements after each loading step
#   - "Closed Loop Phase"
#       - Solve Windkessel model
#           - Solve for new compartment pressures based on volumes from previous timestep
#           - Calculate blood flow between compartments, update volumes
#       - Solve MyoSim
#           - Calcium updated from cell_ion module before MyoSim call
#           - Interpolate cross-bridges based on solution from previous Newton Iteration
#           - Solve ODEs to get new populations
#       - Solve variational problem to get new displacements
#
# Parameters
# ----------
# @param[in] sim_params Dictionary for setting up the simulation
# @param[in] file_inputs Dictionary containing all file information, now namely just mesh casename
# @param[in] output_params Dictionary to specify output path
# @param[in] passive_params Dictionary for passive material law parameters
# @param[in] hs_params Dictionary for all MyoSim parameters
# @param[in] cell_ion_params Dictionary for cell ion module parameters
# @param[in] monodomain_params Dictionary for monodomain parameters (not implemented)
# @param[in] windkessel_params Dictionary for all Windkessel parameters
# @param[out]
def fenics(sim_params,file_inputs,output_params,passive_params,hs_params,cell_ion_params,monodomain_params,windkessel_params):
    global i
    global j

    # We don't do pressure control simulations, probably will get rid of this.
    ispressurectrl = False

    #------------------## Load in all information and set up simulation #-----------------
    ## Assign input/output parameters
    output_path = output_params["output_path"][0]
    casename = file_inputs["casename"][0]

    # Assign parameters for Windkessel
    # will be moving circulatory to its own module and pass in dictionary
    # similar to cell_ion module
    Cao = windkessel_params["Cao"][0]
    Cven = windkessel_params["Cven"][0]
    Vart0 = windkessel_params["Vart0"][0]
    Vven0 = windkessel_params["Vven0"][0]
    Rao = windkessel_params["Rao"][0]
    Rven = windkessel_params["Rven"][0]
    Rper = windkessel_params["Rper"][0]
    V_ven = windkessel_params["V_ven"][0]
    V_art = windkessel_params["V_art"][0]

    # --------  Assign parameters for active force calculation  ----------------
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
    hsl_min_threshold = hs_params["myofilament_parameters"]["passive_l_slack"][0]
    hsl_max_threshold = hs_params["myofilament_parameters"]["hsl_max_threshold"][0]
    xfiber_fraction = hs_params["myofilament_parameters"]["xfiber_fraction"][0]

    ## ---------  Set up information for active force calculation --------------

    # Create x interval for cross-bridges
    xx = np.arange(x_bin_min, x_bin_max + x_bin_increment, x_bin_increment)

    # Define number of intervals cross-bridges are defined over
    no_of_x_bins = np.shape(xx)[0]

    # Define the length of the populations vector
    n_array_length = no_of_attached_states * no_of_x_bins + no_of_detached_states + 2 # +2 for binding sites on/off

    # Need to work out a general way to set this based on the scheme
    n_vector_indices = [[0,0], [1,1], [2,2+no_of_x_bins-1]]

    #------------  Start setting up simulation ---------------------------------
    sim_duration = sim_params["sim_duration"][0]
    save_output = sim_params["save_output"][0]
    step_size = sim_params["sim_timestep"][0]
    loading_number = sim_params["loading_number"][0]
    if sim_params["sim_geometry"][0] == "ventricle" or sim_params["sim_geometry"][0] == "ventricle_lclee_2" or sim_params["sim_geometry"][0] == "ventricle_physloop":
        # For ventricle for now, specify number of cardiac cycles
        cycles = sim_params["sim_type"][1]
        meshfilename = sim_params["sim_type"][2]
    # Cardiac cycle length and number of cycles will be general
    # For now, just including this info in the input file
    BCL = sim_duration # ms

    hsl0 = hs_params["initial_hs_length"][0] # this is now set when creating mesh
    no_of_time_steps = int(cycles*BCL/step_size)
    no_of_cell_time_steps = int(BCL/step_size)

    deg = 4
    parameters["form_compiler"]["quadrature_degree"]=deg
    parameters["form_compiler"]["representation"] = "quadrature"

    # Clear out any old results files
    os.system("rm " + output_path + "*.pvd")
    os.system("rm " + output_path + "*.vtu")

    #--------------- Load in mesh -------------------------------------
    mesh = Mesh()
    f = HDF5File(mpi_comm_world(), meshfilename, 'r')
    f.read(mesh, casename, False)

    if casename == "ellipsoidal":
        #loading_number = 25;
        ugrid = vtk_py.convertXMLMeshToUGrid(mesh)
        ugrid = vtk_py.rotateUGrid(ugrid, sx=0.11, sy=0.11, sz=0.11)
        mesh = vtk_py.convertUGridToXMLMesh(ugrid)

        #don't need to do the vtk_py mesh stuff
    else: #assuming we are using a patient specific mesh
        ugrid = vtk_py.convertXMLMeshToUGrid(mesh)
        ugrid = vtk_py.rotateUGrid(ugrid, sx=0.1, sy=0.1, sz=0.1)
        mesh = vtk_py.convertUGridToXMLMesh(ugrid)

    no_of_int_points = 14 * np.shape(mesh.cells())[0]

    facetboundaries = MeshFunction("size_t", mesh, 2)
    edgeboundaries = MeshFunction("size_t", mesh, 1)

    # ---------  Initialize finite elements  -----------------------------------

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

    # Mixed function space for displacement, pressure, rigid body constraint
    if(ispressurectrl):
        W = FunctionSpace(mesh, MixedElement([Velem,Qelem,VRelem]))
    else:
        W = FunctionSpace(mesh, MixedElement([Velem,Qelem,Relem,VRelem]))

    # V isn't used? Could define function spaces V: Velem, Q:Qelem,VR: VRelem, then W = V*W*VR
    # but below W is explicitly defined using the elements?
    # could define these once and use them for all projections
    #V = VectorFunctionSpace(mesh, 'CG', 2)
    #TF = TensorFunctionSpace(mesh, 'DG', 1)
    #Q = FunctionSpace(mesh,'CG',1)

    f0 = Function(fiberFS)
    s0 = Function(fiberFS)
    n0 = Function(fiberFS)
    hsl0_transmural = Function(Quad)

    f.read(hsl0_transmural, casename+"/"+"hsl0")
    f.read(f0, casename+"/"+"eF")
    f.read(s0, casename+"/"+"eS")
    f.read(n0, casename+"/"+"eN")

    f.read(facetboundaries, casename+"/"+"facetboundaries")
    f.read(edgeboundaries, casename+"/"+"edgeboundaries")

    f.close()
    File(output_path + "facetboundaries.pvd") << facetboundaries
    File(output_path + "edgeboundaries.pvd") << edgeboundaries
    topid = 4
    LVendoid = 2
    epiid = 1
    ##############################################################################

    comm = mesh.mpi_comm()

    File(output_path + "fiber.pvd") << project(f0, VectorFunctionSpace(mesh, "CG", 1))
    File(output_path + "sheet.pvd") << project(s0, VectorFunctionSpace(mesh, "CG", 1))
    File(output_path + "sheet-normal.pvd") << project(n0, VectorFunctionSpace(mesh, "CG", 1))

    ##############################################################################
    isincomp = True#False
    N = FacetNormal (mesh)
    X = SpatialCoordinate (mesh)

    Press = Expression(("P"), P=0.0, degree=0)
    Kspring = Constant(100)
    LVCavityvol = Expression(("vol"), vol=0.0, degree=2)




    #################################################################################



    #bctop1 = DirichletBC(W.sub(0).sub(0), s0(), facetboundaries, topid)
    #bctop2 = DirichletBC(W.sub(0).sub(1), W.sub(0).sub(0), facetboundaries, topid)
    bctop = DirichletBC(W.sub(0).sub(2), Expression(("0.0"), degree = 2), facetboundaries, topid)
    #bctop = DirichletBC(W.sub(0), Expression((("0.0", "0.0", "0.0")), degree = 2), facetboundaries, topid)

    ######### Define an expression for epicardial edge radial displacement field ######################
    class MyExpr(Expression):

        def __init__(self, param, **kwargs):
        	self.param = param

    	def eval(self, values, x):
            nvec = self.param["nvec"]
            rdir = np.array(nvec(x))/np.linalg.norm(nvec(x))

            scalefactor = self.param["scalefactor"]


            values[0] = scalefactor.a*rdir[0]
            values[1] = scalefactor.a*rdir[1]
            values[2] = 0.0


        def value_shape(self):
            return (3,)



    #s0CG.set_allow_extrapolation(True)
    scalefactor_epi = Expression('a', a=0.0, degree=1)
    scalefactor_endo = Expression('a', a=0.0, degree=1)
    #scalefactor_hsl = Expression('a+b*sin(t*c/d+c/2)', a=1.09, b=0.09, c=3.14, d=100, t=0.0, degree=1)
    """param_endo={"nvec": s0CG,
           "scalefactor": scalefactor_endo
            };
    param_epi={"nvec": s0CG,
           "scalefactor": scalefactor_epi
            };"""

    #radexp = MyExpr(param=param, degree=1)
    #endo_exp = MyExpr(param=param_endo, degree=1)
    #epi_exp = MyExpr(param=param_epi, degree=1)
    ###################################################################################################

    # edgeboundaries 1 is epi, 2 is endo
    endoring = pick_endoring_bc(method="cpp")(edgeboundaries, 2)
    epiring = pick_endoring_bc(method="cpp")(edgeboundaries, 1)
    # CKM trying to have radial expansion for entire top
    #endoring = pick_endoring_bc(method="cpp")(facetboundaries, topid)
    #bcedge = DirichletBC(W.sub(0), Expression(("0.0", "0.0", "0.0"), degree = 0), endoring, method="pointwise")
    #bcedge_endo = DirichletBC(W.sub(0), endo_exp, endoring, method="pointwise")
    #bcedge_epi = DirichletBC(W.sub(0), epi_exp, epiring, method="pointwise")
    #bcs = [bctop, bcedge_endo, bcedge_epi]
    bcs = [bctop]

    w = Function(W)
    dw = TrialFunction(W)
    wtest = TestFunction(W)

    if(ispressurectrl):
        du,dp,dc11 = TrialFunctions(W)
    	(u,p,c11) = split(w)
    	(v,q,v11) = TestFunctions(W)
    else:
        du,dp,dpendo,dc11 = TrialFunctions(W)
      	(u,p, pendo,c11) = split(w)
        #(u,p, pendo,c11,lm11) = w.split(True)
      	(v,q, qendo,v11) = TestFunctions(W)


    if(ispressurectrl):
    	pendo = []
    #dt = Expression(("dt"), dt=0.0, degree=1)

    params= {"mesh": mesh,
             "facetboundaries": facetboundaries,
             "facet_normal": N,
             "mixedfunctionspace": W,
             "mixedfunction": w,
             "displacement_variable": u,
             "pressure_variable": p,
             "lv_volconst_variable": pendo,
             "lv_constrained_vol":LVCavityvol,
             "LVendoid": LVendoid,
             "LVendo_comp": 2,
             "fiber": f0,
             "sheet": s0,
             "sheet-normal": n0,
             "incompressible": isincomp,
             "Kappa":Constant(1e5)}

    # Update params from loaded in parameters from json file
    params.update(passive_params)

    uflforms = Forms(params)

    Fmat = uflforms.Fmat()
    Cmat = (Fmat.T*Fmat)
    Emat = uflforms.Emat()
    J = uflforms.J()

    n = J*inv(Fmat.T)*N
    dx = dolfin.dx(mesh,metadata = {"integration_order":2})

    #Ematrix = project(Emat, TF) not used?
    Wp = uflforms.PassiveMatSEF()

    #Active force calculation------------------------------------------------------
    y_vec = Function(Quad_vectorized_Fspace)
    hsl = sqrt(dot(f0, Cmat*f0))*hsl0_transmural
    #hsl = project(sqrt(dot(f0, Cmat*f0))*hsl0_transmural, Quad)
    # Forcing hsl
    #hsl = scalefactor_hsl.a*hsl0_transmural
    hsl_old = Function(Quad)

    delta_hsl = hsl - hsl_old

    cb_force = Constant(0.0)
    y_vec_split = split(y_vec)

    for jj in range(no_of_states):

        f_holder = Constant(0.0)

        if state_attached[jj] == 1:

            cb_ext = cb_extensions[jj]

            for k in range(no_of_x_bins):
                temp_holder = Constant(0.0)

                dxx = xx[k] + delta_hsl * filament_compliance_factor

                n_pop = y_vec_split[n_vector_indices[jj][0] + k]

                temp_holder = n_pop * k_cb_multiplier[jj] * (dxx + cb_ext) * conditional(gt(dxx + cb_ext,0.0), k_cb_pos, k_cb_neg)
                #temp_holder = temp_holder * conditional(gt(abs(dxx),x_bin_max),0.0,1.0)
                f_holder = f_holder + temp_holder
                #f_holder = f_holder + conditional(gt(temp_holder,0.0),temp_holder,0.0)

            f_holder = f_holder * cb_number_density * 1e-9

            f_holder = f_holder * alpha_value

        cb_force = cb_force + f_holder
    cb_force = cb_force * conditional(gt(cb_force,0.0),1.0,0.0)

    #Scalefactor = Constant(0.5)
    cb_force2 = Expression(("f"), f=0, degree=1)
    #Pactive = Scalefactor * cb_force2 * as_tensor(f0[i]*f0[j], (i,j)) #+ 0.25*cb_force * as_tensor(s0[i]*s0[j], (i,j))+ 0.25*cb_force * as_tensor(n0[i]*n0[j], (i,j))
    #Pactive, cb_force = uflforms.TempActiveStress(af_time.af_time)
    Scalefactor2 = Constant(1)
    Pactive = Scalefactor2 * cb_force * as_tensor(f0[i]*f0[j], (i,j)) + xfiber_fraction*cb_force * as_tensor(s0[i]*s0[j], (i,j))+ xfiber_fraction*cb_force * as_tensor(n0[i]*n0[j], (i,j))
    # Automatic differentiation  #####################################################################################################
    F1 = derivative(Wp, w, wtest)*dx
    F2 = inner(Fmat*Pactive, grad(v))*dx
    if(ispressurectrl):
        pressure = Expression(("p"), p=0.0, degree=2)
        F3 = inner(pressure*n, v)*ds(LVendoid)
    else:
        Wvol = uflforms.LVV0constrainedE()
        F3 = derivative(Wvol, w, wtest)
    L4 = inner(as_vector([c11[0], c11[1], 0.0]), u)*dx + \
    	 inner(as_vector([0.0, 0.0, c11[2]]), cross(X, u))*dx + \
    	 inner(as_vector([c11[3], 0.0, 0.0]), cross(X, u))*dx + \
    	 inner(as_vector([0.0, c11[4], 0.0]), cross(X, u))*dx
    F4 = derivative(L4, w, wtest)
    F5 = -Kspring*inner(dot(u,n),dot(v,n))*ds(epiid)  # traction applied as Cauchy stress!, Pactive is 1PK

    # Define circumferential direction
    zaxis = Expression(("0", "0", "1"), degree=2)
    #caxis = cross(s0CG, zaxis)
    #File(output_path + "caxis.pvd") << project(caxis, VectorFunctionSpace(mesh, "CG", 1))

    # Penalty method for enforcing no displacement in circumferential direction
    Kpen = Constant(1000000)
    #F6 = Kpen*inner(dot(u,caxis)*caxis, v)*ds(topid)


    Ftotal = F1 + F2 + F3 + F4 #+ F5#F6 + F4 #+ F5 + F4

    Jac1 = derivative(F1, w, dw)
    Jac2 = derivative(F2, w, dw)
    Jac3 = derivative(F3, w, dw)
    Jac4 = derivative(F4, w, dw)
    Jac5 = derivative(F5, w, dw)
    #Jac6 = derivative(F6, w, dw)

    Jac = Jac1 + Jac2 + Jac3 + Jac4 #+ Jac5 #Jac6 + Jac4#+ Jac5 + Jac4
    ##################################################################################################################################

    solverparams = {"Jacobian": Jac,
                    "F": Ftotal,
                    "w": w,
                    "boundary_conditions": bcs,
                    "Type": 0,
                    "mesh": mesh,
                    "mode": 0
                    }

    # Define solver, for now using default solver. NSolver is from LCLee
    solver= NSolver(solverparams)

    LVCavityvol.vol = uflforms.LVcavityvol()

    print("cavity-vol = ", LVCavityvol.vol)

    # Define paraview files to visualize on mesh
    displacementfile = File(output_path + "u_disp.pvd")
    stressfile = File(output_path + "Stress.pvd")
    pk1file = File(output_path + "pk1_act_on_f0.pvd")
    hsl_file = File(output_path + "hsl_mesh.pvd")
    alpha_file = File(output_path + "alpha_mesh.pvd")
    # Saving pressure/volume data
    if(MPI.rank(comm) == 0):
        fdataPV = open(output_path + "PV_.txt", "w", 0)
        fdataCa = open(output_path + "calcium_.txt", "w", 0)

    tstep = 0
    t = 0

    LVcav_array = [uflforms.LVcavityvol()]
    Pcav_array = [uflforms.LVcavitypressure()*0.0075]

    # Contraction phase
    tarray = []

    # Initialize arrays for saving myosim information
    hslarray = np.zeros((no_of_time_steps+1,no_of_int_points))
    calarray = np.zeros((no_of_time_steps+1,no_of_int_points))
    strarray = np.zeros((no_of_time_steps+1,no_of_int_points))
    pstrarray = np.zeros((no_of_time_steps+1,no_of_int_points))
    gucc_fiber = np.zeros((no_of_time_steps+1,no_of_int_points))
    gucc_trans = np.zeros((no_of_time_steps+1,no_of_int_points))
    gucc_shear = np.zeros((no_of_time_steps+1,no_of_int_points))
    alphaarray = np.zeros((no_of_time_steps+1,no_of_int_points))
    overlaparray = np.zeros((no_of_time_steps+1,no_of_int_points))
    deltahslarray = np.zeros((no_of_time_steps+1,no_of_int_points))
    calcium = np.zeros(no_of_time_steps+1)

    # Get array of cross-bridge populations
    y_vec_array = y_vec.vector().get_local()[:]

    hsl_array = project(sqrt(dot(f0, Cmat*f0))*hsl0, Quad).vector().get_local()[:]

    delta_hsl_array = np.zeros(no_of_int_points)

    for counter in range(0,n_array_length * no_of_int_points,n_array_length):
        # Initializing myosin heads in the off state
        y_vec_array[counter] = 1
        # Initialize all binding sites to off state
        y_vec_array[counter-2] = 1

    Pg, Pff, alpha = uflforms.stress()
    # Pg is guccione stress tensor as first Piola-Kirchhoff

    # Magnitude of bulk passive stress in fiber direction
    Pg_fiber = inner(f0,Pg*f0)
    Pg_transverse = inner(n0,Pg*n0)
    Pg_shear = inner(n0,Pg*f0)

    temp_DG = project(Pff, FunctionSpace(mesh, "DG", 1), form_compiler_parameters={"representation":"uflacs"})
    p_f = interpolate(temp_DG, Quad)
    p_f_array = p_f.vector().get_local()[:]

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

    cb_f_array = project(cb_force, Quad).vector().get_local()[:]

    # Loading phase
    print "memory growth before loading:"
    obg.show_growth()

    print("cavity-vol = ", LVCavityvol.vol)
    for lmbda_value in range(0, loading_number):

        scalefactor_epi.a = lmbda_value*0.002
        scalefactor_endo.a = lmbda_value*0.00345238

        print "Loading phase step = ", lmbda_value

        LVCavityvol.vol += 0.004 #LCL change to smaller value

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
            displacementfile << w.sub(0)
            stresstemp = project(Pactive,TensorFunctionSpace(mesh,'DG',0))
            #stresstemp = project(cb_force,FunctionSpace(mesh,'DG',0))
            stresstemp.rename("Pactive","Pactive")
            stressfile << stresstemp
            pk1temp = project(inner(f0,Pactive*f0),FunctionSpace(mesh,'DG',1))
            pk1temp.rename("pk1temp","pk1temp")
            pk1file << pk1temp
            hsl_temp = project(hsl,FunctionSpace(mesh,'DG',1))
            hsl_temp.rename("hsl_temp","hsl")
            hsl_file << hsl_temp
            alpha_temp = project(alphas,FunctionSpace(mesh,'DG',0))
            alpha_temp.rename("alpha_temp","alpha_temp")
            alpha_file << alpha_temp

        print("cavity-vol = ", LVCavityvol.vol)
        print("p_cav = ", uflforms.LVcavitypressure())


    # Closed-loop phase

    # Initialize the half-sarcomere class. Its methods will be used to solve for cell populations
    hs = half_sarcomere.half_sarcomere(hs_params,1)

    # Initialize cell ion module
    cell_ion = cell_ion_driver.cell_ion_driver(cell_ion_params)

    # Initialize calcium
    calcium[0] = cell_ion.model_class.calculate_concentrations(0,0,fdataCa)

    dumped_populations = np.zeros((no_of_time_steps+1, no_of_int_points, n_array_length))



    counter = 0
    cell_counter = 0
    cycle = 0
    AV_old = 0
    MV_old = 1
    systole = 0

    #print "memory growth before closed loop"
    #obg.show_growth()
    while(cycle < cycles):

        p_cav = uflforms.LVcavitypressure()
        V_cav = uflforms.LVcavityvol()

        tstep = tstep + step_size
        cycle = math.floor(tstep/BCL)
        cell_time = tstep - cycle*BCL


        if(MPI.rank(comm) == 0):

            print "Cycle number = ", cycle, " cell time = ", cell_time, " tstep = ", tstep, " step_size = ", step_size
            #print >>fdataPV, tstep, p_cav*0.0075 , V_cav, Myosim.Get_Ca()


        Part = 1.0/Cao*(V_art - Vart0);
        Pven = 1.0/Cven*(V_ven - Vven0);
        PLV = p_cav;

        if(MPI.rank(comm) == 0):
            print "P_ven = ",Pven;
            print "P_LV = ", PLV;
            print "P_art = ", Part;

        if(PLV <= Part):

            Qao = 0.0;
            AV_new = 0

        else:

            Qao = 1.0/Rao*(PLV - Part);
            AV_new = 1


        if(PLV >= Pven):

            Qmv = 0.0;
            MV_new = 0

        else:

            Qmv = 1.0/Rven*(Pven - PLV);
            MV_new = 1

        Qper = 1.0/Rper*(Part - Pven);

        if(MV_old == 1 and MV_new == 0):
            systole = 1
        if(AV_old == 1 and AV_new == 0):
            systole = 0

        MV_old = MV_new
        AV_old = AV_new

        if(MPI.rank(comm) == 0):

                print "Q_mv = ", Qmv ;
                print "Q_ao = ", Qao ;
                print "Q_per = ", Qper ;
                if(systole == 1):
                    print "********systole**********"
                else:
                    print "***diastole***"

        V_cav_prev = V_cav
        V_art_prev = V_art
        V_ven_prev = V_ven
        p_cav_prev = p_cav

        V_cav = V_cav + step_size*(Qmv - Qao);
        V_art = V_art + step_size*(Qao - Qper);
        V_ven = V_ven + step_size*(Qper - Qmv);

        LVCavityvol.vol = V_cav

        if(MPI.rank(comm) == 0):

                print "V_ven = ", V_ven;
                print "V_LV = ", V_cav;
                print "V_art = ", V_art;


        LVcav_array.append(V_cav)
        Pcav_array.append(p_cav*0.0075)

        if (counter > 0 and (int(counter/no_of_cell_time_steps) == (counter/no_of_cell_time_steps))):
            cell_counter = 0

        cell_counter += 1

        print "cell_counter = ", cell_counter
        """for  i in range(no_of_int_points):

            for j in range(n_array_length):

                dumped_populations[counter, i, j] = y_vec_array[i * n_array_length + j]"""

        # Initialize MyoSim solution holder
        y_vec_array_new = np.zeros(no_of_int_points*n_array_length)

        # Update calcium
        calcium[counter] = cell_ion.model_class.calculate_concentrations(cycle,cell_time,fdataCa) #LCL Commented off

        # Now print out volumes, pressures, calcium
        if(MPI.rank(comm) == 0):
            print >>fdataPV, tstep, p_cav*0.0075 , Part*.0075, Pven*.0075, V_cav, V_ven, V_art, calcium[counter]

        # Quick hack
        if counter == 0:
            overlap_counter = 1
        else:
            overlap_counter = counter
    # Going to try to loop through integration points in python, not in fenics script
        temp_overlap, y_interp, y_vec_array_new = implement.update_simulation(hs, step_size, delta_hsl_array, hsl_array, y_vec_array, p_f_array, cb_f_array, calcium[counter], n_array_length, cell_time, overlaparray[overlap_counter,:])

        for  i in range(no_of_int_points):

            for j in range(n_array_length):

                dumped_populations[counter, i, j] = y_interp[i * n_array_length + j]

        y_vec_array = y_vec_array_new # for Myosim

        #Kurtis moved to here
        y_vec.vector()[:] = y_vec_array # for PDE

        hsl_array_old = hsl_array

        # Kurtis assigning hsl_old function for newton iteration
        hsl_old.vector()[:] = hsl_array_old

        # Hack update cb_force
	"""t_trans = 30
	t0 = 20
	tr = 15
	if(cell_time < t_trans):
        	#cb_force2.f = 70000*sin(cell_time/40.0*3.14)
        	cb_force2.f = 120000*0.5*(1 - cos(pi*cell_time/t0))
		#print cb_force2.f
	else:
		A =  120000*0.5*(1 - cos(pi*t_trans/t0))
		cb_force2.f = A*exp(-1.0*(cell_time - t_trans)/tr)"""

###########################################################################

        #solver.solvenonlinear()
        #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        try:
            solve(Ftotal == 0, w, bcs, J = Jac, form_compiler_parameters={"representation":"uflacs"})
        except:
            print "Newton Iteration non-convergence, saving myosim info"
            np.save(output_path +"dumped_populations", dumped_populations)
            np.save(output_path + "tarray", tarray)
            np.save(output_path + "stress_array", strarray)
            np.save(output_path + "hsl", hslarray)
            np.save(output_path + "overlap", overlaparray)
            np.save(output_path + "gucc_fiber", gucc_fiber)
            np.save(output_path + "gucc_trans", gucc_trans)
            np.save(output_path + "gucc_shear", gucc_shear)
            np.save(output_path + "deltahsl", deltahslarray)
            np.save(output_path + "pstress_array",pstrarray)
            #np.save(output_path + "alpha_array",alphaarray)
            np.save(output_path + "calcium",calarray)
        #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        cb_f_array = project(cb_force, Quad).vector().get_local()[:]

        hsl_array = project(hsl, Quad).vector().get_local()[:]           # for Myosim

        delta_hsl_array = hsl_array - hsl_array_old

        temp_DG = project(Pff, FunctionSpace(mesh, "DG", 1), form_compiler_parameters={"representation":"uflacs"})
        p_f = interpolate(temp_DG, Quad)
        p_f_array = p_f.vector().get_local()[:]

        for ii in range(np.shape(hsl_array)[0]):
            if p_f_array[ii] < 0.0:
                p_f_array[ii] = 0.0

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

        counter += 1
        displacementfile << w.sub(0)
        stresstemp = project(Pactive,TensorFunctionSpace(mesh,'DG',0))
        stresstemp.rename("Pactive","Pactive")
        stressfile << stresstemp
        pk1temp = project(inner(f0,Pactive*f0),FunctionSpace(mesh,'DG',1))
        pk1temp.rename("pk1temp","pk1temp")
        pk1file << pk1temp
        hsl_temp = project(hsl,FunctionSpace(mesh,'DG',1))
        hsl_temp.rename("hsl_temp","hsl")
        hsl_file << hsl_temp
        alpha_temp = project(alphas,FunctionSpace(mesh,'DG',0))
        alpha_temp.rename("alpha_temp","alpha_temp")
        alpha_file << alpha_temp

        tarray.append(tstep)


        hslarray[counter,:] = hsl_array
        strarray[counter,:] = cb_f_array
        pstrarray[counter,:] = p_f_array
        gucc_fiber[counter,:] = pgf_array
        gucc_trans[counter,:] = pgt_array
        gucc_shear[counter,:] = pgs_array
        alphaarray[counter,:] = alpha_array
        overlaparray[counter,:] = temp_overlap
        deltahslarray[counter,:] = delta_hsl_array
        calarray[counter,:] = hs.Ca_conc * np.ones(no_of_int_points)


    if(MPI.rank(comm) == 0):
        fdataPV.close()
        fdataCa.close()

    fluxes, rates = implement.return_rates_fenics(hs)


    # Generate dictionary for output
    outputs = {
    "rates": rates,
    "dumped_populations": dumped_populations,
    "tarray": tarray,
    "strarray": strarray,
    "pstrarray": pstrarray,
    "gucc_fiber": gucc_fiber,
    "gucc_trans": gucc_trans,
    "gucc_shear": gucc_shear,
    "alphaarray": alphaarray,
    "calarray": calarray,
    "hsl": hslarray,
    "overlap": overlaparray
    }


    return(outputs)

    ######################################################################################################
