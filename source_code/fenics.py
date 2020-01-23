from __future__ import division
import sys
sys.path.append("/home/fenics/shared/source_code/dependencies")
import os as os
from dolfin import *
import numpy as np
from matplotlib import pylab as plt
from petsc4py import PETSc
from forms import Forms
from nsolver import NSolver as NSolver
import math
import json
import Python_MyoSim.half_sarcomere.half_sarcomere as half_sarcomere
import Python_MyoSim.half_sarcomere.implement as implement
import vtk_py

## Import parameters from json file
# User provides file name
input_file_name = sys.argv[1]

## Read in the JSON tree structure.
# Check that the file exists, if not , exit
if not os.path.exists(input_file_name):
    print "input file does not exist"
    exit()

with open(input_file_name, 'r') as json_input:
  input_parameters = json.load(json_input)

## Parse out the different types of parameters
sim_params = input_parameters["simulation_parameters"]
file_inputs = input_parameters["file_inputs"]
output_params = input_parameters["output_parameters"]
passive_params = input_parameters["forms_parameters"]["passive_law_parameters"]
hs_params = input_parameters["myosim_parameters"]
cell_ion_params = input_parameters["electrophys_parameters"]["cell_ion_parameters"]
monodomain_params = input_parameters["electrophys_parameters"]["monodomain_parameters"]

## Assign parameters
output_path = output_params["output_path"][0]
input_path = file_inputs["input_directory_path"][0]
#casename = file_inputs["casename"][0]
casename = "ellipsoidal"
casename.type()
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
calcium_path = cell_ion_params["path_to_calcium"][0]
# Check that the output path exists. If it does not, create it and let user know
if not os.path.exists(output_path):
    print "Output path does not exist. Creating it now"
    os.makedirs(output_path)

# Initialize the half-sarcomere class. Its methods will be used to solve for cell populations
hs = half_sarcomere.half_sarcomere(hs_params,1)

# The cell active stress calculation is in the fenics file due to the
# implicit nature of the problem. Active stress and length depend on one another

# As of 1/17/2020, the "3state_with_SRX" scheme is kind of hard coded in Python_MyoSim.
# Want this force calculation in fenics to be for any arbitrary scheme. The input readme
# contains information for each of these parameters


# Set up x vector
xx = np.arange(x_bin_min, x_bin_max + x_bin_increment, x_bin_increment)

# Define number of intervals cross-bridges are defined over
no_of_x_bins = np.shape(xx)[0]

# Define the length of the populations vector
n_array_length = no_of_attached_states * no_of_x_bins + no_of_detached_states + 2
# +2 for binding sites on/off

# Work out a general way to set this based on the scheme
n_vector_indices = [[0,0], [1,1], [2,2+no_of_x_bins-1]]

# Cardiac cycle length and number of cycles will be general
BCL = 100 # ms
cycles = 1

hsl0 = 1000
step_size = 0.5
no_of_time_steps = int(cycles*BCL/step_size)
no_of_cell_time_steps = int(BCL/step_size)
Ca_flag = 4
constant_pCa = 6.5

# Loading calcium from previous simulation
prev_ca = np.load(calcium_path)
prev_ca = prev_ca[:,0]

deg = 4
parameters["form_compiler"]["quadrature_degree"]=deg
parameters["form_compiler"]["representation"] = "quadrature"

os.system("rm ../python_dev/results/fenics_test/*.pvd")
os.system("rm ../python_dev/results/fenics_test/*.vtu")

############################## Insert Mesh ###########################################
#casename = "ellipsoidal"
meshfilename = casename + ".hdf5"
print meshfilename
mesh = Mesh()

f = HDF5File(mpi_comm_world(), meshfilename, 'r')
f.read(mesh, casename, False)

if casename == "ellipsoidal":
    loading_number = 2;
    #don't need to do the vtk_py mesh stuff
else: #assuming we are using a patient specific mesh
    loading_number = 1;
    ugrid = vtk_py.convertXMLMeshToUGrid(mesh)
    ugrid = vtk_py.rotateUGrid(ugrid, sx=0.1, sy=0.1, sz=0.1)
    mesh = vtk_py.convertUGridToXMLMesh(ugrid)

no_of_int_points = 14 * np.shape(mesh.cells())[0]

facetboundaries = MeshFunction("size_t", mesh, 2)
VQuadelem = VectorElement("Quadrature", mesh.ufl_cell(), degree=deg, quad_scheme="default")
VQuadelem._quad_scheme = 'default'
fiberFS = FunctionSpace(mesh, VQuadelem)

f0 = Function(fiberFS)
s0 = Function(fiberFS)
n0 = Function(fiberFS)
f.read(facetboundaries, casename+"/"+"facetboundaries")

f.read(f0, casename+"/"+"eF")
f.read(s0, casename+"/"+"eS")
f.read(n0, casename+"/"+"eN")
f.close()
File("facetboundaries.pvd") << facetboundaries
topid = 4
LVendoid = 2
epiid = 1
##############################################################################

comm = mesh.mpi_comm()

File("fiber.pvd") << project(f0, VectorFunctionSpace(mesh, "CG", 1))
File("sheet.pvd") << project(s0, VectorFunctionSpace(mesh, "CG", 1))
File("sheet-normal.pvd") << project(n0, VectorFunctionSpace(mesh, "CG", 1))

##############################################################################


isincomp = True#False
N = FacetNormal (mesh)
Press = Expression(("P"), P=0.0, degree=0)
Kspring = Constant(100)
LVCavityvol = Expression(("vol"), vol=0.0, degree=2)
#C2 = Constant(0.5)
#C3 = Constant(25.0)
#Cparam = Constant(5.0e2)


V = VectorFunctionSpace(mesh, 'CG', 2)
TF = TensorFunctionSpace(mesh, 'DG', 1)
Q = FunctionSpace(mesh,'CG',1)

Velem = VectorElement("CG", mesh.ufl_cell(), 2, quad_scheme="default")
Velem._quad_scheme = 'default'
Qelem = FiniteElement("CG", mesh.ufl_cell(), 1, quad_scheme="default")
Qelem._quad_scheme = 'default'
Relem = FiniteElement("Real", mesh.ufl_cell(), 0, quad_scheme="default")
Relem._quad_scheme = 'default'
#Quadelem = FiniteElement("Quadrature", mesh.ufl_cell(), degree=deg, quad_scheme="default")
Quadelem = FiniteElement("Quadrature", tetrahedron, degree=deg, quad_scheme="default")
Quadelem._quad_scheme = 'default'

Telem2 = TensorElement("Quadrature", mesh.ufl_cell(), degree=deg, shape=2*(3,), quad_scheme='default')
Telem2._quad_scheme = 'default'
for e in Telem2.sub_elements():
    e._quad_scheme = 'default'
Telem4 = TensorElement("Quadrature", mesh.ufl_cell(), degree=deg, shape=4*(3,), quad_scheme='default')
Telem4._quad_scheme = 'default'
for e in Telem4.sub_elements():
    e._quad_scheme = 'default'
W = FunctionSpace(mesh, MixedElement([Velem,Qelem,Relem]))
Quad = FunctionSpace(mesh, Quadelem)

Quad_vectorized_Fspace = FunctionSpace(mesh, MixedElement(n_array_length*[Quadelem]))

bctop = DirichletBC(W.sub(0).sub(2), Expression(("0.0"), degree = 2), facetboundaries, topid)
bcs = [bctop]

w = Function(W)
dw = TrialFunction(W)
wtest = TestFunction(W)
du,dp,dpendo = TrialFunctions(W)
(u,p, pendo) = split(w)
(v,q, qendo) = TestFunctions(W)

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
#         "C2": C2,
#         "C3": C3,
#         "Cparam": Cparam,
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

Ematrix = project(Emat, TF)
Wp = uflforms.PassiveMatSEF()
Wvol = uflforms.LVV0constrainedE()

#Active force calculation------------------------------------------------------
y_vec = Function(Quad_vectorized_Fspace)
hsl = sqrt(dot(f0, Cmat*f0))*hsl0
hsl_old = Function(Quad)
delta_hsl = hsl - hsl_old

#f_holder = Constant(0.0)
cb_force = Constant(0.0)

y_vec_split = split(y_vec)

for jj in range(no_of_states):

    f_holder = Constant(0.0)

    if state_attached[jj] == 1:

        cb_ext = cb_extensions[jj]

        for k in range(no_of_x_bins):

            dxx = xx[k] + delta_hsl * filament_compliance_factor

            n_pop = y_vec_split[n_vector_indices[jj][0] + k]

            f_holder = f_holder + n_pop * k_cb_multiplier[jj] * (dxx + cb_ext) * conditional(dxx + cb_ext > 0.0, k_cb_pos, k_cb_neg)

        f_holder = f_holder * cb_number_density * 1e-9

        f_holder = f_holder * alpha_value

    cb_force = cb_force + f_holder

Pactive = cb_force * as_tensor(f0[i]*f0[j], (i,j))

# Automatic differentiation  #####################################################################################################
F1 = derivative(Wp, w, wtest)*dx
F2 = inner(Pactive, grad(v))*dx
F3 = derivative(Wvol, w, wtest)
F4 = -Kspring*inner(dot(u,n)*n,v)*ds(epiid)  # traction applied as Cauchy stress!, Pactive is 1PK
Ftotal = F1 + F2 + F3 + F4

Jac1 = derivative(F1, w, dw)
Jac2 = derivative(F2, w, dw)
Jac3 = derivative(F3, w, dw)
Jac4 = derivative(F4, w, dw)
Jac = Jac1 + Jac2 + Jac3 + Jac4
##################################################################################################################################

solverparams = {"Jacobian": Jac,
                "F": Ftotal,
                "w": w,
                "boundary_conditions": bcs,
                "Type": 0,
                "mesh": mesh,
                "mode": 0
                }


solver= NSolver(solverparams)


LVCavityvol.vol = uflforms.LVcavityvol()

print("cavity-vol = ", LVCavityvol.vol)

#displacementfile = File("./test_14/u_disp.pvd")

if(MPI.rank(comm) == 0):
    fdataPV = open("/home/fenics/shared/python_dev/results/pm_test14/PV_.txt", "w", 0)


tstep = 0
t = 0

LVcav_array = [uflforms.LVcavityvol()]
Pcav_array = [uflforms.LVcavitypressure()*0.0075]

# Contraction phase

tarray = []

hslarray = np.zeros((no_of_time_steps+1,no_of_int_points))
calarray = np.zeros((no_of_time_steps+1,no_of_int_points))
strarray = np.zeros((no_of_time_steps+1,no_of_int_points))
pstrarray = np.zeros((no_of_time_steps+1,no_of_int_points))
alphaarray = np.zeros((no_of_time_steps+1,no_of_int_points))

y_vec_array = y_vec.vector().get_local()[:]

#hsl_array = project(sqrt(dot(f0, Cmat*f0))*hsl0, Quad).vector().get_local()[:]
hsl_array = project(hsl, Quad).vector().get_local()[:]

delta_hsl_array = np.zeros(no_of_int_points)

for counter in range(0,n_array_length * no_of_int_points,n_array_length):
    # Initializing myosin heads in the off state
    y_vec_array[counter] = 1
    # Initialize all binding sites to off state
    y_vec_array[counter-2] = 1

Pff, alpha = uflforms.stress()

temp_DG = project(Pff, FunctionSpace(mesh, "DG", 1), form_compiler_parameters={"representation":"uflacs"})
p_f = interpolate(temp_DG, Quad)
p_f_array = p_f.vector().get_local()[:]

temp_DG_1 = project(alpha, FunctionSpace(mesh, "DG", 1), form_compiler_parameters={"representation":"uflacs"})
alphas = interpolate(temp_DG_1, Quad)
alpha_array = alphas.vector().get_local()[:]

cb_f_array = project(cb_force, Quad).vector().get_local()[:]

temp_DG_1 = project(alpha, FunctionSpace(mesh, "DG", 1), form_compiler_parameters={"representation":"uflacs"})
alphas = interpolate(temp_DG_1, Quad)


# Loading phase
print("cavity-vol = ", LVCavityvol.vol)
for lmbda_value in range(0, loading_number):

    print "Loading phase step = ", lmbda_value
    LVCavityvol.vol += 0.005

    p_cav = uflforms.LVcavitypressure()
    V_cav = uflforms.LVcavityvol()

    hsl_array_old = hsl_array

    #solver.solvenonlinear()
    solve(Ftotal == 0, w, bcs, J = Jac, form_compiler_parameters={"representation":"uflacs"})
    hsl_old.vector()[:] = project(hsl, Quad).vector().get_local()[:] # for active stress

    hsl_array = project(hsl, Quad).vector().get_local()[:]           # for Myosim
    #print(np.shape(hsl_array))
    delta_hsl_array = project(sqrt(dot(f0, Cmat*f0))*hsl0, Quad).vector().get_local()[:] - hsl_array_old # for Myosim

    temp_DG = project(Pff, FunctionSpace(mesh, "DG", 1), form_compiler_parameters={"representation":"uflacs"})
    p_f = interpolate(temp_DG, Quad)
    p_f_array = p_f.vector().get_local()[:]

    temp_DG_1 = project(alpha, FunctionSpace(mesh, "DG", 1), form_compiler_parameters={"representation":"uflacs"})
    alphas = interpolate(temp_DG_1, Quad)
    alpha_array = alphas.vector().get_local()[:]

    for ii in range(np.shape(alpha_array)[0]):
        #if (alpha_array[ii] < 1.0):
         #   p_f_array[ii] = 0.0
         if(p_f_array[ii] < 0.0):
             p_f_array[ii] = 0.0

    if(MPI.rank(comm) == 0):

        print >>fdataPV, 0.0, p_cav*0.0075 , V_cav, 0.0
        #displacementfile << w.sub(0)
    print("cavity-vol = ", LVCavityvol.vol)
    print("p_cav = ", uflforms.LVcavitypressure())


# Closed-loop phase
dumped_populations = np.zeros((no_of_time_steps+1, no_of_int_points, n_array_length))


Cao = 0.02/1000.0; #5.0e-5; #Cao = 0.005;
Cven = 2.0/1000.0 #0.02;#Cven = 0.2*10;
Vart0 = 0.0#100.0/1000.0#/2000;
Vven0 = 0.0#1000.0/1000.0#0.001 * 2200.0/1000.0#/2000;
Rao = 10*1000.0*1000.0#Rao = 10*1000.0;
Rven = 1000*1000.0#Rven = 1000.0;
Rper = 200000 * 50 #1000.0;
V_ven = 1200.0/1000.0#/2000;
V_art = 250.0/1000.0#/2000;


counter = 0
cell_counter = 0
cycle = 0
AV_old = 0
MV_old = 1
systole = 0
while(cycle < cycles):

    p_cav = uflforms.LVcavitypressure()
    V_cav = uflforms.LVcavityvol()

    tstep = tstep + step_size
    cycle = math.floor(tstep/BCL)
    cell_time = tstep - cycle*BCL


    if(MPI.rank(comm) == 0):

        print "Cycle number = ", cycle, " cell time = ", cell_time, " tstep = ", tstep, " step_size = ", step_size
        #print >>fdataPV, tstep, p_cav*0.0075 , V_cav, Myosim.Get_Ca()
        print >>fdataPV, tstep, p_cav*0.0075 , V_cav

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
    for  i in range(no_of_int_points):

        for j in range(n_array_length):

            dumped_populations[counter, i, j] = y_vec_array[i * n_array_length + j]


    y_vec.vector()[:] = y_vec_array # for PDE

    # Initialize MyoSim solution holder
    y_vec_array_new = np.zeros(no_of_int_points*n_array_length)
   # Checking size of things passed if __name__ == '__main__':
    #print np.shape(hsl_array)
    #print np.shape(delta_hsl_array)
    #print np.shape(p_f_array)
    #print np.shape(y_vec_array_new)

# Going to try to loop through integration points in python, not in fenics script
    y_vec_array_new = implement.update_simulation(hs, step_size, delta_hsl_array, hsl_array, y_vec_array, p_f_array, cb_f_array, prev_ca[counter], n_array_length)

#    y_vec_array = y_vec_array_new
#    print no_of_int_points
#    for k in range(no_of_int_points):
#        pop_holder = implement.update_simulation(hs, step_size, a_delta_hsl[counter], amir_hsl[counter], y_vec_array[k*n_array_length:(k+1)*n_array_length], p_f_array[k], cb_f_array[k], prev_ca[counter])

#        pop_holder = implement.update_simulation(hs, step_size, delta_hsl_array[k], hsl_array[k], y_vec_array[k*n_array_length:(k+1)*n_array_length], p_f_array[k], cb_f_array[k], prev_ca[counter])



#        y_vec_array_new[k*n_array_length:(k+1)*n_array_length] = pop_holder

    y_vec_array = y_vec_array_new # for Myosim

    hsl_array_old = hsl_array

    #solver.solvenonlinear()
    #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    solve(Ftotal == 0, w, bcs, J = Jac, form_compiler_parameters={"representation":"uflacs"})
    #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    cb_f_array = project(cb_force, Quad).vector().get_local()[:]

    hsl_old.vector()[:] = project(hsl, Quad).vector().get_local()[:] # for active stress

    hsl_array = project(hsl, Quad).vector().get_local()[:]           # for Myosim
    #print(np.shape(hsl_array))
    delta_hsl_array = project(sqrt(dot(f0, Cmat*f0))*hsl0, Quad).vector().get_local()[:] - hsl_array_old # for Myosim

    temp_DG = project(Pff, FunctionSpace(mesh, "DG", 1), form_compiler_parameters={"representation":"uflacs"})
    p_f = interpolate(temp_DG, Quad)
    p_f_array = p_f.vector().get_local()[:]

    temp_DG_1 = project(alpha, FunctionSpace(mesh, "DG", 1), form_compiler_parameters={"representation":"uflacs"})
    alphas = interpolate(temp_DG_1, Quad)
    alpha_array = alphas.vector().get_local()[:]

    for ii in range(np.shape(alpha_array)[0]):
        #if (alpha_array[ii] < 1.0):
         #   p_f_array[ii] = 0.0
         if(p_f_array[ii] < 0.0):
             p_f_array[ii] = 0.0

    calarray[counter,:] = hs.Ca_conc * np.ones(no_of_int_points)

    counter += 1
    #displacementfile << w.sub(0)

    tarray.append(tstep)


    hslarray[counter,:] = hsl_array
    strarray[counter,:] = cb_f_array
    pstrarray[counter,:] = p_f_array
    alphaarray[counter,:] = alpha_array

if(MPI.rank(comm) == 0):
    fdataPV.close()

#rate_constants = np.zeros((no_of_x_bins,no_of_transitions + 1))

fluxes, rates = implement.return_rates_fenics(hs)

#np.save("/home/fenics/shared/test_python_myosim/rates",rate_constants)
np.save(output_path + "rates",rates)
#
np.save(output_path + "dumped_populations",dumped_populations)
#
np.save(output_path + "tarray",tarray)
#
np.save(output_path + "stress_array",strarray)
#
np.save(output_path + "pstress_array",pstrarray)
#
np.save(output_path + "alpha_array",alphaarray)
#
np.save(output_path + "calcium",calarray)
#
np.save(output_path + "HSL",hslarray)

######################################################################################################
