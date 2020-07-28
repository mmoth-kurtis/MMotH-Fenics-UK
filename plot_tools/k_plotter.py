"""
Created on Wed Oct  9 12:40:00 PM 2019

@author: Kurtis Mann

Modifying Amir's plotter to hopefully plot generic LV simulation
Currently, IPython only works using python 2.7
"""
import os as os
import numpy as np
import matplotlib.pyplot as plt
import sys
from matplotlib.animation import FuncAnimation
import time
plt.style.use('seaborn-pastel')

#(Commented out by Kurtis 10/9/2019, not sure what this does or if needed)
"""
no_of_time_steps = 701
n_array_length = 52

no_of_time_steps = 7408
n_array_length = 45
Commented out by Kurtis 10/9/2019, hope to get this info from loaded arrays"""

# base_dir should be shared directory with FEniCS container
# sim_dir should be location of simulation to plot relative to base_dir
#base_dir = '/Users/charlesmann/Academic/UK/MMotH-Fenics-UK/working_directory_untracked/'
#sim_dir = '082315/mouse_wk_params/fixed_fibers/angle_60/more_refined/timestep_small/k3_1000/'
sim_dir = os.getcwd()
lang_flag = 'python'
#sim_dir = 'cpp_stable/test_17/'
#lang_flag = 'c++'

# For now, hard coding bin discretization information
xmin = -10
xmax = 10
bin_width = 1.0
cb_domain = np.arange(xmin,xmax+bin_width,bin_width)
num_bins = np.shape(cb_domain)


# Assuming the dumped npy files from FEniCS always take these names
# Note, these arrays include info for every Gauss point
fenics_pop_file = np.load(sim_dir + '/dumped_populations.npy')
tarray = np.load(sim_dir + '/tarray.npy')
#tarray = np.arange(0,340,0.1)
tarray = tarray[:-1]
stress_array = np.load(sim_dir + '/stress_array.npy')

#print stress_array.ndim
if stress_array.ndim > 1:
    # single cell sims only save for one gauss point, dimension is one less than
    # ventricle simulation
    single_cell_sim_flag = 0
else:
    single_cell_sim_flag = 1
#stress_array = np.load(sim_dir + '/strarray.npy')

calcium = np.load(sim_dir + '/calcium.npy')
#calcium = calcium[:,0]
HSL = np.load(sim_dir + '/hsl.npy')
#HSL = np.load(sim_dir + '/hslarray.npy')
pstress = np.load(sim_dir + '/pstress_array.npy')
if single_cell_sim_flag == 0:
    gucc_fiber = np.load(sim_dir + '/gucc_fiber.npy')
    gucc_trans = np.load(sim_dir + '/gucc_trans.npy')
    gucc_shear = np.load(sim_dir + '/gucc_shear.npy')
overlap = np.load(sim_dir + '/overlap.npy')
# Define number of time steps and array length here
sim_info = fenics_pop_file.shape
num_timesteps = sim_info[0]
num_int_points = sim_info[1]
array_length = sim_info[2]
#print sys.argv[0]
#print sys.argv[1]
gauss_point = int(sys.argv[1])

#gauss_point = 1000
data_range = np.shape(tarray)[0]
#data_range = 50
# Look at how info is dumped from FEniCS. For now, hard code number of detached and attached states, and bins
# Want to be able to visualize distributions, will need this info to set up arrays.
#num_d_states
#num_a_states
#num_bins
#bin_min
#bin_max

fenics_pop_data = np.zeros((num_timesteps,array_length))
for i in range(data_range):

    # Reading in information from just one Gauss point [i = timestep, 0 = gauss point, : is all pop info]
    fenics_pop_data[i,0] = fenics_pop_file[i,gauss_point,0] # state 1 pops
    fenics_pop_data[i,1] = fenics_pop_file[i,gauss_point,1] # state 2 pops
    fenics_pop_data[i,2] = np.sum(fenics_pop_file[i,gauss_point,2:array_length-3]) # state 3
    fenics_pop_data[i,3] = fenics_pop_file[i,0,array_length-1]
    fenics_pop_data[i,4:] = fenics_pop_file[i,0,4:array_length]

"""fenics_rates_file = sim_dir + '/rates.npy'

if lang_flag=='python':
    rates = np.load(fenics_rates_file, allow_pickle=True)
    r1 = rates.item().get('R1')
    r2 = rates.item().get('R2')
    r3 = rates.item().get('R3')
    r4 = rates.item().get('R4')
fenics_rates = np.zeros((array_length-3,5))
#fenics_rates[:,0:5] = np.loadtxt(fenics_rates_file, skiprows = 1, usecols = (0,1,2,3,4))
if lang_flag=='c++':
    fenics_rates[:,0:5] = np.load(fenics_rates_file)
#get_ipython().run_line_magic('matplotlib', 'qt')
#get_ipython().run_line_magic('matplotlib', 'inline')"""

fig = plt.figure()
#------------------------------------------------------------------------------
plt.subplot(422)
#Added -1 so it plots fenics_pop_data[:-1,#] instead of entire array. Time array is one short for some reason
state_1_pops_fenics, = plt.plot(tarray[0:data_range], fenics_pop_data[0:data_range,0],label='SRX')
state_2_pops_fenics, = plt.plot(tarray[0:data_range], fenics_pop_data[0:data_range,1],label='Detached')
state_3_pops_fenics, = plt.plot(tarray[0:data_range], fenics_pop_data[0:data_range,2],label='M Bound')

"""plt.scatter(tarray[::10], myosim_pop_data[::10,0], color = 'k')
plt.scatter(tarray[::10], myosim_pop_data[::10,1], color = 'b')
plt.scatter(tarray[::10], myosim_pop_data[::10,2], color = 'r')"""

plt.legend((state_1_pops_fenics, state_2_pops_fenics, state_3_pops_fenics), ('OFF', 'ON', 'FG'))
plt.title("Myosin Populations")
plt.xlabel('time (s)')
plt.ylabel("Proportions")


#------------------------------------------------------------------------------
plt.subplot(424)
#state_3_pops_fenics, = plt.plot(tarray, np.sum(fenics_pop_data[0:data_range,2:array_length-2]), 'r')
state_3_pops_fenics, = plt.plot(tarray[0:data_range], fenics_pop_data[0:data_range,2])
binding_sites, = plt.plot(tarray[0:data_range], fenics_pop_file[0:data_range,gauss_point,array_length-1])

#plt.scatter(tarray, myosim_pop_data[:,2], color = 'r')
#plt.scatter(tarray[::10], myosim_pop_data[::10,2], color = 'r')

plt.legend((state_3_pops_fenics, binding_sites), ('Xbridges', 'Binding sites'))
plt.xlabel('time (s)')
plt.ylabel("Proportions")
#------------------------------------------------------------------------------
ax3 = plt.subplot(426)
if single_cell_sim_flag > 0:
    plt.plot(tarray[0:data_range],stress_array[0:data_range])
else:
    plt.plot(tarray[0:data_range], stress_array[0:data_range,gauss_point])
#plt.plot(tarray, stress_array[0:data_range])
#plt.scatter(myosim_summary_data[:,0], myosim_summary_data[:,1],color='r')
#plt.scatter(myosim_summary_data[::10,0], myosim_summary_data[::10,1],color='r')
plt.xlabel('time (s)')
plt.ylabel("Stress (Pa)")

#------------------------------------------------------------------------------
plt.subplot(428)
plt.plot(tarray[0:data_range], calcium[0:data_range,0])
#plt.scatter(myosim_summary_data[:,0], myosim_summary_data[:,1],color='r')
plt.xlabel('time (s)')
plt.ylabel("Calcium [M]")
#------------------------------------------------------------------------------
ax2 = plt.subplot(421)
if single_cell_sim_flag > 0:
    plt.plot(tarray[0:data_range],HSL[0:data_range])
else:
    plt.plot(tarray[0:data_range], HSL[0:data_range,gauss_point])
#plt.plot(tarray, HSL[0:data_range])
#plt.scatter(myosim_summary_data[:,0], myosim_summary_data[:,1],color='r')
plt.xlabel('time [s]')
plt.ylabel("hsl (nm)")

#---------------------------------------------------------------------------------
plt.subplot(423)
if single_cell_sim_flag > 0:
    fiber_pstress, = plt.plot(tarray[0:data_range], pstress[0:data_range])
else:
    fiber_pstress, = plt.plot(tarray[0:data_range], pstress[0:data_range,gauss_point])
    gfiber, = plt.plot(tarray[0:data_range], gucc_fiber[0:data_range,gauss_point])
    gtrans, = plt.plot(tarray[0:data_range], gucc_trans[0:data_range,gauss_point])
    gshear, = plt.plot(tarray[0:data_range], gucc_shear[0:data_range,gauss_point])
    plt.legend((fiber_pstress, gfiber, gtrans, gshear), ('fiber', 'G_fiber', 'G_trans', 'G_shear'))

#plt.plot(tarray, pstress[0:data_range])
plt.ylabel('Passive Stress (Pa)')
#------------------------------------------------------------------------------
plt.subplot(425)
plt.plot(tarray[0:data_range], overlap[0:data_range,gauss_point])
plt.ylabel('Overlap')
#plt.subplot(423)
#plt.scatter(myosim_rates[:,0], myosim_rates[:,1],color='k')
"""if lang_flag=='python':
    rate3 = np.zeros(num_bins)
    rate4 = np.zeros(num_bins)
    # Right now, getting these fluxes at last time point
    # J1, J2 are single values, not worth plotting
    # Going to try to divide by populations to get rates instead of fluxes
    #for l in range(num_bins[0]):
#        rate3[l] = j3[l]/(fenics_pop_data[l,1]*(fenics_pop_file[l,gauss_point,array_length-1]-fenics_pop_file[-1,gauss_point,2+l])*bin_width)
#        rate4[l] = j4[l]/fenics_pop_data[l,2]
    plt.subplot(423)
    flux1, = plt.plot(0,r1,'o')
    flux2, = plt.plot(0,r2,'o')
    plt.subplot(425)
    flux3, = plt.plot(cb_domain,r3)
    flux4, = plt.plot(cb_domain,r4)
else:
    plt.subplot(423)
    rate1, = plt.plot(fenics_rates[:,0], fenics_rates[:,1])
    #plt.scatter(myosim_rates[:,0], myosim_rates[:,2],color='r')
    rate2, = plt.plot(fenics_rates[:,0], fenics_rates[:,2])
    plt.legend((rate1, rate2), ('Rate 1', 'Rate 2'))

    #------------------------------------------------------------------------------
    plt.subplot(425)
    #plt.scatter(myosim_rates[:,0], myosim_rates[:,3],color='g')
    rate3, = plt.plot(fenics_rates[:,0], fenics_rates[:,3])
    #plt.scatter(myosim_rates[:,0], myosim_rates[:,4],color='b')
    rate4, = plt.plot(fenics_rates[:,0], fenics_rates[:,4])
    plt.legend((rate3, rate4), ('Attach', 'Detach'))
"""
#------------------------------------------------------------------------------
# Animate cross-bridges during simulation
max_nbound = np.max(fenics_pop_data[:,2])
#print max_nbound
ax1 = plt.subplot(427,xlim=(xmin-1,xmax+1),ylim=(0.00,max_nbound/2))
#ax = plt.axes(xlim=(xmin,xmax),ylim=(0,1))
line1, = ax1.plot([],[],lw=3)
line2, = ax2.plot([],[])
line3, = ax3.plot([],[])
line = [line1, line2, line3]

def init():
    line[0].set_data([],[])
    line[1].set_data([],[])
    line[2].set_data([],[])
    return line

t, m, m2 = [], [], []
y = np.zeros(np.shape(cb_domain))
def animate(i):
    # array_length -1 for cpp, -2 for python
    y = fenics_pop_file[i,gauss_point,2:array_length-2]
    if single_cell_sim_flag > 0:
        m.append(HSL[i])
        m2.append(stress_array[i])
    else:
        m.append(HSL[i,gauss_point])
        m2.append(stress_array[i,gauss_point])
    #print np.shape(cb_domain)
    t.append(tarray[i])
    #print np.shape(y)
    line[0].set_data(cb_domain,y)
    line[1].set_data(t,m)
    line[2].set_data(t,m2)
    time.sleep(0.5)
    return line


anim = FuncAnimation(fig, animate, init_func=init, frames = num_timesteps-1, interval = 1, blit=True)

#mng = plt.get_current_fig_manager()
#mng.frame.Maximize(True)
#plt.figure()

plt.get_current_fig_manager().full_screen_toggle()
plt.show()
