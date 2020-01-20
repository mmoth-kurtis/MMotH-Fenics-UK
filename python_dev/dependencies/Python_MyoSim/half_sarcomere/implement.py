# Functions for implementating half-sarcomere class
import numpy as np
import pandas as pd


def update_simulation(self, time_step, delta_hsl, hsl, y0, pf, cbf, calcium, n_array_length, set_data = 0):

    # Need to do some kinetics stuff
    time_step = time_step/1000
    # Update calcium
    #self.membr.evolve_kinetics(time_step, activation)
    #self.Ca_conc = self.membr.myofilament_Ca_conc

    # Update calcium from values loaded in from previous simulation
    # Update all properties of self passed in from FE sim
    self.Ca_conc = calcium

    # Going to loop through integration points here, neglecting myosin isoforms for now
    num_int_points = np.shape(hsl)
    num_int_points = num_int_points[0]
    y_pops = np.zeros(np.size(y0))

    #print y0[0:53]
    for i in range(num_int_points):
        self.hs_length = hsl[i]
        self.myof.cb_force = cbf[i]
        self.myof.pas_force = pf[i]
        self.hs_force = self.myof.cb_force+self.myof.pas_force
        self.myof.y[0:n_array_length] = y0[i*n_array_length:(i+1)*n_array_length]
#        print "myof_y" + str(np.shape(self.myof.y))
        #if i==1:
            #print self.myof.y
        # Myofilaments
        self.myof.evolve_kinetics(time_step, self.Ca_conc)
        #if i==1:
            #print self.myof.y
        if (np.abs(delta_hsl[i]) > 0.0):
            # Need to move some things
            self.myof.move_cb_distributions(delta_hsl[i])
        #self.hs_length = self.hs_length + delta_hsl

        # Assign int point's population vector to larger y vector
        y_pops[i*n_array_length:(i+1)*n_array_length] = self.myof.y[0:n_array_length]

    # Update forces
    # Don't need these, will be reset at beginning of this fcn
    #self.myof.set_myofilament_forces()
    #self.hs_force = self.myof.total_force
    #print y_pops[0:53]
    return y_pops

def return_rates_fenics(self):
    fluxes, rates = self.myof.return_fluxes(self.myof.y, self.Ca_conc)

    return fluxes, rates

def update_data_holder(self, dt, activation):
    # Update data struct for half-sarcomere
    self.data_buffer_index = self.data_buffer_index + 1
    self.hs_time = self.hs_time + dt
    self.hs_data.at[self.data_buffer_index, 'hs_time'] = self.hs_time
    self.hs_data.at[self.data_buffer_index, 'activation'] = activation
    self.hs_data.at[self.data_buffer_index, 'Ca_conc'] = self.Ca_conc

    self.hs_data.at[self.data_buffer_index, 'hs_length'] = self.hs_length

    self.hs_data.at[self.data_buffer_index, 'hs_force'] = self.hs_force
    self.hs_data.at[self.data_buffer_index, 'cb_force'] = self.myof.cb_force
    self.hs_data.at[self.data_buffer_index, 'pas_force'] = self.myof.pas_force

    if (self.myof.kinetic_scheme == '3state_with_SRX'):
        self.hs_data.at[self.data_buffer_index, 'M_OFF'] = \
            self.myof.y[0]
        self.hs_data.at[self.data_buffer_index, 'M_ON'] = \
            self.myof.y[1]
        self.hs_data.at[self.data_buffer_index, 'M_bound'] = \
            np.sum(self.myof.y[2 + np.arange(self.myof.no_of_x_bins)])
        self.hs_data.at[self.data_buffer_index, 'n_off'] = \
            np.sum(self.myof.y[-2])
        self.hs_data.at[self.data_buffer_index, 'n_on'] = \
            np.sum(self.myof.y[-1])

        # Update fluxes
        fluxes = self.myof.return_fluxes(self.myof.y, self.Ca_conc)
        self.hs_data.at[self.data_buffer_index, 'J1'] = fluxes['J1']
        self.hs_data.at[self.data_buffer_index, 'J2'] = fluxes['J2']
        self.hs_data.at[self.data_buffer_index, 'J3'] = np.sum(fluxes['J3'])
        self.hs_data.at[self.data_buffer_index, 'J4'] = np.sum(fluxes['J4'])
        self.hs_data.at[self.data_buffer_index, 'Jon'] = fluxes['Jon']
        self.hs_data.at[self.data_buffer_index, 'Joff'] = fluxes['Joff']

    if (self.membr.kinetic_scheme == "Ten_Tusscher_2004"):
        # Ten Tusscher membrane voltage is in mV
        self.hs_data.at[self.data_buffer_index, 'membrane_voltage'] = \
            0.001*self.membr.y[0]

    self.hs_data.at[self.data_buffer_index, 'cb_number_density'] = \
        self.cb_number_density
