This document explains each input in the "fenics_input.json" file. Terms in "" represent the parameter name, and terms in [] represent the possible values.
Note: All inputs must be in brackets [], as the code loads in "key"[value][0]. This was to be able to include units next to input values, and standardize input.

Author: Kurtis Mann  
Date Edited: 1/20/2020  


### **Simulation Parameters**
**"sim_geometry"**: Choose between [single_cell], [sheet], or [ventricle].  
[single_cell]: Simulates a unit cube meant to represent a single half sarcomere.  
[sheet]: A 2D thin sheet of sarcomeres.  
[ventricle]: Full ventricle simulations. Can choose specific meshes or use ellipsoidal.  

**"sim_type"**: For the given geometry, choose what type of simulation.  

**Single Cell** | Description   
----------------|-------  
[isometric]:| Unit cube where length is fixed.  
[ktr]:| Tension recovery experiment.  
[force]:| Unit cube where a given force boundary condition is enforced.  
[single_cell_custom]:| Allows for a combination of length and force controlled simulations. Also allows for length changes.  

Tissue Sheet | Description
-------------|---------------
[sheet]:| Place-holder for now. Will use this for fiber splay simulations or spiral dynamics.  

Left Ventricle | Description
---------------|--------------
[beat, "path/to/mesh"]:| Supply the (relative?) path to the desired mesh and it will "beat" for "sim_duration"  
[vena_cava_occlusion, "path/to/mesh"]:| Use the specified mesh and simulate a vena cava occlusion. Useful for calculating ESPVR.  

**"sim_duration":** Enter an *integer* time for the simulation in seconds.  

**"sim_timestep":** Time-step to be used in implicit finite element solving.  

### **Output Parameters**  

### **Forms Parameters**
**"passive_law":** Select one of the available passive tissue laws.  
* [guccione_transverse_isotropy]: Phenomenological model. See "". Must define C, bf, bt, bfs.  
    * [c]:  
    * [bf]:  
    * [bt]:  
    * [bfs]:  
* [semi_structural]: Separates out the myofiber passive stress from bulk response. Requires Guccione parameters AND c2, c3. See above and "Microstructure-based finite element model of left ventricle passive inflation" by Xi, et al.  
    * [c]:  
    * [c2]:  
    * [c3]:  
    * [bf]:  
    * [bt]:  
    * [bfs]:  
* [full_structural]: Plan to implement a fully structural model.
  "Kappa"

### **MyoSim Parameters**  
**""max_rate":** [5000,"s^-1"]  
**"temperature":** [288, "Kelvin"]  
**"cb_number_density":** [7.67e16, "number of cb's/m^2"]  
**"initial_hs_length":** [1000, "nm"]  

# **Myofilament Parameters**  
  **"kinetic_scheme":** ["3state_with_SRX"]  
  **"num_states":** [3]  
  **"num_attached_states":** [1]  
  **"num_transitions":** [4]  
  **"cb_extensions":** [[0.0, 0.0, 4.75642], "power-stroke distance in nm"]  
  **"state_attached":** [[0, 0, 1]]  
  **"k_cb_multiplier":** [[1.0, 1.0, 1.0]]  
  **"k_cb_pos":** [0.001, "N*m^-1"]  
  **"k_cb_neg":** [0.001, "N*m^-1"]  
  **"alpha":**[1.0]  
  **"k_1":** [9.623166, "s^-1"]  
  **"k_force":** [1.96345e-4, "(N^-1)(m^2)"]  
  **"k_2":** [1000, "s^-1"]  
  **"k_3":** [5435.531288, "(nm^-1)(s^-1)"]  
  **"k_4_0":** [2543.864648, "s^-1"]  
  **"k_4_1":** [0.18911849, "nm^-4"]  
  **"k_cb":** [0.001, "N*m^-1"]  
  **"x_ps":** [4.75642, "nm"]  
  **"k_on":** [1.5291356e8, "(M^-1)(s^-1)"]  
  **"k_off":** [100, "s^-1"]  
  **"k_coop":** [6.38475]  
  **"bin_min":** [-12, "nm"]  
  **"bin_max":** [12, "nm"]  
  **"bin_width":** [0.5, "nm"]  
  **"filament_compliance_factor":** [0.5]  
  **"thick_filament_length":** [815, "nm"]  
  **"thin_filament_length":** [1120, "nm"]  
  **"bare_zone_length":** [80, "nm"]  
  **"k_falloff":** [0.0024]  
  **"passive_mode":** ["exponential"]  
  **"passive_exp_sigma":** [500]  
  **"passive_exp_L":** [80]  
  **"passive_l_slack":** [900, "nm"]  
