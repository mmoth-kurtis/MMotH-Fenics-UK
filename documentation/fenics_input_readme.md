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

**Tissue Sheet | Description**
---------------|---------------
[sheet]:| Place-holder for now. Will use this for fiber splay simulations or spiral dynamics.  
  
**Left Ventricle | Description**  
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

### **Myosim Parameters**
  
