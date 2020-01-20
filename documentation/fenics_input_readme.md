This document explains each input in the "fenics_input.json" file. Terms in "" represent the parameter name, and terms in [] represent the possible values.
Note: All inputs must be in brackets [], as the code loads in "key"[value][0]. This was to be able to include units next to input values, and standardize input.

Author: Kurtis Mann  
Date Edited: 1/16/2010  


<h1>Simulation Parameters"</h1>  
These define what drives the simulation.  

 **"sim_geometry": [single_cell], [sheet], [ventricle]**  
 Choose between [single_cell], [sheet], or [ventricle].
    [single_cell]: Simulates a unit cube meant to represent a single half sarcomere.  
    [sheet]: A 2D thin sheet of sarcomeres.  
    [ventricle]: Full ventricle simulations. Can choose specific meshes or use ellipsoidal.  

  **"sim_type"**    
  For the given geometry, choose what type of simulation.  
  
    Single Cell Simulations  
    
    [isometric]: Unit cube where length is fixed.  
    [ktr]: Tension recovery experiment.  
    [force]: Unit cube where a given force boundary condition is enforced.  
    [single_cell_custom]: Allows for a combination of length and force controlled simulations. Also allows for length changes.  

    Tissue Sheet Simulations  
    Place-holder for now. Will use this for fiber splay simulations or spiral dynamics.  
  
    Left Ventricle Simulations  
      
    [beat, "path/to/mesh"]: Supply the (relative?) path to the desired mesh and it will "beat" for "sim_duration"  
    [vena_cava_occlusion, "path/to/mesh"]: Use the specified mesh and simulate a vena cava occlusion. Useful for calculating ESPVR.  
  
 <h2>"sim_duration": [*integer*]</h2>  
 Enter an *integer* time for the simulation in seconds.  
  
 <h2>"sim_timestep": [*real*]</h2>  
 Time-step to be used in implicit finite element solving.  

<h1>"output_parameters"</h1>  
These control what information gets output for a given simulation.  

<h1>"forms_parameters"</h1>  
Input all parameters needed for the "forms.py" file. Forms is a class object that contains all of the methods needed to calculate necessary continuum tensors, as well as return stresses from constitutive equations.  

  "passive_law": Select one of the available passive tissue laws.  
    [guccione_transverse_isotropy]: Phenomenological model. See "". Must define C, bf, bt, bfs.  
      [c]:  
      [bf]:  
      [bt]:  
      [bfs]:  
    [semi_structural]: Separates out the myofiber passive stress from bulk response. Requires Guccione parameters AND c2, c3. See above and "Microstructure-based finite element model of left ventricle passive inflation" by Xi, et al.  
      [c]:  
      [c2]:  
      [c3]:  
      [bf]:  
      [bt]:  
      [bfs]:  
    [full_structural]: Plan to implement a fully structural model.  

  "Kappa"

<h1>"myosim_parameters"</h1>  
These are the active contraction parameters to be used by the embedded python version of MyoSim.  
