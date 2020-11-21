---
title: Cell Isometric Twitch
parent: Demos
grand_parent: Getting Started
nav_order: 1
---

Summary
-------
A unit cube mesh consisting of six tetrahedral elements is used to model a single muscle cell tension response to a twitch calcium transient. The instruction file is included in the repository, and can be downloaded <a href="https://github.com/mmoth-kurtis/MMotH-Vent/blob/master/demos/cell_isometric_twitch_demo/cell_isometric_twitch_demo.json" >here</a>

Simulation Protocol
-------------------
The muscle cell is stretched 11.5% over 5 ms so that the half-sarcomere length increases from 950 nm to 1100 nm. This stretch is maintained for 15 ms at a calcium concentration of 1e-7 M to allow the cross-bridges to reach steady state. Then the cell is activated with a skeletal muscle calcium transient approximation from (citation) using the [two-compartment calcium](../../../model_formulations/calcium_models/two_compartment_model/two_compartment_model.md) model, and stretch held fixed.

Boundary Conditions & Assumptions
---------------------------------
- Left face displacement is fixed in the x-direction.
- Right face displacement is prescribed to give the desired magnitude of stretch.
- A single point on the left face is completely fixed to prevent rigid body translation.
- The edges on the ends are fixed in either y or z to allow expansion/compression due to incompressibility while keeping the cross-section area square and prevent rigid-body rotation.

Results
-------
If plotted using the k_plotter.py file, the cell-level model results are below:

<video width="800" height="500" controls>
  <source src="test_animation.mp4" type="video/mp4">
</video>


And the cube deformation is seen in Paraview, shown below. The color shows the magnitude of active stress generated in the f0 direction (the direction shown by the vectors in the bottom image):

<video width="800" height="500" controls>
  <source src="displacement_animation.m4v" type="video/mp4">
</video>
<img src="https://github.com/mmoth-kurtis/MMotH-Vent/blob/master/docs/pages/getting_started/running_a_simulation/cell_isometric_demo_page/f0_cell_isometric_demo_2.png?raw=true" width="800" height="500">
