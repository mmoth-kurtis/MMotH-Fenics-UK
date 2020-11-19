---
title: Cell Isometric Twitch
parent: Demos
grand_parent: Getting Started
nav_order: 1
---

Summary
-------
A unit cube mesh consisting of six tetrahedral elements is used to model a single muscle cell tension response to a twitch calcium transient.

Simulation Protocol
-------------------
The muscle cell is stretched 11.5% over 5 ms so that the half-sarcomere length increases from 950 nm to 1100 nm. This stretch is maintained for 10 ms at a calcium concentration of 1e-7 M to allow the cross-bridges to reach steady state. Then the cell is activated with a skeletal muscle calcium transient approximation using the [two-compartment calcium](../../../model_formulations/calcium_models/two_compartment_model/two_compartment_model.md) model.

Boundary Conditions & Assumptions
---------------------------------
- Left face displacement is fixed in the x-direction.
- Right face displacement is prescribed to give the desired magnitude of stretch.
- A single point on the left face is completely fixed to prevent rigid body translation.
- The edges on the ends are fixed in either y or z to allow expansion/compression due to incompressibility while keeping the cross-section area rectangular.

Results
-------
If plotted using the k_plotter.py file, the cell-level model results are below:

* Insert plots

And the cube deformation is seen in Paraview here:
