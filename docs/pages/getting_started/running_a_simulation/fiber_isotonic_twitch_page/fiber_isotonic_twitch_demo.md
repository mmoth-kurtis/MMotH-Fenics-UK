---
title: Fiber Isotonic Twitch With Compliance
parent: Demos
grand_parent: Getting Started
nav_order: 4
---

Summary
-------
A muscle fiber twitch simulation. A cylinder geometry in which the ends represent non-contracting, fibrous tissue with a nonlinear stiffness and the quadrature points represent contracting myofibrils.  

Simulation Protocol
-------------------
The right end of the fiber is displaced 

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
