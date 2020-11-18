---
title: Demos
parent: Getting Started
nav_order: 2
has_children: true
---
Overview
========

MMotH-Vent relies on user-specified information in the form of instruction files. These instruction files use the JSON format and define dictionaries of information to be parsed. Note, the convention is for input values to be specified as elements of a list, e.g. "keyword": [value] rather than "keyword": value. This page gives an overview of the expected keywords of an instruction file and their options. Demos are provided to show some of the current simulations that can be performed.  

Make sure Docker is running and enter the command line of your container by following the instructions at the bottom of the [previous page](../installation/installation.md#enter-container-command-line).

How to run an instruction files
------------------------------

Navigate to ```/home/fenics/shared/demos/```, where you will see two directories:  
* ```single_cell_length_control_demo/```
* ```ellipsoidal_ventricle_demo/```

# Single Cell Demo
Enter the ```single_cell_length_control_demo/``` directory, which contains the instruction file ```singlecell_demo.json```. In general, the syntax to execute MMotH-Vent is
```
python [path to fenics_driver.py] [path to instruction file]
```
Therefore to execute this instruction file, use the following command:  

```
python ../../source_code/fenics_driver.py singlecell_demo.json
```

The default option in this instruction file is to save all of the output in the same directory from which MMotH-Vent is called. MMotH-Vent outputs a number of numpy arrays that contain cell-level information. The cell data at any gauss point can be plotted by executing
```
python [path to k_plotter] [integer gauss point]
```
Assuming the path hierarchy is unchanged from the master repository, after running the demo the MyoSim data can be visualized by
```
python /home/fenics/shared/plot_tools/k_plotter 0
```
which should yield the following plot:
<img src="https://github.com/mmoth-kurtis/MMotH-Vent/blob/master/docs/assets/images/Screen%20Shot%202020-07-01%20at%205.03.18%20PM.png?raw=true" alt="titlepage" width="800"/>  

Finally, the mesh and displacement solutions are saved as ParaView files. It is recommended to download Paraview to view these solutions.

# Ventricle Demo
Under construction

<a href="../installation/installation.html" class="btn btn--primary"><< Installation</a>
<a href="../creating_input_files/fenics_input_readme.html" class="btn btn--primary">Building a Mesh >></a>
