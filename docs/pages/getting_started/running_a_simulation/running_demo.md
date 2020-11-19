---
title: Demos
parent: Getting Started
nav_order: 2
has_children: true
has_toc: false
---

MMotH-Vent relies on user-specified information in the form of instruction files. These instruction files use the JSON format and define dictionaries of information to be parsed. Note, the convention is for input values to be specified as elements of a list, e.g. "keyword": [value] rather than "keyword": value. This page gives an overview of the expected keywords of an instruction file and their options. Demos are provided to show some of the current simulations that can be performed.  

How to Run An Instruction File
------------------------------
Make sure Docker is running and enter the command line of your container by following the instructions at the bottom of the [previous page](../installation/installation.md#enter-container-command-line). Then navigate to the directory that contains the instruction file. In general, the syntax for running MMotH-Vent is:

```
python <path to fenics_driver.py> <path to instruction file>
```

Assuming the file structure is the default from the repository, you can run an instruction file by changing the directory to the one containing the instruction file and executing the following:  

```
python /home/fenics/shared/source_code/fenics_driver.py <name of instruction file>
```

<a href="../installation/installation.html" class="btn btn--primary"><< Installation</a>
<a href="../creating_input_files/fenics_input_readme.html" class="btn btn--primary">Building a Mesh >></a>
