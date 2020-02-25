This document shows planned changes/improvements to the FEniCS framework. This will be updated at least weekly to show progress.  
Author: Kurtis Mann  
Date Edited:  1/22/2020  

**Documentation**
- [ ] Finish .md file for JSON input (need to add descriptions for myosim parameters)
- [ ] Create a .md manual/documentation for fenics script
- [ ] Update .md file for creating meshes
- [ ] Revise/improve instructions for running docker/fenics
- [ ] Comment ALL code and figure out Doxygen  
- [ ] Put all markdown files in doxygen generated documentation

**FEniCS Functionality**  
- [ ] ~~Organize fenics script so that all existing cases (singe cell isometric, single cell force bc, ellipsoid LV, patient specific mesh LV)  are in one file.** Currently, the user has to specify the cardiac period and number of cycles to simulate. This will work for the sensitivity analyis, but moving forward I need to implement a way to specify a general simulation time and allow for variable period?~~
- [x] Write "fenics_driver" where this script runs the appropriate file depending on the simulation you want to run. There are too many differences between running single cell vs. ventricle to put this all in one script.
- [x] Complete json input file so that the above cases can be switched between without changing anything in the fenics script.  
- [ ] **Have small json input file for mesh generation, be able to specify output directory and mesh name**
    - [ ] **Test this now**

**Optimization**
- [ ] **Start coding up optimization wrapper**
  - Have messaged Dylan about his optimization repository.
  - Dylan has shared his repository with me. This will be a great help. Speaking to him about running jobs in parallel
  - Looking at three approaches:  
    1) Use NLopt. This is free/open source and allows the user to easily switch between different optimization algorithms.
    2) Make FEniCS talk to LS-Opt. We are familiar with this program and it is powerful and flexible. However, it is not open source. Could decrease time to get optimiztaion functional though.
    3) Particle Swarm
- [ ] Generalize Python Myosim schemes  
  - Some notes: We need optimization to be relatively quick. Is there a way to leverage optimizations at the cellular level to inform us about parameters at organ level? For example, optimize cell level parameters to fix some, then optimize fewer parameters at organ level? Also get steady state information so we can optimize to one cardiac cycle at organ level? There should be a small range of what the stress and such should look like at the cell level. Also, use sensitivity analysis to inform us about the jump from single cell to organ. For example, fiber angle will probably play a large role. Will we assign the same parameters to each cell, or different based on transmural location? Will we need to optimize cell level stuff in such a way to account for fiber orientation effects on organ pressure? Do we include hill coefficient in objective function or is this forcing the solution too much?

**Sobol Sensitivity**
- [x] Get MRI data from 3 state study
  * The raw data is in the "austin_backup" directory on the cluster
- [x] Learn how to create meshes from Hossein (Friday)
- [x] Create mesh for rat LV's for three state sensitivity analysis
    - First one created, ~~has some issues~~ resolved.
- [ ] Decide on calcium transient to be used (match 3 state paper?)
- [ ] Optimize parameters for sensitivity analysis
    - Insert parameters from 3 state paper, then tune windkessel
    - Then iterate between optimizing passive and active parameters
- [ ] Code up Sobol' indices sensitivity framework
- [ ] Finalize model inputs and outputs for analysis  

**Miscellaneous**
- [x] Set up GitHub repository
- [x] Initialize Doxygen documentation
  * To utilize the special commands of Doxygen, ## instead of document strings """ must be used. Going to take a bit, but want to document every file.
- [ ] **Finish paper review**
- [ ] **Notes/discuss 3 state paper revisions**
- [ ] Notes on chapter 2 of nonlinear book  

**Cell Ion Models**
- [x] Modularize cell ion/calcium transient  
- [ ] Test Shannon model (eventually)

**Monodomain**
- [ ] Look at monodomain (already in FEniCS)  
