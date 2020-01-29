This document shows planned changes/improvements to the FEniCS framework. This will be updated at least weekly to show progress.  
Author: Kurtis Mann  
Date Edited:  1/22/2020  

**Documentation**
- [ ] Finish .md file for JSON input
- [ ] Create a .md manual/documentation for fenics script
- [ ] Update .md file for creating meshes
- [ ] Revise/improve instructions for running docker/fenics
- [ ] Comment ALL code and figure out Doxygen  
- [ ] Put all markdown files in doxygen generated documentation

**FEniCS Functionality**  
- [ ] **Organize fenics script so that all existing cases (singe cell isometric, single cell force bc, ellipsoid LV, patient specific mesh LV)  are in one file.** Currently, the user has to specify the cardiac period and number of cycles to simulate. This will work for the sensitivity analyis, but moving forward I need to implement a way to specify a general simulation time and allow for variable period?
- [ ] **Complete json input file so that the above cases can be switched between without changing anything in the fenics script.**  
- [ ] **Have small json input file for mesh generation, be able to specify output directory and mesh name**
- [ ] **Start coding up optimization wrapper**
  - Have messaged Dylan about his optimization repository.
  - Dylan has shared his repository with me. This will be a great help. Speaking to him about running jobs in parallel
- [ ] Generalize Python Myosim schemes  


**Sobol Sensitivity**
- [x] Get MRI data from 3 state study
  * The raw data is in the "austin_backup" directory on the cluster
- [x] Learn how to create meshes from Hossein (Friday)
- [ ] **Create meshes for rat LV's for three state sensitivity analysis**
    - First one created, has some issues
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
- [ ] Notes on 3.1-3.4 of nonlinear book  

**Cell Ion Models**
- [ ] Modularize cell ion/calcium transient  

**Monodomain**
- [ ] Look at monodomain (already in FEniCS)  
