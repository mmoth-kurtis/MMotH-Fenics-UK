This document shows planned changes/improvements to the FEniCS framework. This will be updated at least weekly to show progress.  
Author: Kurtis Mann  
Date Edited:  1/22/2020  
  
**Documentation**
- [ ] Finish .md file for JSON input
- [ ] Create a .md manual/documentation for fenics script
- [ ] Create .md file for creating meshes
- [ ] Revise/improve instructions for running docker/fenics
- [ ] Comment ALL code and figure out Doxygen  

**Miscellaneous**
- [x] Set up GitHub repository
- [x] Initialize Doxygen documentation  
  * To utilize the special commands of Doxygen, ## instead of document strings """ must be used. Going to take a bit, but want to document every file.  

**FEniCS Functionality**  
- [ ] Organize fenics script so that all existing cases (singe cell isometric, single cell force bc, ellipsoid LV, patient specific mesh LV)  are in one file. Also just organize in general
- [ ] Complete json input file so that the above cases can be switched between without changing anything in the fenics script.
**Sobol Sensitivity**
- [ ] Learn how to create meshes from Hossein
- [ ] Create meshes for rat LV's for three state sensitivity analysis
- [ ] Start coding up optimization wrapper
- [ ] Optimize parameters for sensitivity analysis
- [ ] Code up Sobol' indices sensitivity framework

