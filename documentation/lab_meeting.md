Kurtis Mann  
Edited: 2/12/2020

## Goals for 2/24/2020 - 3/2/2020  
- [ ] Get three state revisions to Dr. Wenk by Friday
- [ ] Run rat mesh with fewer elements with new windkessel parameters (today)
- [ ] Follow up e-mail to Russel about Singularity (today)
- [ ] Documentation for 30 minutes daily
- [ ] Add stopping criteria for optimization
- [ ] Implement way to output data from each opt simulation (useful for SOBOL)

## Goals for 2/17/2020 - 2/24/2020
- [x] Make projects on Github for fenics development and sobol sensitivity
- [x] Really work on 3 state revisions
- [ ] Run rat specific mesh with fewer elements, especially on lab machine
- [ ] Tune windkessel parameters for rat mesh
- [ ] Keep in touch with Russel about singularity on cluster
- [x] Work on PSO generalization (specifics on Github)
- [ ] Work on documentation (Doxygen) and instruction markdowns

## Goals for 2/12/2020 - 2/17/2020
- [ ] Add fenics development plans into issue tracker
    * Show lab this for PSO during lab meeting!
- [x] Code own python swarm optimization?
  - [x] Look more into SMT (python surrogate modeling toolbox)- derivatives
- [ ] Run rat specific mesh with fewer elements
- [ ] Get containers working on Amir's old computer/get Newton working
- [x] Work with Russel to get Singularity on cluster
- [x] Submit paper review tonight!!
- [x] Work on 3 state paper revisions
- [x] Notes on nonlinear book
- [ ] Update instruction markdowns

## Goals for 2/5/2020 - 2/12/2020:
- [ ] Add all fenics development plans into issue tracker
- [x] Investigate github project planner capability
- [x] Update some documentation daily
- [x] Begin coding optimization framework (Ask them about LS Opt vs NLopt)
    * LS Opt is familiar to our lab
    * LS Opt has SOBOL sensitivity already built in
    * LS Opt would probably quicker to get working
    * NLopt is open source
    * Not easy to optimize to entire traces, don't know if we can do this in NLopt
    * See notes about optimization in development_tasklist.md
- [x] Finish journal paper review
- [x] Code fenics driver
  -[x] Run elliptical model with driver
  -[x] Implement cell ion drive into single cell script
  -[x] Run single cell with driver script
- [ ] Run simulation with 3 state mesh, tune Windkessel parameters (with fewer parameters)
- [x] Follow up with Russel about singularity

## Progress from 1/29/2020 - 2/5/2020:
- [x] Added some issues to github issue tracker
- [x] Scanned with Hossein
- [x] Finished building a mesh for sensitivity analysis
- [x] Decided to write "driver" script instead of condensing all cases to one fenics file
    * There are just too many differences from case to case
- [x] Created calcium module/framework to connect different cell ion modules
- [x] Looked into optimization framework
- [x] Attempted to run one complete simulation with rat mesh
    * Need to reduce the number of elements, having trouble running complete simulation on laptop
- [x] Worked on 3 state revisions
- [x] Worked on paper review

## Progress from 1/22/2020 - 1/29/2020:  
- [x] Most of the parameters are read in from JSON (had to convert from unicode)  
- [x] Learned mesh generation (first mesh finished, need to speak to Hossein about one step)  
- [x] Worked on fenics script organization  
    - [x] Learned github issue functionality  
    - Code is functional for three-state sensitivity  
    - Moving forward, need to switch from specifying cardiac period and number of cycles to just specifying a simulation duration (makes code general, helps with optimizing passive parameters to end diastole)  
- [x] First pass at reviewing paper  
- [x] Worked on revisions to three state manuscript  
