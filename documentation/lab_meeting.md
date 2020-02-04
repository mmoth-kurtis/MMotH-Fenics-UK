2/5/2020

Goals for 2/5/2020 - 2/12/2020:'
- [ ] Add all fenics development plans into issue tracker
- [ ] Investigate github project planner capability
- [ ] Update some documentation daily
- [ ] Begin coding optimization framework (Ask them about LS Opt vs NLopt)
    * LS Opt is familiar to our lab
    * LS Opt has SOBOL sensitivity already built in
    * LS Opt would probably quicker to get working
    * NLopt is open source
    * Not easy to optimize to entire traces, don't know if we can do this in NLopt
- [ ] Finish journal paper review
- [ ] Code fenics driver 
- [ ] Run simulation with 3 state mesh, tune Windkessel parameters

Progress from 1/29/2020 - 2/5/2020:
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

Progress from 1/22/2020 - 1/29/2020:  
- [x] Most of the parameters are read in from JSON (had to convert from unicode)  
- [x] Learned mesh generation (first mesh finished, need to speak to Hossein about one step)  
- [x] Worked on fenics script organization  
    - [x] Learned github issue functionality  
    - Code is functional for three-state sensitivity  
    - Moving forward, need to switch from specifying cardiac period and number of cycles to just specifying a simulation duration (makes code general, helps with optimizing passive parameters to end diastole)  
- [x] First pass at reviewing paper  
- [x] Worked on revisions to three state manuscript  



Goals for next week:  
- [ ] Add all "feature tracking" to github issue tracker  
- [ ] Investigate Github project planning capability  
- [ ] Build all meshes for three state sensitivity  
- [ ] Update all documentation  
- [ ] Have made significant progress on coding optimization framework  
- [ ] Significant progress on three state revisions  
- [ ] Finish paper review  
- [ ] Start thinking about calcium transient for sensitivity paper  
  
Questions:  
- [x] Hossein: can you help me with the transform world matrix mesh generation step?  
- [x] General: Tuning Windkessel parameters?   
