import sys
import json
import os
import dependencies
from dependencies import recode_dictionary
from dependencies import load_parameters
import numpy as np
sys.path.append("/home/fenics/shared/source_code/")
## This should be running in a FEniCS container, make it easy to import necessary source code

## This is the simulation driver.
#
#Loads in parameters and selects appropriate fenics script to run
# User provides file name
input_file_name = sys.argv[1]

# Check that the file exists, if not , exit
if not os.path.exists(input_file_name):
    print "input file does not exist"
    exit()

# Load in JSON dictionary
with open(input_file_name, 'r') as json_input:
  input_parameters = json.load(json_input)

# Convert any unicode values to python strings so they work with some cpp libraries
recode_dictionary.recode(input_parameters)

## Parse out the different types of parameters
sim_params = input_parameters["simulation_parameters"]
file_inputs = input_parameters["file_inputs"]
output_params = input_parameters["output_parameters"]
passive_params = input_parameters["forms_parameters"]["passive_law_parameters"]
hs_params = input_parameters["myosim_parameters"]
cell_ion_params = input_parameters["electrophys_parameters"]["cell_ion_parameters"]
monodomain_params = input_parameters["electrophys_parameters"]["monodomain_parameters"]
windkessel_params = input_parameters["windkessel_parameters"]

## Assign input/output parameters
output_path = output_params["output_path"][0]
input_path = file_inputs["input_directory_path"][0]

# This may only be needed in ventricle simulations
casename = file_inputs["casename"][0]

# Check that the output path exists. If it does not, create it and let user know
if not os.path.exists(output_path):
    print "Output path does not exist. Creating it now"
    os.makedirs(output_path)

# Figure out which script needs to be executed
if sim_params["sim_geometry"][0] == "ventricle":
    fenics_script = "fenics_LV"
    print fenics_script

elif sim_params["sim_geometry"][0] == "single_cell":
    fenics_script = "fenics_singlecell_isometric"

# For now, going to have to import appropriate script as a module and wrap the whole
# thing in a function so that inputs and output can be passed
script_name = __import__(fenics_script)

# Call the "fenics" function within the script
output_dictionary = script_name.fenics(sim_params,file_inputs,output_params,passive_params,hs_params,cell_ion_params,monodomain_params,windkessel_params)

# Save the appropriate output information
np.save(output_path + "rates",output_dictionary["rates"])
np.save(output_path + "dumped_populations",output_dictionary["dumped_populations"])
np.save(output_path + "tarray",output_dictionary["tarray"])
np.save(output_path + "stress_array",output_dictionary["strarray"])
np.save(output_path + "pstress_array",output_dictionary["pstrarray"])
np.save(output_path + "alpha_array",output_dictionary["alphaarray"])
np.save(output_path + "calcium",output_dictionary["calarray"])
np.save(output_path + "hsl",output_dictionary["hsl"])

# If user wants to visualize, do that here
#if output_params["visualize_flag"] > 0:

    # call plotting script
 #   print "going to visualize"
