import sys
sys.path.append("/Users/charlesmann/Academic/UK/fenics/source_code/dependencies")

from cell_ion_module import cell_ion_driver
import json
# Different path because not running in fenics yet
import recode_dictionary
# Something is not right when using this
# Test it out
input_file_name = sys.argv[1]
with open(input_file_name, 'r') as json_input:
  input_parameters = json.load(json_input)
recode_dictionary.recode(input_parameters)
cell_ion_params = input_parameters["electrophys_parameters"]["cell_ion_parameters"]

cell_ion = cell_ion_driver.cell_ion_driver(cell_ion_params)
temp = cell_ion.model.calculate_concentrations(1,.05)
print temp
