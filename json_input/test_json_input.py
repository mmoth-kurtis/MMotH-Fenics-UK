import json
import sys
sys.path.append("/Users/charlesmann/Academic/UK/fenics/source_code/dependencies")
import recode_json_strings as rc
import recode_dictionary
def main():
  # Define the file name. If this was in another directory, I would use the os.path.join() function.
  # It's pretty handy.
  input_file_name = sys.argv[1]


  # Read in the JSON tree structure.
  with open(input_file_name, 'r') as f:
    json_tree = json.load(f)
    #sim_params = json_tree["simulation_parameters"]
  # This is another way to read it in but I find the above to be simpler.
  # f = open("example_1.json", 'r')
  # json_tree = json.load(f)
  # f.close()
  recode_dictionary.recode(json_tree)
  
  """for v in sim_params.values():

    #for value in v.values():


    counter = 0
    for j in v:
        print v
        print counter
        print j
        if type(j) is unicode:
            rcj = rc._byteify(j)
            print counter
            v[counter] = rcj
            print v[counter]
            print type(v[counter])
        counter +=1
                #value[counter] = rcj
                #print value
            #counter += 1"""
          #file_inputs[key]=[value]
  #print file_inputs
  #print file_inputs
  #for key, value in file_inputs.items():
      #for j in value:
          #print type(j)
  """for value in file_inputs.values():
      for j in value:
          print type(j)"""
  #casename.type()
  #recoded_casename = rc._byteify(casename)
  #recoded_casename.type()

  # Add another optional parameter.
  #json_tree["optional_parameters"]["super_optional"] = "dumb"

  # Define the output file name. Again, if you wanted this in a directory called, "outputs", I would
  # Do the following:
  #  output_file_name = os.path.join("outputs", "example_1_output.json")
  #  output_file_name = "example_1_output.json"

  # Dump the JSON file. We use the `indent` argument to specify we want the JSON object to have 2
  # spaces between objects. This makes it far more readable.
  #with open(output_file_name, 'w') as f:
    #json.dump(json_tree, f, indent=2)


if __name__ == "__main__":
  main()
