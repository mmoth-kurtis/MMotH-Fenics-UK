import numpy as np
import math
simulation_dir = '/Users/charlesmann/Academic/UK/MMotH-Fenics-UK/working_directory_untracked/082315/coarse/output/'
stress = np.load(simulation_dir + 'stress_array.npy')
nan_vector = np.zeros(np.shape(stress))
checksum = 0
for i in range(0,1701):
    for j in range(0,20174):
        nan_vector[i,j]=math.isnan(stress[i,j])
        checksum = checksum + nan_vector[i,j]

print checksum
