import numpy as np

import matplotlib.pyplot as plt

refinement = 5

x1_partition = np.linspace(0,10,refinement)
x2_partition = np.linspace(15,30,refinement)
x3_partition = np.linspace(10,15,refinement)


init_point_dictionary = {"x1": [[0, 10], 0, 0],
    "x2": [[15, 30], 0, 0]
    }

index_and_counter_dict = {"x1": [0, 0], "x2": [0,0]}

iteration_counter = 0

num_particles = 25

swarm = []

for i in np.arange(num_particles):

    temp_working_dict = init_point_dictionary

    for var in init_point_dictionary.keys():

        if iteration_counter < 1:

            if index_and_counter_dict[var][1] < (refinement-1) and index_and_counter_dict[var][0] < (refinement-1):

                index_and_counter_dict[var][0] +=1

                iteration_counter += 1

    for var in temp_working_dict.keys():

        temp_working_dict[var][1] = x1_partition[index_and_counter_dict[var][0]]

    iteration_counter = 0

    swarm.append(temp_working_dict)

x = []
y = []

for i in range(len(swarm)):
    x.append(swarm[i]["x1"][1])
    y.append(swarm[i]["x2"][2])

#print x, y

xnew, ynew, znew = np.meshgrid(x1_partition,x2_partition,x3_partition)

print xnew
#print ynew

#print np.shape(xnew)
#plt.plot(xnew,ynew,marker="o")
xnew = xnew.reshape((np.prod(xnew.shape)))
print xnew

#print np.shape(xnew)
ynew = ynew.reshape((np.prod(ynew.shape)))
print ynew
znew = znew.reshape((np.prod(znew.shape)))

coords = zip(xnew,ynew,znew)

#print coords
#print np.shape(coords)

plt.plot(xnew,ynew,marker="o")

plt.show()
