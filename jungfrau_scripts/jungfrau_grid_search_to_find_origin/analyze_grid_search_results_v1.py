# In order to use this script, first need to grep for the following from the log files dumped by grid_search code
# grep  "Spots within critical pixel RMSD cutoff = 3" out_grid_search_v3/log_rank0* >> v3_results_3_spots.dat
# Then feed that file to this script
# Plot histograms as need be through ipython session
 

import matplotlib.pyplot as plt
all_dx = []
all_dy = []
all_dz = []
# First read in the file which has all the grepped information regarding the Spots within critical pixel RMSD cutoff
with open('v3_results_2_spots.dat') as fin:
  for line in fin:
    if line !='/n':
      ax = line.split()
      dx = float(ax[-3]) 
      dy = float(ax[-2]) 
      dz = float(ax[-1]) 
      all_dx.append(dx)
      all_dy.append(dy)
      all_dz.append(dz)
      
x = plt.hist(all_dx, bins=200, range=(-15, 15))
y = plt.hist(all_dy, bins=200, range=(-15, 15))
z = plt.hist(all_dz, bins=200, range=(-15, 15))

print ('Max dx value = ', x[1][np.argmax(x[0])], )
print ('Max dy value = ', y[1][np.argmax(y[0])], )
print ('Max dz value = ', z[1][np.argmax(z[0])], )

from IPython import embed; embed(); exit()
