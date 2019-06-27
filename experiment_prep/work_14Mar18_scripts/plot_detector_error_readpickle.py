# Script to plot dr(spots) vs the error in detector dimentsion (dx,dy)
import os
import numpy as np
from dials.array_family import flex

# Detector configuration 
x0=-364.892 # mm
y0=38.4 # mm
z0=-331.0 # mm
pixel_size=0.11 #mm
nx=32*2
ny=32*2

error_cutoff = 10.0

from libtbx import easy_pickle

data = easy_pickle.load('error_data.pickle')

import matplotlib.pyplot as plt
from matplotlib import cm
import numpy as np

from scitbx.math import five_number_summary
d = data.reshape(nx*ny)
plt.hist(d, bins=200)
plt.xticks(xrange(0,800,10))
print 'FIVE NUMBER SUMMARY = ',five_number_summary(d)

fig,ax = plt.subplots()
cax = ax.imshow(data, cmap=cm.coolwarm, interpolation='nearest', origin='lower', extent = [-int(nx/2),int(nx/2),-int(ny/2), int(ny/2)])
ax.set_title('Error in detector configuration after minimization')
ax.set_xlabel('Detector perturbation in x-axis[mm]')
ax.set_ylabel('Detector perturbation in y-axis[mm]')
#from IPython import embed; embed()
#ax.set_xticks(range(-32,33))
#ax.set_yticks(range(-32,33))
cbar = fig.colorbar(cax, ticks=[0.0, 1.0, 2.0, 3.0,4.0,5.0,6.0,7.0,10.0])
cbar.ax.set_yticklabels(['0.1', '0.5mm', '2.0 mm','3','4','5','6','7','10']) 
plt.show()






