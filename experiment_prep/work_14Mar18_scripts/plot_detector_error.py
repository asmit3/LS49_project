# Script to plot dr(spots) vs the error in detector dimentsion (dx,dy)
import os
import numpy as np
from dials.array_family import flex

directories = [x[0] for x in os.walk('../') if 'dir' in x[0]]
# Detector configuration 
x0=-364.892 # mm
y0=38.4 # mm
z0=-331.0 # mm
pixel_size=0.11 #mm
nx=32*2
ny=32*2

error_cutoff = 10.0

dx=[]
dy=[]
data = np.zeros((nx,ny)) + error_cutoff
for directory in directories:
#  dx.append(float(directory.split("_")[-3]))
#  dy.append(float(directory.split("_")[-2]))
  dx = int(directory.split("_")[-3])
  dy = int(directory.split("_")[-2])
  try:
    fin = open(os.path.join(directory,'out.log'))
  except Exception:
    continue
  for line in fin:
    if line !='\n' and 'AVERAGE_ERROR_IN_SPOT_DIST' not in line:
      dR=line
  dR = dR[1:-2]
  try:
    x1=float(dR.split(',')[0])
    y1=float(dR.split(',')[1])
#    z1=float(dR.split(',')[2])
  except:
    #from IPython import embed; embed(); exit()
    print 'SOMETHING WRONG in ',directory
    continue
  detector_error_arr = flex.double([x1-x0, y1-y0])
  detector_error = flex.sum(detector_error_arr**2)
  #from IPython import embed; embed(); exit()

  #print 'dR = ',dR,dx,dy
  #data.append((dx,dy,dR))
  if np.sqrt(detector_error) > error_cutoff:
    data[dx+int(nx/2),dy+int(ny/2)] = error_cutoff
  else:
    data[dx+int(nx/2),dy+int(ny/2)]= np.sqrt(detector_error)
#exit()
from libtbx import easy_pickle

easy_pickle.dump('error_data.pickle', data)

exit()
import matplotlib.pyplot as plt
from matplotlib import cm
import numpy as np

#print data
#from IPython import embed; embed(); exit()
from scitbx.math import five_number_summary
d = data.reshape(nx*ny)
plt.hist(d, bins=200)
plt.xticks(xrange(0,800,10))
print 'FIVE NUMBER SUMMARY = ',five_number_summary(d)

fig,ax = plt.subplots()
cax = ax.imshow(data, cmap=cm.coolwarm, interpolation='nearest', origin='lower', extent = [-int(nx/2),int(nx/2),-int(ny/2), int(ny/2)])
ax.set_title('Error in detector configuration after minimization')
#from IPython import embed; embed()
#ax.set_xticks(range(-32,33))
#ax.set_yticks(range(-32,33))
cbar = fig.colorbar(cax, ticks=[0.0, 0.2, 0.5])
cbar.ax.set_yticklabels(['0.1', '0.5mm', '2.0 mm','3','4','5','6','7','10']) 
#cbar = fig.colorbar(cax, ticks=[0.0, 1.0, 2.0, 3.0,4.0,5.0,6.0,7.0,10.0])
#cbar.ax.set_yticklabels(['0.1', '0.5mm', '2.0 mm','3','4','5','6','7','10']) 
plt.show()






