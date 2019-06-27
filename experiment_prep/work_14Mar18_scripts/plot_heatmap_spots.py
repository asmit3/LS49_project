# Script to plot dr(spots) vs the error in detector dimentsion (dx,dy)
import os
import numpy as np

directories = [x[0] for x in os.walk('../') if 'dir' in x[0]]
print directories
dx=[]
dy=[]
data = np.zeros((32,32)) + 5.0
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
    if line !='\n' and 'AVERAGE_ERROR_IN_SPOT_DIST' in line:
      dR=float(line.split("=")[-1])

  #print 'dR = ',dR,dx,dy
  #data.append((dx,dy,dR))
  data[dx,dy]=int(dR/.11)*0.11
#exit()

import matplotlib.pyplot as plt
import numpy as np

#print data
#from IPython import embed; embed(); exit()

fig,ax = plt.subplots()
cax = ax.imshow(data, cmap='hot', interpolation='nearest', origin='lower')
ax.set_title('Average error between predicted and observed spot distance after minimization')
cbar = fig.colorbar(cax, ticks=[0, 1.0, 3.0])
cbar.ax.set_yticklabels(['0', '1.0mm', '3.0 mm']) 
plt.show()


