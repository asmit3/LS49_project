# Needs a spot_count.pickle file which is essentially the sum of all the spot counts from the 
# individual pickle files dumped by the grid search code. 
# Some assumptions made here about the grid search that was done i.e 100 bins, 0.15 mm grid size etc so please adapt if used for other cases
# Please see work in /global/cscratch1/sd/asmit/LS49/LS49_SAD_v3/diffBragg_refinement/jungfrau_grid_plot_v1

from libtbx.easy_pickle import load
import numpy as np

spot_count=np.array(load('spot_count.pickle'))
# Set up the grid that was used in the grid search
dx=[0.0]
for ix in range(1,100+1):
  dx.append(ix*0.15)
  dx.append(-ix*0.15)

all_r=[]
for x in dx:
  for y in dx:
    for z in dx:
      all_r.append((x,y,z))

print (len(all_r))

t = np.argsort(spot_count)[-100:]
x=[]
y=[]
z=[]
n=[]

for tt in t:
  x.append(all_r[tt][0])
  y.append(all_r[tt][1])
  z.append(all_r[tt][2])
  n.append(spot_count[tt])

n2 = n/max(n)


print ('Best grid point = ', x[-1], y[-1], z[-1], n[-1])
print ('2nd Best grid point = ', x[-2], y[-2], z[-2], n[-2])
print ('3rd Best grid point = ', x[-3], y[-3], z[-3], n[-3])

from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as  plt
fig=plt.figure()
ax=fig.gca(projection='3d') 
cmap=plt.cm.Oranges

ax.scatter3D(x,y,z,c=n, cmap=cmap)
plt.xlabel('X')
plt.ylabel('Y')
plt.show()

#from IPython import embed; embed(); exit()
