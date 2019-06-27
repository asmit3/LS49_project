import numpy as np

x=np.random.rand(20)
y=np.random.rand(20)
z=np.random.rand(20)
a = np.random.random((16, 16))

import matplotlib.pyplot as plt

plt.contourf(x,y,z,cmap=plt.cm.rainbow)
plt.show()
