import numpy as np
import matplotlib.pyplot as plt

d = np.genfromtxt('out.csv')
plt.plot(d[0,:],d[1,:])
plt.show()
