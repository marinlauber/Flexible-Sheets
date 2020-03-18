import numpy as np
import matplotlib.pyplot as plt

d = np.genfromtxt('out.csv')
plt.plot(d[1,:],d[0,:])
plt.show()
