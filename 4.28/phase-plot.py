import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt('cm3-1.txt', skiprows=1)
plt.plot(data[:, 0], data[:, 1], label='21')
plt.plot(data[:, 0], data[:, 2], label='22')
plt.plot(data[:, 0], data[:, 2] - data[:, 1], label='diff')
plt.legend()
plt.show()
