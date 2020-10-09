import pypdm

import numpy as np
import matplotlib.pyplot as plt

t = np.linspace(1, 20, 200)
y = np.sin(t)
s = np.zeros(200)

freq, theta = pypdm.pdm(200, t, y, s, 0.01, 1, 0.001)

plt.plot(freq, theta, 'k')
plt.axvline(1/2/np.pi, color='red')
for i in range(2,4):
    plt.axvline(1/2/i/np.pi, color='red', ls='--')
plt.xlabel('Frequency')
plt.ylabel('Theta')
plt.show()
