import pypdm

import numpy as np
import matplotlib.pyplot as plt

t = np.linspace(1, 20, 200)
y = np.sin(t)
s = np.zeros(t.size)

freq, theta = pypdm.pdm(t, y, s, 0.01, 1, 0.001, 10)

fig, (ax1, ax2) = plt.subplots(2, 1, constrained_layout=True)

ax1.plot(t, y, 'k.')
ax1.set_xlabel('Time')
ax1.set_ylabel('Mag')

ax2.plot(freq, theta, 'k')
ax2.axvline(1/2/np.pi, color='red')
for i in range(2,4):
    ax2.axvline(1/2/i/np.pi, color='red', ls='--')
ax2.set_xlabel('Frequency')
ax2.set_ylabel('Theta')
plt.show()
# plt.savefig('Py-PDM-example.png')
