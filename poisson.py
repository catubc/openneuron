import numpy as np
import matplotlib.pyplot as plt

n_spikes = 200
spikes = np.random.poisson(10, n_spikes)+(np.random.randint(200))
print spikes

bin_width = 2  #100ms bins
y = np.histogram(spikes, bins = np.arange(0,200,bin_width))
plt.plot(y[1][:-1], y[0], linewidth=3, color='blue')
plt.xlim(0,250)
plt.ylim(0,100)
plt.show()
