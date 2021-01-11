import matplotlib.pyplot as plt
import numpy as np
import pylab
t, x, y = np.loadtxt('generatedOutputTest1.txt', delimiter=' ',unpack=True)

pylab.plot(t, x, '-b', label='tempMin')
pylab.plot(t, y, '-r', label='tempMax')
plt.xlabel('iteration')
plt.ylabel('tempMin/tempMax[s]')
plt.title('Wykres temperatur')
plt.legend()
plt.show()