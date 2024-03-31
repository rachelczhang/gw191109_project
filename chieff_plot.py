import pandas as pd 
import matplotlib.pyplot as plt 
from scipy.integrate import simps

filename = pd.read_csv('chieff.csv', header=None, delimiter=',')
print('filename', filename)

x = filename[0]
y = filename[1]

plt.plot(x, y)
plt.margins(x=0, y=0)
plt.ylim(0, 1.8)
plt.show()
