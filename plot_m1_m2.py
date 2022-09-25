import numpy as np 
import matplotlib.pyplot as plt 
from matplotlib import rcParams
#'rcParams['text.usetex'] = True

with open('all_catalog_bhmergers.txt') as f:
	next(f)
	lines = f.readlines()
m1_1g2g = []
m2_1g2g = []
m1_2g2g = []
m2_2g2g = []
for line in lines:
	x = line.split('\t')
	if float(x[10][0:7]) > 0.01 and float(x[10][-8:-1]) > 0.01:
		if float(x[8]) >= float(x[9]):
			m1_2g2g.append(float(x[8]))
			m2_2g2g.append(float(x[9]))
		else:
			m1_2g2g.append(float(x[9]))
			m2_2g2g.append(float(x[8]))
	elif (float(x[10][0:7]) > 0.01 and float(x[10][-8:-1]) < 0.01) or (float(x[10][-8:-1]) > 0.01 and float(x[10][0:7]) < 0.01):
		#if float(x[10][0:7]) > 0.01 and float(x[10][-8:-1])<0.01:
			#print('1g', (x[10][-8:-1]))
		#else:
			#print('1g', x[10][0:7])
		if float(x[8]) >= float(x[9]):
			m1_1g2g.append(float(x[8]))
			m2_1g2g.append(float(x[9]))
		else:
			m1_1g2g.append(float(x[9]))
			m2_1g2g.append(float(x[8]))
print('1g2g count', len(m1_1g2g))
print('2g2g count', len(m1_2g2g))
plt.scatter(m1_1g2g, m2_1g2g, s=5, c='b', alpha=0.5,  label='1G-2G')
plt.scatter(m1_2g2g, m2_2g2g, s=5, c='r', alpha=0.5, label='2G-2G')
plt.scatter(65, 47, s=5, c='black', alpha=0.4, marker=(5, 1), label='GW191109_010717')
plt.scatter(19.3, 14, s=5, c='cyan', alpha=0.4,  marker=(5, 1), label='GW200225_060421')
#plt.axvspan(xmin = 54, xmax = 76, ymin=34, ymax=62)
plt.errorbar(19.3, 14, yerr=[[3], [5]], xerr=[[3.5], [2.8]], ecolor='cyan', alpha = 0.4, elinewidth=3, capsize=8)
plt.errorbar(65, 47, yerr=[[13], [15]], xerr=11, ecolor='black', alpha = 0.4, elinewidth=3, capsize=8)
plt.xlabel(r'Primary mass [$M_{\rm sol}$]')
plt.ylabel(r'Secondary mass [$M_{\rm sol$}]')
plt.legend()
plt.xscale('log')
plt.yscale('log')
plt.grid(visible=True, which='both')
plt.savefig('m1m2.png')

