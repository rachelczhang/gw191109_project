import numpy as np 
import matplotlib.pyplot as plt 
import pandas as pd 
import csv 

bcm_data = pd.read_hdf('bcm.h5', 'bcm')
bpp_data = pd.read_hdf('bpp.h5', 'bpp')

range47 = bcm_data[(bcm_data['tphys'] == 13700.0) & (bcm_data['mass0_1'] >= 34) & (bcm_data['mass0_1'] <= 62)]
print('range47', range47)

og_masses = []
for i in np.unique(range47.index.values):
	og_masses.append(bcm_data.loc[i][(bcm_data.loc[i]['tphys'] == 0.0)]['mass0_1'])

print('og masses', og_masses)
