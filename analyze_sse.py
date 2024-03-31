import numpy as np 
import matplotlib.pyplot as plt 
import pandas as pd 
import csv 

bcm_data = pd.read_hdf('bcm_0.001.h5', 'bcm')
bpp_data = pd.read_hdf('bpp_0.001.h5', 'bpp')

############# VERSION 2 ##################
preSN_list = []
diff_mass = []
# save list of all binaries that go through SN
for i in np.unique(bpp_data.index.values):
	# print(bpp_data.loc[3][['tphys', 'mass_1', 'evol_type']])
	# print(bpp_data.loc[3][['tphys', 'mass_1', 'evol_type']].iloc[0])
	# print('len', len(bpp_data.loc[3].index))
	if bpp_data.loc[i][bpp_data.loc[i]['evol_type'] == 15.0].size != 0:
		# check if post-SN mass is between 34 and 62 solar masses
		postSN_mass = bpp_data.loc[i][(bpp_data.loc[i]['evol_type'] == 15.0)]['mass_1']
		# print('post', postSN_mass)
		# print('post2', postSN_mass.item())
		if postSN_mass.item() <= 62 and postSN_mass.item() >= 34:
			# find the index at which the supernova happens
			for j in range(len(bpp_data.loc[i])):
				if bpp_data.loc[i]['evol_type'].iloc[j] == 15.0:
					break
			# index before supernova gives pre-SN mass
			preSN_mass = bpp_data.loc[i].iloc[j-1]['mass_1']
			preSN_list.append(preSN_mass)
			diff = preSN_mass.item()-postSN_mass.item()
			diff_mass.append(diff)
			# print('preSN mass', preSN_mass)
			# print('postSN mass', postSN_mass)
params = {
        # latex
        'text.usetex': True,
        # fonts
        'font.family': 'serif', 

        # figure and axes
        'figure.figsize': (14,7),
        'figure.titlesize': 35, 
        'axes.grid': False, 
        'axes.titlesize':25,
        #'axes.labelweight': 'bold', 
        'axes.labelsize': 25,

        # tick markers
        'xtick.direction': 'in', 
        'ytick.direction': 'in', 
        'xtick.labelsize': 22,
        'ytick.labelsize': 22, 
        'xtick.major.size': 10.0,
        'ytick.major.size': 10.0,
        'xtick.minor.size': 3.0,
        'ytick.minor.size': 3.0,

        # legend
        'legend.fontsize': 20,
        'legend.frameon': True,
        #'legend.framealpha':1.0,

        # colors
        'image.cmap': 'viridis',

        # saving figures
        'savefig.dpi': 300
        }

plt.rcParams.update(params)
plt.rcParams['font.serif'] = ['Computer Modern', 'Times New Roman']
plt.rcParams['font.family'] = ['serif', 'STIXGeneral']
plt.rcParams['savefig.bbox'] = 'tight'
plt.rcParams['text.latex.preamble'] = r'\usepackage{amsmath}' #for \text command
plt.rcParams.update({'font.size': 22})
plt.grid(visible=True, which='both')
#  for A1
# plt.xlim(0, 1.8)
# plt.hist(diff_mass, density=True, bins=500)
#  for A2
print('diff mass', diff_mass)
plt.hist(diff_mass, density=True, bins=15)
plt.xlabel(r'Mass loss [M$_{\odot}$]')
plt.ylabel(r'Relative number of stars')
plt.show()

############# VERSION 1 ###################
range47 = bcm_data[(bcm_data['tphys'] == 13700.0) & (bcm_data['mass0_1'] >= 34) & (bcm_data['mass0_1'] <= 62)]
# print('range47', range47)

mold_list = []
mnew_list = []
diff_mass = []
for i in np.unique(range47.index.values):
	timestamp = bpp_data.loc[i][(bpp_data.loc[i]['evol_type'] == 15.0)]['tphys']
	if timestamp.size != 0:
		sn_data = bpp_data.loc[i][(bpp_data.loc[i]['tphys'] == timestamp.iloc[0])]
		# print('sn data', sn_data)
		m_old = sn_data.iloc[0]['mass_1']
		m_new = sn_data.iloc[1]['mass_1']
		mold_list.append(m_old)
		mnew_list.append(m_new)
		diff_mass.append(m_old-m_new)
		# if diff_mass[-1]>40:
		# 	print('full', bpp_data.loc[i])
		# 	print('sn data', sn_data)
		# 	print('mold', sn_data.iloc[0])
		# 	print('mnew', sn_data.iloc[1])
		# 	print('old evoltype', sn_data.iloc[0]['evol_type'])
		# 	print('new evoltype', sn_data.iloc[1]['evol_type'])
		# 	print('diff', diff_mass)
# print('mold dist', mold_list)
# print('mnew dist', mnew_list)
# print('diff mass', diff_mass)
plt.hist(diff_mass, density=True)
plt.xlabel(r'Mass loss [M$_{\rm sol}$]')
plt.ylabel(r'Relative number of stars')
# plt.show()

plt.hist(mold_list, density=True)
plt.xlabel(r'Pre-SN Mass [M$_{\odot}$]')
plt.ylabel(r'Relative number of stars')
# plt.show()

plt.hist(mnew_list, density=True)
plt.xlabel(r'Post-SN Mass [M$_{\odot}$]')
plt.ylabel(r'Relative number of stars')
# plt.show()