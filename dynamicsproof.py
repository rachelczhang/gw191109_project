import matplotlib.pyplot as plt 
import numpy as np 

def chieff(m1, m2, chi1, chi2):
	# print('chi1', chi1)
	# print('chi2', chi2)
	chi = []
	for i in range(len(m1)):
		chi.append((m1[i]*chi1[i]+m2[i]*chi2[i])/(m1[i]+m2[i]))
	return chi 

n = 10**5
costheta1 = np.random.uniform(-1, 1, size=n)
costheta2 = np.random.uniform(-1, 1, size=n)

# case 1: m1=47, m2=65, chi1=0, chi2=0.7
m1 = [47]*n
m2 = [65]*n
chi1_1g = [0]*n
chi1_2g = [0.7*i for i in costheta1]
chi2_2g = [0.7*i for i in costheta2]
# x = [chi1[i]+chi2[i] for i in range(n)]
chi12 = chieff(m1, m2, chi1_1g, chi2_2g)
chi22 = chieff(m1, m2, chi1_2g, chi2_2g)
# print('chieff', chi)

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
plt.rcParams['text.latex.preamble'] = [r'\usepackage{amsmath}'] #for \text command
plt.rcParams.update({'font.size': 22})

plt.hist(chi12, density=True, histtype = 'step', linewidth=3, label='1G-2G')
plt.hist(chi22, density=True, histtype='step', linewidth=3, color='r', label='2G-2G')
plt.xlim(-0.7, 0.7)
plt.xlabel(r'$\rm \chi_{eff}$')
plt.ylabel('Relative number of binary BHs')
plt.legend()
plt.grid(visible=True, which='both')
plt.show()

# plt.scatter(x, chi)
# plt.show()