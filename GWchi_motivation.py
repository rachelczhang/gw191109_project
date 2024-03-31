import matplotlib.pyplot as plt
import pesummary
from pesummary.io import read
print(pesummary.__version__)
import h5py
import numpy as np 
from numpy import trapz
from scipy import stats

# %config InlineBackend.figure_format='retina'

file_name = './GW191109_010717.h5'
file_name1 = './GW200225_060421.h5'
data = read(file_name)
data1 = read(file_name1)
print('Found run labels:')
print(data.labels)
print(data1.labels)

samples_dict = data.samples_dict
samples_dict1 = data1.samples_dict
posterior_samples = samples_dict['C01:Mixed']
posterior_samples1 = samples_dict1['C01:Mixed']
parameters = sorted(list(posterior_samples.keys()))
parameters1 = sorted(list(posterior_samples1.keys()))
print(parameters)
print(parameters1)
print('chieff', posterior_samples['chi_eff'])
print('chieff1', posterior_samples1['chi_eff'])

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

# fig = posterior_samples.plot('mass_1_source', type='hist')
# posterior_samples.plot('mass_1_source', type='hist',fig=fig, kde=True)
# fig.show()

# fig = posterior_samples.plot('mass_2_source', type='hist')
# posterior_samples.plot('mass_2_source', type='hist',fig=fig, kde=True)
# fig.show()

# fig = posterior_samples.plot('chi_eff', type='hist')
# posterior_samples.plot('chi_eff', type='hist',fig=fig, kde=True)
# fig.show()

# fig = posterior_samples.plot('spin_1x', type='hist')
# posterior_samples.plot('spin_1x', type='hist',fig=fig, kde=True)
# fig.show()

# fig = posterior_samples.plot('spin_1y', type='hist')
# posterior_samples.plot('spin_1y', type='hist',fig=fig, kde=True)
# fig.show()

# fig = posterior_samples.plot('spin_1z', type='hist')
# posterior_samples.plot('spin_1z', type='hist',fig=fig, kde=True)
# fig.show()

# fig = posterior_samples.plot('spin_2x', type='hist')
# posterior_samples.plot('spin_2x', type='hist',fig=fig, kde=True)
# fig.show()

# fig = posterior_samples.plot('spin_2y', type='hist')
# posterior_samples.plot('spin_2y', type='hist',fig=fig, kde=True)
# fig.show()

# fig = posterior_samples.plot('spin_2z', type='hist')
# posterior_samples.plot('spin_2z', type='hist',fig=fig, kde=True)
# fig.show()

# fig = posterior_samples.plot(type='corner', parameters=['mass_1_source','mass_2_source', 'chi_eff'])
chieff = posterior_samples['chi_eff']
chieff1 = posterior_samples1['chi_eff']
kde = stats.gaussian_kde(chieff)
kde1 = stats.gaussian_kde(chieff1)
plt.hist(chieff, density=True, bins=100, alpha=0.3)
plt.hist(chieff1, density=True, bins=100, alpha=0.3)
plt.show()

x = np.linspace(-1, 0.6, 1000)
plt.plot(x, kde(x), linewidth=4, label='GW191109')
plt.axvline(x=0, color='y', linestyle='--', linewidth=4)
# plt.text(-0.36, 0.75, r'90.6\%')
# plt.text(0.1, 0.06, r'9.4\%')
plt.margins(x=0, y=0)
plt.ylim(0, 1.9)
plt.xlabel(r'$\chi_{\rm eff}$')
plt.ylabel(r'Probability density function')

x = np.linspace(-1, 0.6, 1000)
plt.plot(x, kde1(x), linewidth=4, color='black', label='GW200225')
# plt.axvline(x=0, color='r', linestyle='--', linewidth=4)
# plt.text(-0.15, 1.25, r'86.7\%')
# plt.text(0.02, 0.09, r'13.3\%')
plt.margins(x=0, y=0)
plt.legend()
plt.ylim(0, 3.25)
plt.xlabel(r'$\chi_{\rm eff}$')
plt.ylabel(r'Probability density function')
plt.show()

n, bins, patches = plt.hist(chieff1, density=True, bins=100)
plt.show()

binscent = np.array([0.5 * (bins[i] + bins[i+1]) for i in range(len(bins)-1)])
plt.plot(binscent, n)
print('binslength', len(binscent))
print('nlength', len(n))

positive_binscent = []
positive_n = []
for i in range(len(binscent)):
	if binscent[i] > 0:
		positive_binscent.append(binscent[i])
		positive_n.append(n[i])

area = trapz(n, x = binscent)
print("area = ", area)

pos_area = trapz(positive_n, x = positive_binscent)
print('positive area = ', pos_area)

# plt.axvline(x=0, color='r', linestyle='--')
# plt.text(-0.36, 0.75, r'90.6\%')
# plt.text(0.1, 0.06, r'9.4\%')
# # fig = posterior_samples.plot('chi_eff', type='hist', kde=True)
# # posterior_samples.plot('chi_eff', type='hist',fig=fig, kde=True)
# # fig.show()
# plt.margins(x=0, y=0)
# plt.ylim(0, 1.9)
# plt.xlabel(r'$\chi_{\rm eff}$')
# plt.ylabel(r'Probability density function')
# plt.show()