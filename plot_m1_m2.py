import numpy as np 
import matplotlib.pyplot as plt 
from matplotlib import rcParams
import pesummary
from pesummary.io import read
print(pesummary.__version__)
import h5py
from numpy import trapz
from scipy import stats
from matplotlib import transforms
from scipy.stats import gaussian_kde as kde
import pickle
import os


file_name = './GW191109_010717.h5'
# file_name1 = './GW200225_060421.h5'
data = read(file_name)
# data1 = read(file_name1)
# print('Found run labels:')
# print(data.labels)
# print(data1.labels)

samples_dict = data.samples_dict
# samples_dict1 = data1.samples_dict
posterior_samples = samples_dict['C01:Mixed']
# posterior_samples1 = samples_dict1['C01:Mixed']
samps = dict(GW191109=posterior_samples)
parameters = sorted(list(posterior_samples.keys()))
print('parameters', parameters)
# parameters1 = sorted(list(posterior_samples1.keys()))
mass_1 = posterior_samples['mass_1_source']
mass_2 = posterior_samples['mass_2_source']

plt.hist(mass_1)
plt.hist(mass_2)
plt.show()
# plt.hist2d(mass_1,mass_2)
# plt.show()

with open('all_catalog_bhmergers.txt') as f:
	next(f)
	lines = f.readlines()
m1_1g1g = []
m2_1g1g = []
m1_1g2g = []
m2_1g2g = []
m1_2g2g = []
m2_2g2g = []
for line in lines:
	x = line.split('\t')
	if float(x[10][0:7]) > 0.6 and float(x[10][-8:-1]) > 0.6:
		if float(x[8]) >= float(x[9]):
			m1_2g2g.append(float(x[8]))
			m2_2g2g.append(float(x[9]))
		else:
			m1_2g2g.append(float(x[9]))
			m2_2g2g.append(float(x[8]))
	elif (float(x[10][0:7]) > 0.6 and float(x[10][-8:-1]) < 0.1) or (float(x[10][-8:-1]) > 0.6 and float(x[10][0:7]) < 0.1):
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
	elif float(x[10][0:7]) < 0.1 and float(x[10][-8:-1]) < 0.1:
		if float(x[8]) >= float(x[9]):
			m1_1g1g.append(float(x[8]))
			m2_1g1g.append(float(x[9]))
		else:
			m1_1g1g.append(float(x[9]))
			m2_1g1g.append(float(x[8]))

def count(m1_list, m2_list):
	count = 0
	for i in range(len(m1_list)):
		if m1_list[i] <= 76 and m1_list[i] >= 54 and m2_list[i] <= 62 and m2_list[i] >= 32:
			count += 1
	return count
	
print('1g2g count', count(m1_1g2g, m2_1g2g))
print('2g2g count', count(m1_2g2g, m2_2g2g))

m_1g = m1_1g1g + m2_1g1g + m2_1g2g
m_2g = m1_2g2g + m2_2g2g + m1_1g2g 

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

plt.hist(m_1g, label='1g', bins=80, density=True, histtype='step')
plt.hist(m_2g, label='2g', bins=50, density=True, histtype='step')
plt.xlim(0, 100)
plt.xlabel(r'Mass [$M_\odot$]')
plt.ylabel('Relative Number')
plt.legend()
plt.show()

# # gw19_m1 = posterior_samples['mass_1_source']
# # gw20_m1 = posterior_samples1['mass_1_source']
# # kde = stats.gaussian_kde(gw19_m1)
# # kde1 = stats.gaussian_kde(gw20_m1)
# # gw19m1_range = np.linspace(50, 80, 1000)
# # gw20m1_range = np.linspace(10, 25, 1000)
# # # plt.plot(m1_range, kde(m1_range), linewidth=4, label='GW191109_010717')
# # # plt.plot(m2_range, kde1(m2_range), linewidth=4, label='GW200225_060421')
# # gw19_m2 = posterior_samples['mass_2_source']
# # gw20_m2 = posterior_samples1['mass_2_source']
# # kde_m2 = stats.gaussian_kde(gw19_m2)
# # kde1_m2 = stats.gaussian_kde(gw20_m2)
# # gw19m2_range = np.linspace(30, 70, 1000)
# # gw20m2_range = np.linspace(10, 20, 1000)
# # base = plt.gca().transData
# # rot = transforms.Affine2D().rotate_deg(90)

# # def scatter_hist(x, y, x1, y1, ax, ax_histx, ax_histy):
# #     # no labels
# # 	ax_histx.tick_params(axis="x", labelbottom=False)
# # 	ax_histy.tick_params(axis="y", labelleft=False)

# #     # the scatter plot:
# # 	ax.scatter(x, y, s=5, alpha=0.5,  label='1G-2G')
# # 	ax.scatter(x1, y1, s=5, c='r', alpha=0.5, label='2G-2G')
# #     # now determine nice limits by hand:
# # 	# binwidth = 0.25
# # 	# xymax = max(np.max(np.abs(x)), np.max(np.abs(y)))
# # 	# lim = (int(xymax/binwidth) + 1) * binwidth
# # 	# print('lim', lim)
# # 	# bins = np.arange(-lim, lim + binwidth, binwidth)
# # 	# ax_histx.hist(x, bins=bins, density=True)
# # 	# ax_histy.hist(y, bins=bins, orientation='horizontal', density=True)
# # 	ax_histx.hist(x, bins=np.logspace(np.log10(10), np.log10(200), 70), alpha=0.5, density=True)
# # 	ax_histx.hist(x1, bins=np.logspace(np.log10(10), np.log10(200), 70), color='r', alpha=0.5, density=True)
# # 	ax_histx.plot(gw19m1_range, kde(gw19m1_range), linewidth=2, color='black', alpha=0.5, label='GW191109_010717')
# # 	ax_histx.plot(gw20m1_range, kde1(gw20m1_range), linewidth=2, color='orange', alpha=0.5, label='GW200225_060421')
# # 	ax_histy.hist(y, bins=np.logspace(np.log10(10), np.log10(200), 70), alpha=0.5, orientation='horizontal', density=True)
# # 	ax_histy.hist(y1, bins=np.logspace(np.log10(10), np.log10(200), 70), color='r', alpha=0.5, orientation='horizontal', density=True)
# # 	ax_histy.plot(gw19m2_range, kde(gw19m2_range), linewidth=2, color='black', alpha=0.5, transform=base+rot, label='GW191109_010717')
# # 	ax_histy.plot(gw20m2_range, kde1(gw20m2_range), linewidth=2, color='orange', alpha=0.5, transform=base+rot, label='GW200225_060421')

# # # definitions for the axes
# # left, width = 0.1, 0.65 # moves to right, stretches image horizontally
# # bottom, height = 0.1, 0.65 # moves up, stretches image vertically
# # spacing = 0.005 # space between 2D and 1D plots

# # rect_scatter = [left, bottom, width, height]
# # rect_histx = [left, bottom + height + spacing, width, 0.2]
# # rect_histy = [left + width + spacing, bottom, 0.2, height]
# # fig = plt.figure()
# # ax = fig.add_axes(rect_scatter)
# # ax_histx = fig.add_axes(rect_histx, sharex=ax)
# # ax_histy = fig.add_axes(rect_histy, sharey=ax)

# # # use the previously defined function
# # scatter_hist(m1_1g2g, m2_1g2g, m1_2g2g, m2_2g2g, ax, ax_histx, ax_histy)
# # # ax.scatter(m1_2g2g, m2_2g2g, s=5, c='r', alpha=0.5, label='2G-2G')
# # ax.scatter(65, 47, s=70, c='black', alpha=0.4, marker=(5, 1), label='GW191109_010717')
# # ax.scatter(19.3, 14, s=70, c='orange', alpha=0.8,  marker=(5, 1), label='GW200225_060421')
# # ax.errorbar(19.3, 14, yerr=[[3], [5]], xerr=[[3.5], [2.8]], ecolor='orange', alpha = 0.6, elinewidth=3, capsize=10)
# # ax.errorbar(65, 47, yerr=[[13], [15]], xerr=11, ecolor='black', alpha = 0.4, elinewidth=3, capsize=10)
# # ax.set_xlabel(r'Primary mass [$\rm M_{\rm sol}$]')
# # ax.set_ylabel(r'Secondary mass [$\rm M_{\rm sol}$]')
# # ax.legend(loc='upper left')
# # ax.set_xscale('log')
# # ax.set_yscale('log')
# # ax.grid(visible=True, which='both')
# # plt.show()
m1_less = []
m2_less = []
m1_more = []
m2_more = []
for i in range(len(m1_1g2g)):
	if m2_1g2g[i] >= 40.5:
		m1_more.append(m1_1g2g[i])
		m2_more.append(m2_1g2g[i])
	else:
		m1_less.append(m1_1g2g[i])
		m2_less.append(m2_1g2g[i])

# print('mass1', mass_1)
# h, xedges, yedges, image = plt.hist2d(mass_1, mass_2)
# print('h', h) # [[], []] each list inside the overall list goes from bottom left to vertically up and continues to a new x for each new list
# print('xedges', xedges)
# print('yedges', yedges)
# plt.show()

def get_midpoints_edges(edges):
	mid_edges = []
	for e in range(len(edges)-1):
		mid_edges.append((edges[e]+edges[e+1])/2)
	return mid_edges

def pickle_contour_data(event, post_name, contour_dir, cdata):
    outfile = os.path.join(contour_dir, '{}_{}_contour_data.pkl'.format(event, post_name))
    with open(outfile, 'wb') as outp:
        pickle.dump(cdata, outp)
        
class Bounded_2d_kde(kde):
    r"""Represents a two-dimensional Gaussian kernel density estimator
    for a probability distribution function that exists on a bounded
    domain."""

    def __init__(self, pts, xlow=None, xhigh=None, ylow=None, yhigh=None, *args, **kwargs):
        """Initialize with the given bounds.  Either ``low`` or
        ``high`` may be ``None`` if the bounds are one-sided.  Extra
        parameters are passed to :class:`gaussian_kde`.

        :param xlow: The lower x domain boundary.

        :param xhigh: The upper x domain boundary.

        :param ylow: The lower y domain boundary.

        :param yhigh: The upper y domain boundary.
        """
        pts = np.atleast_2d(pts)

        assert pts.ndim == 2, 'Bounded_kde can only be two-dimensional'

        super(Bounded_2d_kde, self).__init__(pts.T, *args, **kwargs)

        self._xlow = xlow
        self._xhigh = xhigh
        self._ylow = ylow
        self._yhigh = yhigh

    @property
    def xlow(self):
        """The lower bound of the x domain."""
        return self._xlow

    @property
    def xhigh(self):
        """The upper bound of the x domain."""
        return self._xhigh

    @property
    def ylow(self):
        """The lower bound of the y domain."""
        return self._ylow

    @property
    def yhigh(self):
        """The upper bound of the y domain."""
        return self._yhigh

    def evaluate(self, pts):
        """Return an estimate of the density evaluated at the given
        points."""
        pts = np.atleast_2d(pts)
        assert pts.ndim == 2, 'points must be two-dimensional'

        x, y = pts.T
        pdf = super(Bounded_2d_kde, self).evaluate(pts.T)
        if self.xlow is not None:
            pdf += super(Bounded_2d_kde, self).evaluate([2*self.xlow - x, y])

        if self.xhigh is not None:
            pdf += super(Bounded_2d_kde, self).evaluate([2*self.xhigh - x, y])

        if self.ylow is not None:
            pdf += super(Bounded_2d_kde, self).evaluate([x, 2*self.ylow - y])

        if self.yhigh is not None:
            pdf += super(Bounded_2d_kde, self).evaluate([x, 2*self.yhigh - y])

        if self.xlow is not None:
            if self.ylow is not None:
                pdf += super(Bounded_2d_kde, self).evaluate([2*self.xlow - x, 2*self.ylow - y])

            if self.yhigh is not None:
                pdf += super(Bounded_2d_kde, self).evaluate([2*self.xlow - x, 2*self.yhigh - y])

        if self.xhigh is not None:
            if self.ylow is not None:
                pdf += super(Bounded_2d_kde, self).evaluate([2*self.xhigh - x, 2*self.ylow - y])
            if self.yhigh is not None:
                pdf += super(Bounded_2d_kde, self).evaluate([2*self.xhigh - x, 2*self.yhigh - y])

        return pdf

    def __call__(self, pts):
        pts = np.atleast_2d(pts)
        out_of_bounds = np.zeros(pts.shape[0], dtype='bool')

        if self.xlow is not None:
            out_of_bounds[pts[:, 0] < self.xlow] = True
        if self.xhigh is not None:
            out_of_bounds[pts[:, 0] > self.xhigh] = True
        if self.ylow is not None:
            out_of_bounds[pts[:, 1] < self.ylow] = True
        if self.yhigh is not None:
            out_of_bounds[pts[:, 1] > self.yhigh] = True

        results = self.evaluate(pts)
        results[out_of_bounds] = 0.
        return results

def ms2q(pts):
    """Transformation function from component masses to chirp mass and mass ratio"""
    pts = np.atleast_2d(pts)

    m1 = pts[:, 0]
    m2 = pts[:, 1]
    assert (m1 > 0).all() and (m2 > 0).all()
    mc = np.power(m1 * m2, 3./5.) * np.power(m1 + m2, -1./5.)
    if not np.isfinite(mc).all():
        idx_bad = np.where(np.logical_not(np.isfinite(mc)))
        print('ms2q(): infinite or nan value encountered in chirp mass')
        print('idx:', idx_bad)
        print('mc:', mc[idx_bad])
        print('m1:', m1[idx_bad])
        print('m2:', m2[idx_bad])
    q = m2/m1
    return np.column_stack([mc, q])

def estimate_2d_post(params, post, data2=None,
                     xlow=None, xhigh=None, ylow=None, yhigh=None, transform=None,
                     gridsize=500, **kwargs):
    x = post[params[0]]
    y = post[params[1]]

    if transform is None:
        transform = lambda x: x

    deltax = x.max() - x.min()
    deltay = y.max() - y.min()
    x_pts = np.linspace(x.min() - .1*deltax, x.max() + .1*deltax, gridsize)
    y_pts = np.linspace(y.min() - .1*deltay, y.max() + .1*deltay, gridsize)

    xx, yy = np.meshgrid(x_pts, y_pts)

    positions = np.column_stack([xx.ravel(), yy.ravel()])

    # Calculate the KDE
    pts = np.array([x, y]).T
    selected_indices = np.random.choice(len(pts), len(pts)//2, replace=False)
    kde_sel = np.zeros(len(pts), dtype=bool)
    kde_sel[selected_indices] = True
    kde_pts = transform(pts[kde_sel])
    untransformed_den_pts = pts[~kde_sel]
    den_pts = transform(untransformed_den_pts)

    Nden = den_pts.shape[0]

    post_kde=Bounded_2d_kde(kde_pts, xlow=xlow, xhigh=xhigh, ylow=ylow, yhigh=yhigh)
    den = post_kde(den_pts)
    inds = np.argsort(den)[::-1]
    den = den[inds]

    z = np.reshape(post_kde(transform(positions)), xx.shape)

    return {'xx':xx, 'yy':yy, 'z':z, 'kde':den, 'kde_sel':kde_sel}

def generate_contour_data(events_dict_O, contour_dir, name, gs=50):
    """Save contour data for all events to `contour_dir`.
    """
    post_name = name
    yhigh=1.0
    print('\nGenerating contour data for %s at gridsize %d ...'%(post_name, gs))
    for event_name, pos in events_dict_O.items():
        print(event_name)
        contour_data = estimate_2d_post(['mass_1_source', 'mass_2_source'], pos, 
                                        transform=ms2q, yhigh=yhigh, gridsize=gs)
        pickle_contour_data(event_name, post_name, contour_dir, contour_data)

def getCL(Z,CL):
    #Z = P(vars) with vars on equally spaced grid
    zf = Z.flatten()
    #get total probability in Z (up to grid spacing normalization factor)
    ptot = np.sum(zf)

    # Get indices of sorted array of Zs
    i_sort = np.argsort(zf)[::-1]
    levels=[]
    # Iterate until you've accumulated 90% of the probability
    for cl in CL:
        zsum = 0
        zlast = 0
        j=0
        while zsum/ptot<=cl:
            i_max = i_sort[j]
            j+=1
            zlast=zf[i_max]
            zsum+=zf[i_max]
        levels.append(zlast)
    return levels  

generate_contour_data(samps, '.','with_reflection',gs=100)
with open('GW191109_with_reflection_contour_data.pkl','rb') as f:
    x = pickle.load(f)
levels = getCL(x['z'],CL=[.9])
plt.contour(x['xx'],x['yy'],x['z'],levels=levels, colors='black')

# # this worked without accounting for the reflection boundary 
# kde_post = stats.gaussian_kde([mass_1, mass_2])
# m1 = m2 = np.linspace(0, 150, 100)
# grid_m1, grid_m2 = np.meshgrid(m1,m2)
# Zpost = np.reshape(kde_post(np.vstack([grid_m1.ravel(),grid_m2.ravel()])).T, grid_m1.shape)
# plt.contour(grid_m1, grid_m2, Zpost, levels=getCL(Zpost, [0.9, 0.5]), colors='black')

# xmidedges = get_midpoints_edges(xedges)
# ymidedges = get_midpoints_edges(yedges)
# h_transp = np.array(h.T.tolist())
# print('h_transp', h_transp)
# plt.scatter(m1_1g1g, m2_1g1g, s=5, c='orange', alpha=0.5,  label='1G-1G')
# plt.scatter(m1_1g2g, m2_1g2g, s=5, c='b', alpha=0.5,  label='1G-2G')
plt.scatter(m1_less, m2_less, s=10, c='b', alpha=1,  label='1G-2G')
plt.scatter(m1_more, m2_more, s=10, c='orange', alpha=1,  label='1G(star collisions)-2G')
plt.scatter(m1_2g2g, m2_2g2g, s=10, c='r', alpha=1, label='2G-2G')
plt.scatter(65, 47, s=70, c='black', alpha=1, marker=(5, 1), label='GW191109_010717')
# plt.hist2d(mass_1,mass_2, alpha=0.5, cmap='Greys')
# plt.scatter(19.3, 14, s=70, c='cyan', alpha=0.4,  marker=(5, 1), label='GW200225_060421')
#plt.axvspan(xmin = 54, xmax = 76, ymin=34, ymax=62)
# plt.errorbar(19.3, 14, yerr=[[3], [5]], xerr=[[3.5], [2.8]], ecolor='cyan', alpha = 0.4, elinewidth=3, capsize=8)
# plt.errorbar(65, 47, yerr=[[13], [15]], xerr=11, ecolor='black', alpha = 0.4, elinewidth=3, capsize=8)
plt.xlabel(r'Primary mass [$\rm M_{\odot}$]')
plt.ylabel(r'Secondary mass [$\rm M_{\odot}$]')
# print('xmidedges', xmidedges)
# print('ymidedges', ymidedges)
# plt.contour(xmidedges, ymidedges, h_transp, levels=np.linspace(np.min(h), np.max(h), 9), cmap='Greys')
plt.legend()
plt.xscale('log')
plt.yscale('log')
plt.grid(visible=True, which='both')
# plt.ylim(5, 150)
plt.show()
plt.savefig('m1m2.png')

