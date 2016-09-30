from matplotlib import pyplot as plt
from astropy.table import Table
from scipy import optimize
from scipy.stats import distributions as dist
import numpy as np

# Load data

data = Table.read("../data/ferengi_debiasable_data.fits")

# Use only galaxies with surface brightness/redshift ranges that are considered "debiasable"


data = data#[data['Correctable_Category']=='correctable']

# Limit to galaxies that have data at z_sim = 0.3, since that's what we're normalizing to.
unique_galaxies = set(data['sdss_id'])
z0ind = np.zeros(len(data),dtype=bool)
for ug in unique_galaxies:
    ind = (data['sdss_id'] == ug)
    if data[ind]['sim_redshift'].min() < 0.301:
        z0ind[ind] = True

data_z0 = data[z0ind]

def fzeta_exp(p,x):

    #y = p[0] * np.exp(-1 * (x-p[1])/p[2])

    y = np.exp(-1 * (x-0.3)/p[0])

    return y

def fzeta_lin(p,x):

    y = p[0] + p[1] * x

    return y

def fzeta(p,x):
    # results are qualitatively the same for both lin and exp versions
    return fzeta_exp(p,x)

def errfunc(p,x,y,s):

    err = (y - fzeta(p,x))/s

    return err

def errfunc_lin(p,x,y,s):

    err = (y - fzeta_lin(p,x))/s

    return err

def error_bars(k,n=40,c=0.683):

    f_gal_lower = dist.beta.ppf((1-c)/2.,k+1,n-k+1)
    f_gal_upper = dist.beta.ppf(1-(1-c)/2.,k+1,n-k+1)
    f_gal_err = (f_gal_upper - f_gal_lower) / 2.0

    return f_gal_err

def common_labels(fig,xlabel=None,ylabel=None,xfontsize=16,yfontsize=40,
                  xlabelpad=None, ylabelpad=None):

    # Set common labels
    cax = fig.add_subplot(111)    # The big subplot
    cax.set_axis_bgcolor('none')
    cax.spines['top'].set_color('none')
    cax.spines['bottom'].set_color('none')
    cax.spines['left'].set_color('none')
    cax.spines['right'].set_color('none')
    cax.tick_params(labelcolor='none', top='off', bottom='off', left='off', right='off')
    cax.set_xlabel(xlabel,fontsize=xfontsize, labelpad=xlabelpad)
    cax.set_ylabel(ylabel,fontsize=yfontsize, labelpad=ylabelpad)

    return cax

p_guess = np.array([0.5])

nrows = 4
ncols = 5

# Set up plot
fig,axarr = plt.subplots(nrows=nrows,ncols=ncols,sharex=True,sharey=True,figsize=(18,14))
bigax = common_labels(fig,'Redshift',r'$\frac{1-f_{{\rm features},z=0.3}}{1-f_{\rm features}}$',20,28, 12, 12)
zarr = np.linspace(0,1,50)

# For examples, only plot galaxies with an evolution correction of zero.
evol = 0.0
e0 = data_z0[np.absolute(data_z0['sim_evolution'] - evol) < 0.001]

e0_z0 = e0[e0['sim_redshift'] < 0.35]
#e0_z0 = data_z0

unique_galaxies = []
plist = np.linspace(0.1, 1.0, nrows*ncols+1)
nchoose = 2

for p1, p2 in zip(plist[:-1], plist[1:]):
    p_match = (e0_z0['p_features'] > p1) & (e0_z0['p_features'] <= p2)
    if p_match.sum() > nchoose:
        #print(p_match.sum(), p_match.nonzero()[0], nchoose)
        p_match = np.random.choice(p_match.nonzero()[0], nchoose, replace=False)
        nchoose = 1
    elif p_match.sum() <= nchoose:
        #print(p_match.sum(), p_match.nonzero()[0], nchoose)
        nchoose = 1 + nchoose - p_match.sum()
        p_match = p_match.nonzero()[0]
    p_match = p_match[np.argsort(e0_z0['p_features'][p_match])]
    unique_galaxies.extend(e0_z0['sdss_id'][p_match])

for ax in axarr.ravel():
    if len(unique_galaxies) == 0: break

    slen = 0
    # Make sure there are enough points to fit a function
    while slen < (len(p_guess)+1):
        ind = (e0['sdss_id'] == unique_galaxies.pop())
        slen = sum(ind)

    galaxy1 = e0[ind]
    galaxy1.sort('sim_redshift')

    z_gal = galaxy1['sim_redshift']
    f_gal = galaxy1['p_features']

    #  ADD ERROR BARS
    n = 40  # assume 40 classifications per galaxy; it'd be better to use true value, though
    f_gal_err = error_bars(f_gal*n,n)
    f_gal_norm = (1-f_gal[0]) /(1- f_gal)
    f_gal_norm_err = np.sqrt((f_gal_err/f_gal)**2 + (f_gal_err[0]/f_gal[0])**2) * f_gal_norm

    # Values must be explicitly cast as double-type precision for optimization to work. Incredibly frustrating.
    # Fix: http://stackoverflow.com/questions/12473406/scipy-optimize-leastsq-returns-best-guess-parameters-not-new-best-fit

    p, cov, infodict, mesg, ier = optimize.leastsq(errfunc,p_guess,args=(z_gal.astype(np.float64),
                                                                         f_gal_norm.astype(np.float64),
                                                                         f_gal_norm_err.astype(np.float64)),
                                                   full_output=1)
    ax.plot(z_gal,f_gal_norm,lw=2)
    ax.errorbar(z_gal,f_gal_norm, f_gal_norm_err)
    ax.plot(zarr,fzeta_exp(p,zarr),'--',lw=1)
    zeta = '={:.2f}'.format(p[0]) if p[0] <= 10.0 else '>10'
#    ax.set_title('$f_{z=0}={:.2f}\; \zeta{:s}$'.format(f_gal[0], zeta), y=1.01, fontsize=16)
    ax.set_title('$f_{features,z=0.3}=%s\; \zeta%s$'%(round(f_gal[0],2),zeta),y=1.01,fontsize=16)
    ax.set_xlim(0.11,1.05)
    ax.set_ylim(0,2)

fig.subplots_adjust(hspace=0.20, wspace=0.05)

fig.savefig('../writeup/figures/zeta_examples_sorted.pdf')
