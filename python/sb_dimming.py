"""
Is the surface brightness of galaxies in the FERENGI images dimming as expected (1+z)^4?

References:

    http://www.astro.caltech.edu/~george/ay127/Ay127_CosmoTests.pdf
    Calvi et al. (2014) http://arxiv.org/pdf/1410.2281v1.pdf

"""

from astropy.io import fits
from matplotlib import pyplot as plt
import numpy as np

path = '/Users/willettk/Dropbox/gzhubble'

data = fits.getdata('%s/ferengi_all_weighted_and_meta.fits' % path,1)

objids = set(data['objid'])

fig = plt.figure()
ax = fig.add_subplot(111)

def plot_gals(ax,objid,label=None):
    gals = (data['objid'] == objid) & (data['sim_evolution'] == 0.)
    if gals.sum() == 8:
        # Normalize everything to start at the same surface brightess
        mu0 = data[gals]['mu_max_i'][0]
        ax.plot(data[gals]['sim_redshift'],data[gals]['mu_max_i'] - mu0,
            color='black',
            label = label,
            alpha=min(1.,0.3 + (22. - mu0)/4. * 0.7))


for objid in objids:
    plot_gals(ax,objid)

plot_gals(ax,objids.pop(),label='FERENGI')



# Overplot the expected level of cosmological dimming from Tolman

def cosmo_dim(mu_0,z,z0=0.):

    return mu_0 + 10. * np.log10(1. + (z - z0))


zarr = np.linspace(0.3,1.,100)

ax.plot(zarr,cosmo_dim(0.0,zarr,z0=0.3),color='red',lw=3,label=r'$(1+z)^4$')

ax.set_xlabel('redshift',fontsize=16)
ax.set_ylabel(r'$\Delta\mu_{max,i}$',fontsize=20)

ax.legend(loc='best')

ax.set_xlim(0,1)

plt.show()
