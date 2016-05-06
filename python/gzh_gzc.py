# Import any necessary packages

from matplotlib import pyplot as plt
from matplotlib.colors import LogNorm
from astropy.io import fits
import numpy as np

import warnings

warnings.filterwarnings('ignore', category=RuntimeWarning, append=True)

data = fits.getdata('/Users/willettk/Dropbox/gzhubble/gzh_gzcandels.fits',1)

pd = [
      {'gzc_col':'t00_smooth_or_featured_a1_features_frac',
       'gzh_col':'t01_smooth_or_features_a02_features_or_disk_fraction',
       'label':'features (raw)'},
      {'gzc_col':'t00_smooth_or_featured_a1_features_weighted_frac',
       'gzh_col':'t01_smooth_or_features_a02_features_or_disk_weighted_fraction',
       'label':'features (weighted)'},
      {'gzc_col':'t00_smooth_or_featured_a1_features_weighted_frac',
       'gzh_col':'t01_smooth_or_features_a02_features_or_disk_best_fraction',
       'label':'features (best)'},
      {'gzc_col':'t09_disk_edge_on_a0_yes_weighted_frac',
       'gzh_col':'t02_edgeon_a01_yes_weighted_fraction',
       'label':'edge-on',
       'constraints':{'t02_edgeon_total_weight':20,'t12_clumpy_a02_no_weighted_fraction':0.3}},
      {'gzc_col':'t12_spiral_pattern_a0_yes_weighted_frac',
       'gzh_col':'t04_spiral_a01_spiral_weighted_fraction',
       'label':'spiral',
       'constraints':{'t04_spiral_total_weight':20,'t02_edgeon_a02_no_weighted_fraction':0.25}},
      {'gzc_col':'t16_merging_tidal_debris_a0_merging_frac',
       'gzh_col':'t08_odd_feature_a06_merger_weighted_fraction',
       'label':'merger',
       'constraints':{'t01_smooth_or_features_a03_star_or_artifact_weighted_fraction':0.10,'t06_odd_total_weight':20}},
       ]


# Make a Simmons-style plot of the averaged vote fractions, akin to Figure 7 in the GZC paper.

fig,axarr = plt.subplots(2,3,figsize=(16,10))

# Colors
blue = '#377EB8'
red = '#E41A1C'
green = '#4DAF4A'

for p,ax in zip(pd,axarr.ravel()):
    if p.has_key('constraints'):
        for k in p['constraints']:
            if k == 't01_smooth_or_features_a03_star_or_artifact_weighted_fraction':
                data = data[data[k] < p['constraints'][k]]
            else:
                data = data[data[k] >= p['constraints'][k]]

    gzc = data[p['gzc_col']]
    gzh = data[p['gzh_col']]

    # Bin by the _other_ axis

    nbins = 10
    bins = np.linspace(0,1,nbins)
    db = bins[1] - bins[0]

    gzc_avg,gzc_std = np.zeros(nbins-1),np.zeros(nbins-1)
    gzh_avg,gzh_std = np.zeros(nbins-1),np.zeros(nbins-1)
    for i,b in enumerate(bins[:-2]):
        gzc_avg[i] = np.mean(gzc[(gzh >= b) & (gzh < b+db)])
        gzh_avg[i] = np.mean(gzh[(gzc >= b) & (gzc < b+db)])
        gzc_std[i] = np.std(gzc[(gzh >= b) & (gzh < b+db)])
        gzh_std[i] = np.std(gzh[(gzc >= b) & (gzc < b+db)])
    gzc_avg[-1] =  np.mean(gzc[(gzh > bins[-2]) & (gzh <= bins[-1])])
    gzh_avg[-1] =  np.mean(gzh[(gzc > bins[-2]) & (gzc <= bins[-1])])
    gzc_std[-1] =  np.std(gzc[(gzh > bins[-2]) & (gzh <= bins[-1])])
    gzh_std[-1] =  np.std(gzh[(gzc > bins[-2]) & (gzc <= bins[-1])])

    hex1 = ax.hexbin(gzc,gzh,gridsize=50,cmap=plt.cm.gray_r,norm=LogNorm())
    cb = plt.colorbar(hex1,ax=ax)
    #cb.set_label(r'$N$',fontsize=16)

    plotbins = bins[:-1] + db/2.
    ax.errorbar(gzc_avg,plotbins,xerr=gzc_std,color=blue,marker='o',ls='-',markersize=10)
    ax.errorbar(plotbins,gzh_avg,yerr=gzh_std,color=red,marker='o',ls='-',markersize=10)

    ax.set_xlim(-0.05,1.05)
    ax.set_ylim(-0.05,1.05)
    ax.set_xlabel(r'GZC {0} vote fraction'.format(p['label']),fontsize=16)
    ax.set_ylabel(r'GZH {0} vote fraction'.format(p['label']),fontsize=16)

    # Super kluge
    if ax == axarr.ravel()[2]:
        ax.set_xlabel(r'GZC features (weighted) vote fraction',fontsize=16)


fig.tight_layout()
plt.show()
