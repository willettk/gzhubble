from astropy.io import fits
from matplotlib import pyplot as plt

data = fits.getdata('/Users/willettk/Astronomy/Research/GalaxyZoo/fits/gz_hst_table.fits',1)

sample = data['sample']

original = (sample == 'original_initial') | (sample == 'original_flagmissed')
cosmos = (sample == 'cosmos_initial') | (sample == 'cosmos_flagmissed')
fake_agn = (sample == 'fake_agn_v1') | (sample == 'fake_agn_v1')
stripe82 =  (sample == 'stripe82') | (sample == 'stripe82_coadd')

fig,ax = plt.subplots(1,1,figsize=(8,8))
ax.hist(data[original]['t01_smooth_or_features_total_count'],bins=range(170),histtype='step',lw=2,label='AEGIS,GEMS,GOODS')
ax.hist(data[cosmos]['t01_smooth_or_features_total_count'],bins=range(170),histtype='step',lw=2,label='COSMOS')
ax.hist(data[fake_agn]['t01_smooth_or_features_total_count'],bins=range(170),histtype='step',lw=2,label='simulated AGN')
ax.hist(data[stripe82]['t01_smooth_or_features_total_count'],bins=range(170),histtype='step',lw=2,label='SDSS Stripe 82')

ax.set_xlabel('classification count',fontsize=16)
ax.set_ylabel('count',fontsize=16)

import numpy as np
print 'original',np.median(data[original]['t01_smooth_or_features_total_count'])
print 'cosmos',np.median(data[cosmos]['t01_smooth_or_features_total_count'])
print 'fake_agn',np.median(data[fake_agn]['t01_smooth_or_features_total_count'])
print 'stripe82',np.median(data[stripe82]['t01_smooth_or_features_total_count'])

ax.legend(fontsize=16)

plt.show()
