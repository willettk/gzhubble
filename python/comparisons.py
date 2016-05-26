from PIL import Image
from astropy.io import fits
from astropy import units as u
from matplotlib import cm
from matplotlib import pyplot as plt
from matplotlib.colors import LogNorm
from numpy import ma
from scipy.ndimage.filters import gaussian_filter
import matplotlib.gridspec as gridspec
import numpy as np
import random
import scipy.stats.distributions as dist
import warnings
from astropy.cosmology import WMAP9
    

warnings.filterwarnings('ignore', category=RuntimeWarning, append=True)

gzh_path = '/Users/willettk/Astronomy/Research/GalaxyZoo/gzhubble'

# Abbreviate some of the GZH morphology categories
pfeat = 't01_smooth_or_features_a02_features_or_disk_best_fraction'
psmooth = 't01_smooth_or_features_a01_smooth_best_fraction'
pclumpy = 't12_clumpy_a01_yes_weighted_fraction'
podd = 't06_odd_a01_yes_weighted_fraction'
pstar = 't01_smooth_or_features_a03_star_or_artifact_weighted_fraction'
pspiral = 't04_spiral_a01_spiral_weighted_fraction'
pirr = 't08_odd_feature_a04_irregular_fraction'
pmerger = 't08_odd_feature_a06_merger_fraction'
pdustlane = 't08_odd_feature_a07_dust_lane_fraction'
plens = 't08_odd_feature_a02_lens_or_arc_fraction'

# Colors
blue = '#377EB8'
red = '#E41A1C'
green = '#4DAF4A'

def get_image(survey_id,goods=False,ns='n'):

    if goods:
        filename = "/Volumes/3TB/gz4/GOODS_full/jpg/standard/goods_full_{0}_g{1}_standard.jpg".format(ns,survey_id)
    else:
        filename = "/Volumes/3TB/gzh/jpg/{0}.jpg".format(survey_id)
    img = Image.open(filename)
    return img

def get_dlist(aegis,gems,goods,cassata,tasca,zest):

    dlist = [
             {'data':aegis,
              'survey':'AEGIS',
              'citation':'Lotz+08',
              'magcol':'MAG_BEST_HI',
              'maglim':23.5,
              'classcol':'TYPE'},
             {'data':gems,
              'survey':'GEMS',
              'citation':'Haeussler+07',
              'magcol':'MAG_BEST_HI',
              'maglim':23.5,
              'classcol':'TYPE'},
             {'data':goods,
              'survey':'GOODS',
              'citation':'Bundy+05',
              'magcol':'ZMAG',
              'maglim':22.5,
              'classcol':'TYPE'},
             {'data':cassata,
              'survey':'COSMOS',
              'citation':'Cassata+07',
              'magcol':'MAG_BEST_HI',
              'maglim':23.5,
              'classcol':'AUTOCLASS'},
             {'data':tasca,
              'survey':'COSMOS',
              'citation':'Tasca+11',
              'magcol':'MAG_BEST_HI',
              'maglim':23.5,
              'classcol':'CLASS_INT'},
             {'data':zest,
              'survey':'COSMOS',
              'citation':'Scarlata+07',
              'magcol':'MAG_BEST_HI',
              'maglim':23.5,
              'classcol':'TYPE'},
              ]

    return dlist

def cameron_err(c,k,n):
    f_lower = dist.beta.ppf((1-c)/2., k+1,n-k+1)
    f_upper = dist.beta.ppf(1-(1-c)/2., k+1,n-k+1)

    return f_lower,f_upper

def vl_plot():

    # Plot the volume-limited samples split by survey

    tasca = fits.getdata('{0}/comparisons/COSMOS/cosmos_tasca_gzh.fits'.format(gzh_path),1)
    cassata = fits.getdata('{0}/comparisons/COSMOS/cosmos_cassata_gzh.fits'.format(gzh_path),1)
    zest = fits.getdata('{0}/comparisons/COSMOS/cosmos_zest_gzh.fits'.format(gzh_path),1)
    aegis = fits.getdata('{0}/comparisons/AEGIS/aegis_gzh.fits'.format(gzh_path),1)
    gems = fits.getdata('{0}/comparisons/GEMS/gems_gzh.fits'.format(gzh_path),1)
    goods = fits.getdata('{0}/comparisons/GOODS/goods_gzh.fits'.format(gzh_path),1)


    dlist = get_dlist(aegis,gems,goods,cassata,tasca,zest)

    fig,axarr = plt.subplots(3,2,sharex=True,sharey=True,figsize=(12,12))

    for s,ax in zip(dlist,axarr.ravel()):
        d = s['data']
        newdata = volume_limit(d,1.0,s['magcol'],s['maglim'])
        olddata = d[d['Z_BEST'] > 0.]
        absmag_old = olddata[s['magcol']] - WMAP9.distmod(olddata['Z_BEST']).value
        absmag_new = newdata[s['magcol']] - WMAP9.distmod(newdata['Z_BEST']).value
        appmag_old = olddata[s['magcol']]
        appmag_new = newdata[s['magcol']]

        s['data'] = newdata
        ax.scatter(olddata['Z_BEST'],absmag_old,c='gray')
        ax.scatter(newdata['Z_BEST'],absmag_new,c='red')
        ax.set_xlim(-0.5,4)
        ax.set_ylim(-15,-25)
        ax.set_xlabel(r'$z$',fontsize=20)
        ax.set_ylabel(r'$M$',fontsize=20)
        ax.set_title("{0} ".format(s['survey'])+r'$m_\mathrm{I|i|z}<$'+'{0:.1f}'.format(s['maglim']),fontsize=20)


    fig.tight_layout()

    plt.show()

    plt.close()

    return None

def gems():

    do_all = False

    data = fits.getdata('{0}/comparisons/GEMS/gems_gzh.fits'.format(gzh_path),1)

    # Only consider data that didn't run up against the GALFIT constraints?

    discs = (data['n'] < 2.5)
    bulges   = (data['n'] > 2.5)

    constrained = (data['flag_constraint'] == 1)
    unconstrained = (data['flag_constraint'] == 0)

    if do_all:
        fig = plt.figure(figsize=(10,6))

        ax1 = fig.add_subplot(121)
        ax1.scatter(data[constrained & discs]['n'],
            data[constrained & discs][pfeat],
            marker='.',
            color='#377EB8',
            label=r'$n<2.5$')
        ax1.scatter(data[constrained & bulges]['n'],
            data[constrained & bulges][pfeat],
            marker='.',
            color='#E41A1C',
            label=r'$n<2.5$')
        ax1.set_title('Good GALFIT parameters ({0})'.format(constrained.sum()))

        ax2 = fig.add_subplot(122)
        ax2.scatter(data[unconstrained & discs]['n'],
            data[unconstrained & discs][pfeat],
            marker='.',
            color='#377EB8',
            label=r'$n<2.5$')
        ax2.scatter(data[unconstrained & bulges]['n'],
            data[unconstrained & bulges][pfeat],
            marker='.',
            color='#E41A1C',
            label=r'$n<2.5$')
        ax2.set_title('Bad GALFIT parameters ({0})'.format(unconstrained.sum()))

        for ax in (ax1,ax2):
            ax.set_xlabel(r'$n$',fontsize=25)
            ax.set_ylabel(r'$f_{features,best}$',fontsize=20)
            ax.set_xlim(-0.5,8.5)
            ax.set_ylim(-0.1,1.1)

        fig.tight_layout()

    gzh_feat = (data[pfeat] >= 0.8)
    gzh_int = (data[pfeat] < 0.8) & (data[pfeat] > 0.2)
    gzh_smooth = (data[pfeat] <= 0.2)

    if do_all:
        print "\n All galaxies"
        print "-------------------------------------------"

        print "\n Average f_features,best for discs: {0:.2f}+-{1:.2f}".format(np.mean(data[constrained & discs][pfeat]),np.std(data[constrained & discs][pfeat]))
        print "\n Average f_features,best for bulges: {0:.2f}+-{1:.2f}".format(np.mean(data[constrained & bulges][pfeat]),np.std(data[constrained & bulges][pfeat]))

        print "\n      f_feat >= 0.8: {0:.4f} discs (n<2.5), {1:.4f} bulges (n>2.5)".format(sum(constrained & gzh_feat & discs)*1./sum(constrained & gzh_feat),sum(constrained & gzh_feat & bulges)*1./sum(constrained & gzh_feat))
        print "\n 0.2 < f_feat < 0.8: {0:.4f} discs (n<2.5), {1:.4f} bulges (n>2.5)".format(sum(constrained & gzh_int & discs)*1./sum(constrained & gzh_int),sum(constrained & gzh_int & bulges)*1./sum(constrained & gzh_int))
        print "\n      f_feat <= 0.2: {0:.4f} discs (n<2.5), {1:.4f} bulges (n>2.5)".format(sum(constrained & gzh_smooth & discs)*1./sum(constrained & gzh_smooth),sum(constrained & gzh_smooth & bulges)*1./sum(constrained & gzh_smooth))

    # Maglim

    lim850 = 21
    maglim = (data['MAG_BEST_HI'] < lim850)

    if do_all:
        print "\n Brightness-limited (m < {0})".format(lim850)
        print "-------------------------------------------"

        print "\n Average f_features,best for discs: {0:.3f}".format(np.mean(data[constrained & maglim & discs][pfeat]))
        print "\n Average f_features,best for bulges: {0:.3f}".format(np.mean(data[constrained & maglim & bulges][pfeat]))

        print "\n      f_feat >= 0.8: {0:.4f} discs (n<2.5), {1:.4f} bulges (n>2.5)".format(sum(constrained & maglim & gzh_feat & discs)*1./sum(constrained & maglim & gzh_feat),sum(constrained & maglim & gzh_feat & bulges)*1./sum(constrained & maglim & gzh_feat))
        print "\n 0.2 < f_feat < 0.8: {0:.4f} discs (n<2.5), {1:.4f} bulges (n>2.5)".format(sum(constrained & maglim & gzh_int & discs)*1./sum(constrained & maglim & gzh_int),sum(constrained & maglim & gzh_int & bulges)*1./sum(constrained & maglim & gzh_int))
        print "\n      f_feat <= 0.2: {0:.4f} discs (n<2.5), {1:.4f} bulges (n>2.5)".format(sum(constrained & maglim & gzh_smooth & discs)*1./sum(constrained & maglim & gzh_smooth),sum(constrained & maglim & gzh_smooth & bulges)*1./sum(constrained & maglim & gzh_smooth))


    #plt.show()

    # Consider bright galaxies with bulge-like Sersic indices, but high f_features (roughly 23%). What are these galaxies,
    # and why don't the two morphological methods agree?

    subsample = (constrained & maglim & bulges & gzh_feat)
    main = (constrained & maglim)
    main_bulge = (constrained & maglim & bulges)
    main_feat = (constrained & maglim & gzh_feat)

    if do_all:
        print "\n Exploring GALFIT bulges but GZH features ({0} galaxies)".format(sum(subsample))
        print "*************************************************************"
        print "*"
        print "* (V-I) color: {0:.2f}+-{1:.2f}, {2:.2f}+-{3:.2f}, {4:.2f}+-{5:.2f}, {6:.2f}+-{7:.2f}, ".format(np.mean(data[subsample]['MAG_BEST_LOW'] - data[subsample]['MAG_BEST_HI']),np.std(data[subsample]['MAG_BEST_LOW'] - data[subsample]['MAG_BEST_HI']),np.mean(data[main]['MAG_BEST_LOW'] - data[main]['MAG_BEST_HI']),np.std(data[main]['MAG_BEST_LOW'] - data[main]['MAG_BEST_HI']),np.mean(data[main_bulge]['MAG_BEST_LOW'] - data[main_bulge]['MAG_BEST_HI']),np.std(data[main_bulge]['MAG_BEST_LOW'] - data[main_bulge]['MAG_BEST_HI']),np.mean(data[main_feat]['MAG_BEST_LOW'] - data[main_feat]['MAG_BEST_HI']),np.std(data[main_feat]['MAG_BEST_LOW'] - data[main_feat]['MAG_BEST_HI']))
        print "\n* Apparent size (Re): {0:.0f}+-{1:.0f}, {2:.0f}+-{3:.0f}, {4:.0f}+-{5:.0f}, {6:.0f}+-{7:.0f}, ".format(np.mean(data[subsample]['Re']),np.std(data[subsample]['Re']),np.mean(data[main]['Re']),np.std(data[main]['Re']),np.mean(data[main_bulge]['Re']),np.std(data[main_bulge]['Re']),np.mean(data[main_feat]['Re']),np.std(data[main_feat]['Re']))
        print "\n* Redshift: {0:.2f}+-{1:.2f}, {2:.2f}+-{3:.2f}, {4:.2f}+-{5:.2f}, {6:.2f}+-{7:.2f}, ".format(np.mean(data[subsample]['Z_BEST']),np.std(data[subsample]['Z_BEST']),np.mean(data[main]['Z_BEST']),np.std(data[main]['Z_BEST']),np.mean(data[main_bulge]['Z_BEST']),np.std(data[main_bulge]['Z_BEST']),np.mean(data[main_feat]['Z_BEST']),np.std(data[main_feat]['Z_BEST']))
        print "\n* f_clumpy: {0:.2f}+-{1:.2f}, {2:.2f}+-{3:.2f}, {4:.2f}+-{5:.2f}, {6:.2f}+-{7:.2f}, ".format(np.mean(data[subsample][pclumpy]),np.std(data[subsample][pclumpy]),np.mean(data[main][pclumpy]),np.std(data[main][pclumpy]),np.mean(data[main_bulge][pclumpy]),np.std(data[main_bulge][pclumpy]),np.mean(data[main_feat][pclumpy]),np.std(data[main_feat][pclumpy]))
        print "\n* f_odd: {0:.2f}+-{1:.2f}, {2:.2f}+-{3:.2f}, {4:.2f}+-{5:.2f}, {6:.2f}+-{7:.2f}, ".format(np.mean(data[subsample][podd]),np.std(data[subsample][podd]),np.mean(data[main][podd]),np.std(data[main][podd]),np.mean(data[main_bulge][podd]),np.std(data[main_bulge][podd]),np.mean(data[main_feat][podd]),np.std(data[main_feat][podd]))
        print "\n* f_merger: {0:.2f}+-{1:.2f}, {2:.2f}+-{3:.2f}, {4:.2f}+-{5:.2f}, {6:.2f}+-{7:.2f}, ".format(np.mean(data[subsample][pmerger]),np.std(data[subsample][pmerger]),np.mean(data[main][pmerger]),np.std(data[main][pmerger]),np.mean(data[main_bulge][pmerger]),np.std(data[main_bulge][pmerger]),np.mean(data[main_feat][pmerger]),np.std(data[main_feat][pmerger]))
        print "\n* Sersic index: {0:.1f}+-{1:.1f}, {2:.1f}+-{3:.1f}, {4:.1f}+-{5:.1f}, {6:.1f}+-{7:.1f}, ".format(np.mean(data[subsample]['n']),np.std(data[subsample]['n']),np.mean(data[main]['n']),np.std(data[main]['n']),np.mean(data[main_bulge]['n']),np.std(data[main_bulge]['n']),np.mean(data[main_feat]['n']),np.std(data[main_feat]['n']))
        print "\n* f_features: {0:.1f}+-{1:.1f}, {2:.1f}+-{3:.1f}, {4:.1f}+-{5:.1f}, {6:.1f}+-{7:.1f}, ".format(np.mean(data[subsample][pfeat]),np.std(data[subsample][pfeat]),np.mean(data[main][pfeat]),np.std(data[main][pfeat]),np.mean(data[main_bulge][pfeat]),np.std(data[main_bulge][pfeat]),np.mean(data[main_feat][pfeat]),np.std(data[main_feat][pfeat]))

        slen = sum(subsample)
        figimg = plt.figure(figsize=(12,12))
        n = 5
        gs=gridspec.GridSpec(n,n)
        gs.update(wspace=0.01)
        gs.update(hspace=0.01)
        for i in range(n*n):
            ax=plt.subplot(gs[i/n,i%n])
            gal = data[subsample][int(random.random()*slen)]
            plt.imshow(get_image(gal['survey_id']))
            plt.tick_params(labelbottom='off',labelleft='off')
            ax.annotate('$\mathrm{n = %s}$' % gal['n'],fontsize=15,xy=(0.02,.97),
                xycoords='axes fraction',verticalalignment='top',color='white')
        
            ax.annotate('$\mathrm{f_{feat} = %s}$'%round(gal[pfeat],2),fontsize=15,xy=(0.02,.02),
                xycoords='axes fraction',color='white')
        
        
        figimg.text(0.5,0.92,'GEMS: GZH featured but GALFIT bulgy',fontsize=20,ha='center')

    notstar = (data[pstar] < 0.2)
    subsample2 = (constrained & maglim & notstar & discs & gzh_smooth)
    main = (constrained & maglim & notstar)
    main_discs = (constrained & maglim & notstar & discs)
    main_smooth = (constrained & maglim & notstar & gzh_smooth)

    if do_all:
        print "\n Exploring GALFIT discs but GZH smooth ({0} galaxies)".format(sum(subsample2))
        print "*************************************************************"
        print "*"
        print "* (V-I) color: {0:.2f}+-{1:.2f}, {2:.2f}+-{3:.2f}, {4:.2f}+-{5:.2f}, {6:.2f}+-{7:.2f}, ".format(np.mean(data[subsample2]['MAG_BEST_LOW'] - data[subsample2]['MAG_BEST_HI']),np.std(data[subsample2]['MAG_BEST_LOW'] - data[subsample2]['MAG_BEST_HI']),np.mean(data[main]['MAG_BEST_LOW'] - data[main]['MAG_BEST_HI']),np.std(data[main]['MAG_BEST_LOW'] - data[main]['MAG_BEST_HI']),np.mean(data[main_discs]['MAG_BEST_LOW'] - data[main_discs]['MAG_BEST_HI']),np.std(data[main_discs]['MAG_BEST_LOW'] - data[main_discs]['MAG_BEST_HI']),np.mean(data[main_smooth]['MAG_BEST_LOW'] - data[main_smooth]['MAG_BEST_HI']),np.std(data[main_smooth]['MAG_BEST_LOW'] - data[main_smooth]['MAG_BEST_HI']))
        print "\n* Apparent size (Re): {0:.0f}+-{1:.0f}, {2:.0f}+-{3:.0f}, {4:.0f}+-{5:.0f}, {6:.0f}+-{7:.0f}, ".format(np.mean(data[subsample2]['Re']),np.std(data[subsample2]['Re']),np.mean(data[main]['Re']),np.std(data[main]['Re']),np.mean(data[main_discs]['Re']),np.std(data[main_discs]['Re']),np.mean(data[main_smooth]['Re']),np.std(data[main_smooth]['Re']))
        print "\n* Redshift: {0:.2f}+-{1:.2f}, {2:.2f}+-{3:.2f}, {4:.2f}+-{5:.2f}, {6:.2f}+-{7:.2f}, ".format(np.mean(data[subsample2]['Z_BEST']),np.std(data[subsample2]['Z_BEST']),np.mean(data[main]['Z_BEST']),np.std(data[main]['Z_BEST']),np.mean(data[main_discs]['Z_BEST']),np.std(data[main_discs]['Z_BEST']),np.mean(data[main_smooth]['Z_BEST']),np.std(data[main_smooth]['Z_BEST']))
        print "\n* f_clumpy: {0:.2f}+-{1:.2f}, {2:.2f}+-{3:.2f}, {4:.2f}+-{5:.2f}, {6:.2f}+-{7:.2f}, ".format(np.mean(data[subsample2][pclumpy]),np.std(data[subsample2][pclumpy]),np.mean(data[main][pclumpy]),np.std(data[main][pclumpy]),np.mean(data[main_discs][pclumpy]),np.std(data[main_discs][pclumpy]),np.mean(data[main_smooth][pclumpy]),np.std(data[main_smooth][pclumpy]))
        print "\n* f_odd: {0:.2f}+-{1:.2f}, {2:.2f}+-{3:.2f}, {4:.2f}+-{5:.2f}, {6:.2f}+-{7:.2f}, ".format(np.mean(data[subsample2][podd]),np.std(data[subsample2][podd]),np.mean(data[main][podd]),np.std(data[main][podd]),np.mean(data[main_discs][podd]),np.std(data[main_discs][podd]),np.mean(data[main_smooth][podd]),np.std(data[main_smooth][podd]))
        print "\n* f_merger: {0:.2f}+-{1:.2f}, {2:.2f}+-{3:.2f}, {4:.2f}+-{5:.2f}, {6:.2f}+-{7:.2f}, ".format(np.mean(data[subsample2][pmerger]),np.std(data[subsample2][pmerger]),np.mean(data[main][pmerger]),np.std(data[main][pmerger]),np.mean(data[main_discs][pmerger]),np.std(data[main_discs][pmerger]),np.mean(data[main_smooth][pmerger]),np.std(data[main_smooth][pmerger]))
        print "\n* Sersic index: {0:.1f}+-{1:.1f}, {2:.1f}+-{3:.1f}, {4:.1f}+-{5:.1f}, {6:.1f}+-{7:.1f}, ".format(np.mean(data[subsample2]['n']),np.std(data[subsample2]['n']),np.mean(data[main]['n']),np.std(data[main]['n']),np.mean(data[main_discs]['n']),np.std(data[main_discs]['n']),np.mean(data[main_smooth]['n']),np.std(data[main_smooth]['n']))
        print "\n* f_smooth: {0:.1f}+-{1:.1f}, {2:.1f}+-{3:.1f}, {4:.1f}+-{5:.1f}, {6:.1f}+-{7:.1f}, ".format(np.mean(data[subsample2][pfeat]),np.std(data[subsample2][pfeat]),np.mean(data[main][pfeat]),np.std(data[main][pfeat]),np.mean(data[main_discs][pfeat]),np.std(data[main_discs][pfeat]),np.mean(data[main_smooth][pfeat]),np.std(data[main_feat][pfeat]))

        slen = sum(subsample2)
        figimg = plt.figure(figsize=(12,12))
        n = 5
        gs=gridspec.GridSpec(n,n)
        gs.update(wspace=0.01)
        gs.update(hspace=0.01)
        for i in range(n*n):
            ax=plt.subplot(gs[i/n,i%n])
            gal = data[subsample2][int(random.random()*slen)]
            try:
                plt.imshow(get_image(gal['survey_id']))
                plt.tick_params(labelbottom='off',labelleft='off')
                ax.annotate('$\mathrm{n = %s}$' % gal['n'],fontsize=15,xy=(0.02,.97),
                    xycoords='axes fraction',verticalalignment='top',color='white')
        
                ax.annotate('$\mathrm{f_{feat} = %s}$'%round(gal[pfeat],2),fontsize=15,xy=(0.02,.02),
                    xycoords='axes fraction',color='white')
            except IOError:
                print "Could not find {0}".format(survey_id)
        
        
        figimg.text(0.5,0.92,'GEMS: GZH smooth but GALFIT discs',fontsize=20,ha='center')


    # What are clumpy galaxies like?

    clumpy_sample = (constrained & maglim) & (data[pfeat] > 0.23) & ((data['t12_clumpy_a01_yes_weighted_fraction'] *data['t12_clumpy_total_weight']) >= 20) & (data[pclumpy] > 0.5)

    print "Clumpy galaxies  n: {0:.2f}+-{1:.2f}".format(np.mean(data[clumpy_sample]['n']),np.std(data[clumpy_sample]['n']))
    print "Clumpy galaxies Re: {0:.0f}+-{1:.0f}".format(np.mean(data[clumpy_sample]['Re']),np.std(data[clumpy_sample]['Re']))

    # Plot Figure 16 from Haeussler+07


    fig = plt.figure(figsize=(10,6))

    ax1 = fig.add_subplot(121)
    ax1.scatter(data[constrained & discs]['mag_850lp'],
        np.log10(data[constrained & discs]['Re']*0.03),
        marker='.',
        alpha=0.3,
        color='#377EB8',
        label=r'$n>2.5$')
    ax1.scatter(data[constrained & bulges]['mag_850lp'],
        np.log10(data[constrained & bulges]['Re']*0.03),
        marker='.',
        alpha=0.3,
        color='#E41A1C',
        label=r'$n<2.5$')
    ax1.scatter(data[clumpy_sample]['mag_850lp'],
        np.log10(data[clumpy_sample]['Re']*0.03),
        marker='o',
        color='#4DAF4A',
        label='clumpy')
    ax1.legend(loc='upper right',fontsize=12)
    
    ax1.set_xlabel(r'$m_\mathrm{850LP}$ [mag]',fontsize=20)
    ax1.set_ylabel(r'$\log{R_e}$ [arcsec]',fontsize=20)
    ax1.set_xlim(15,27)
    ax1.set_ylim(-2,1.5)

    ax2 = fig.add_subplot(122)
    ax2.scatter(data[constrained & discs]['MU_HI'],
        data[constrained & discs]['n'],
        marker='.',
        alpha=0.3,
        color='#377EB8',
        label=r'$n>2.5$')
    ax2.scatter(data[constrained & bulges]['MU_HI'],
        data[constrained & bulges]['n'],
        marker='.',
        alpha=0.3,
        color='#E41A1C',
        label=r'$n<2.5$')
    ax2.scatter(data[clumpy_sample]['MU_HI'],
        data[clumpy_sample]['n'],
        marker='o',
        color='#4DAF4A',
        label='clumpy')
    ax2.legend(loc='upper left',fontsize=12)
    
    ax2.set_xlabel(r'$\mu_\mathrm{850LP}$ [mag arcsec$^{-2}$]',fontsize=20)
    ax2.set_ylabel(r'$n$',fontsize=25)
    ax2.set_xlim(10,28)
    ax2.set_ylim(0,8)
    
    fig.tight_layout()

    print "{0:.2f} percent have n > 2.5".format(sum(data[clumpy_sample]['n'] > 2.5) * 1./sum(clumpy_sample))

    plt.show()

def aegis():

    data = fits.getdata('{0}/comparisons/AEGIS/aegis_gzh.fits'.format(gzh_path),1)

    # Limit to galaxies with good morphology measurements from Lotz+08

    flag_i = (data['FLAG_I'] == 0)
    flag_v = (data['FLAG_V'] == 0)
    sn_v = (data['SN_V'] >= 3.)
    sn_i = (data['SN_I'] >= 3.)

    goodmorph = (flag_i & flag_v & sn_v & sn_i)

    gzh_feat = (data[pfeat] >= 0.8)
    gzh_int = (data[pfeat] < 0.8) & (data[pfeat] > 0.2)
    gzh_smooth = (data[pfeat] <= 0.2)

    # Plot G vs. M20 (use I-band data)

    # Compare to Figure 4 in Snyder+15

    fig,((ax1,ax2),(ax3,ax4)) = plt.subplots(2,2,sharex='col',sharey='row',figsize=(10,10))

    bins_x = np.linspace(-3,0,40)
    bins_y = np.linspace(0,1,40)
    levels=np.linspace(3,50,8)


    hr1 = np.histogram2d(data[goodmorph & gzh_feat]['M20_I'],data[goodmorph & gzh_feat]['G_I'],bins=(bins_x,bins_y))
    fi = gaussian_filter(hr1[0].T,0.5)
    CS = ax1.contour(bins_x[1:],bins_y[1:],fi,levels,colors=blue,linewidths=1.5)
    CS.collections[0].set_label('features ({0})'.format(sum(goodmorph & gzh_feat)))

    hr2 = np.histogram2d(data[goodmorph & gzh_int]['M20_I'],data[goodmorph & gzh_int]['G_I'],bins=(bins_x,bins_y))
    fi = gaussian_filter(hr2[0].T,0.5)
    CS = ax2.contour(bins_x[1:],bins_y[1:],fi,levels,colors=green,linewidths=1.5)
    CS.collections[0].set_label('intermediate ({0})'.format(sum(goodmorph & gzh_int)))

    hr3 = np.histogram2d(data[goodmorph & gzh_smooth]['M20_I'],data[goodmorph & gzh_smooth]['G_I'],bins=(bins_x,bins_y))
    fi = gaussian_filter(hr3[0].T,0.5)
    CS = ax3.contour(bins_x[1:],bins_y[1:],fi,levels,colors=red,linewidths=1.5)
    CS.collections[0].set_label('smooth ({0})'.format(sum(goodmorph & gzh_smooth)))

    hr4 = np.histogram2d(data[goodmorph]['M20_I'],data[goodmorph]['G_I'],bins=(bins_x,bins_y))
    fi = gaussian_filter(hr4[0].T,0.5)
    CS = ax4.contour(bins_x[1:],bins_y[1:],fi,np.linspace(3,200,20),colors='grey',linewidths=1.5)
    CS.collections[0].set_label('all ({0})'.format(sum(goodmorph)))

    for ax in (ax1,ax2,ax3,ax4):
        ax.set_xlim(-0.4,-3.2)
        ax.set_ylim(0.29,0.75)
        ax.legend(loc='upper left',fontsize=15)

        # Overplot rough early-type/late-type and merger separation lines from L08

        xarr = np.linspace(-4,0,100)

        m1 = 0.14
        b1 = 0.80

        m2 = -0.14
        b2 = 0.33

        M20_0 = (b2-b1) / (m1-m2)
        G_0 = m1*M20_0 + b1
        xtrim = np.linspace(-4,M20_0,100)

        ax.plot(xtrim,m1*xtrim + b1,linestyle='--',color='black')
        ax.plot(xarr,m2*xarr + b2,linestyle='--',color='green')

    [ax.set_xlabel(r'$M_{20}$',fontsize=30) for ax in (ax3,ax4)]
    [ax.set_ylabel(r'$G$',fontsize=30) for ax in (ax3,ax1)]

    fig.tight_layout()

    # What are the percentages?

    gm20_mergers_I = (data['G_I'] > (m2*data['M20_I'] + b2)) & goodmorph
    gm20_es0sa_I = (data['G_I'] <= (m2*data['M20_I'] + b2)) & (data['G_I'] > (m1*data['M20_I'] + b1)) & goodmorph
    gm20_sbscir_I = (data['G_I'] <= (m2*data['M20_I'] + b2)) & (data['G_I'] <= (m1*data['M20_I'] + b1)) & goodmorph

    gzh_es0sa_I = data[gm20_es0sa_I][pfeat]
    gzh_sbscir_I = data[gm20_sbscir_I][pfeat]

    gm20_mergers_V = (data['G_V'] > (m2*data['M20_V'] + b2)) & goodmorph
    gm20_es0sa_V = (data['G_V'] <= (m2*data['M20_V'] + b2)) & (data['G_V'] > (m1*data['M20_V'] + b1)) & goodmorph
    gm20_sbscir_V = (data['G_V'] <= (m2*data['M20_V'] + b2)) & (data['G_V'] <= (m1*data['M20_V'] + b1)) & goodmorph

    gzh_es0sa_V = data[gm20_es0sa_V][pfeat]
    gzh_sbscir_V = data[gm20_sbscir_V][pfeat]

    print "\nG-M20 (I-band) E/S0/Sa: {0:.2f} f_feat < 0.5, {1:.2f} f_feat >= 0.5.".format(sum(gzh_es0sa_I < 0.5) * 1./sum(gm20_es0sa_I),sum(gzh_es0sa_I >= 0.5) * 1./sum(gm20_es0sa_I))
    print "\nG-M20 (I-band) Sb/Sc/Ir: {0:.2f} f_feat < 0.5, {1:.2f} f_feat >= 0.5.".format(sum(gzh_sbscir_I < 0.5) * 1./sum(gm20_sbscir_I),sum(gzh_sbscir_I >= 0.5) * 1./sum(gm20_sbscir_I))
    print "\nG-M20 (V-band) E/S0/Sa: {0:.2f} f_feat < 0.5, {1:.2f} f_feat >= 0.5.".format(sum(gzh_es0sa_V < 0.5) * 1./sum(gm20_es0sa_V),sum(gzh_es0sa_V >= 0.5) * 1./sum(gm20_es0sa_V))
    print "\nG-M20 (V-band) Sb/Sc/Ir: {0:.2f} f_feat < 0.5, {1:.2f} f_feat >= 0.5.".format(sum(gzh_sbscir_V < 0.5) * 1./sum(gm20_sbscir_V),sum(gzh_sbscir_V >= 0.5) * 1./sum(gm20_sbscir_V))

    fig,(ax1,ax2) = plt.subplots(1,2,figsize=(10,5))

    ax1.hist(data[gm20_sbscir_V][pfeat],bins=20,histtype='stepfilled',color=blue,label='Sb/Sc/Ir')
    ax1.hist(data[gm20_es0sa_V][pfeat],bins=20,histtype='stepfilled',color=red,label='E/S0/Sa')
    ax1.hist(data[gm20_mergers_V][pfeat],bins=20,histtype='step',lw=3,color=green,label='mergers')
    ax1.set_title("V-band")

    ax1.legend(loc='upper right',fontsize=16)

    ax2.hist(data[gm20_sbscir_I][pfeat],bins=20,histtype='stepfilled',color=blue,label='Sb/Sc/Ir')
    ax2.hist(data[gm20_es0sa_I][pfeat],bins=20,histtype='stepfilled',color=red,label='E/S0/Sa')
    ax2.hist(data[gm20_mergers_I][pfeat],bins=20,histtype='step',lw=3,color=green,label='mergers')
    ax2.set_title("I-band")

    for ax in (ax1,ax2):
        ax.set_xlim(0,1)
        ax.set_xlabel(r'$f_\mathrm{features,best}$',fontsize=20)
        ax.set_ylabel('Count',fontsize=20)

    fig.tight_layout()

    # Mergers
    fig,(ax1,ax2) = plt.subplots(1,2,figsize=(10,5))

    ax1.hist(data[gm20_es0sa_V][pmerger],bins=20,normed=True,histtype='step',lw=3,color=red,label='G-M20 E/S0/Sa')
    ax1.hist(data[gm20_sbscir_V][pmerger],bins=20,normed=True,histtype='step',lw=3,color=blue,label='G-M20 Sb/Sc/Ir')
    ax1.hist(data[gm20_mergers_V][pmerger],bins=20,normed=True,histtype='step',lw=3,color=green,label='G-M20 mergers')
    ax1.set_title("V-band")

    ax1.legend(loc='upper right',fontsize=16)

    ax2.hist(data[gm20_es0sa_I][pmerger],bins=20,normed=True,histtype='step',lw=3,color=red,label='G-M20 E/S0/Sa')
    ax2.hist(data[gm20_sbscir_I][pmerger],bins=20,normed=True,histtype='step',lw=3,color=blue,label='G-M20 Sb/Sc/Ir')
    ax2.hist(data[gm20_mergers_I][pmerger],bins=20,normed=True,histtype='step',lw=3,color=green,label='G-M20 mergers')
    ax2.set_title("I-band")

    for ax in (ax1,ax2):
        ax.set_xlim(0,1)
        ax.set_xlabel(r'$f_\mathrm{merger}$',fontsize=20)
        ax.set_ylabel('Count',fontsize=20)

    fig.tight_layout()

    # What percentage of mergers in GZH are in the region marked?

    gzh_mergers = (data[pmerger] > 0.5) & ((data['t08_odd_feature_a06_merger_weighted_fraction'] * data['t08_odd_feature_total_weight']) >= 20)

    print "{0} total mergers in GZH (f_merger > 0.5) with good morphologies.".format(sum(gzh_mergers & goodmorph))
    print "Fraction of GZH mergers in G-M20 space (V-band): {0:.0f}%".format(sum(gzh_mergers & gm20_mergers_V) * 1./ sum(gzh_mergers & goodmorph) * 100.)
    print "Fraction of GZH mergers in G-M20 space (I-band): {0:.0f}%".format(sum(gzh_mergers & gm20_mergers_I) * 1./ sum(gzh_mergers & goodmorph) * 100.)

    def F(G,M20):
        return abs(-0.693*M20_I + 4.95*G_I - 3.85)

    M20_I = data[gzh_mergers]['M20_I']
    G_I = data[gzh_mergers]['G_I']
    absF_I = [F(G,M20) if G >= (0.14*M20 + 0.778) else -1*F(G,M20) for G,M20 in zip(G_I,M20_I)]
    M20_V = data[gzh_mergers]['M20_V']
    G_V = data[gzh_mergers]['G_V']
    absF_V = [F(G,M20) if G >= (0.14*M20 + 0.778) else -1*F(G,M20) for G,M20 in zip(G_V,M20_V)]

    print "F value (V-band) for GZH mergers: {0:.1f} +- {1:.1f}".format(np.mean(absF_V),np.std(absF_V))
    print "F value (I-band) for GZH mergers: {0:.1f} +- {1:.1f}".format(np.mean(absF_I),np.std(absF_I))

    # What number of Lotz mergers have f_merger > 0.2?

    print "{0}/{1} have pmerger < 0.2 (I-band). Avg = {2:.2f}+-{3:.2f}".format(sum(data[gm20_mergers_I][pmerger] < 0.2),sum(gm20_mergers_I),np.mean(data[gm20_mergers_I][pmerger]),np.std(data[gm20_mergers_I][pmerger]))
    print "{0}/{1} have pmerger < 0.2 (V-band). Avg = {2:.2f}+-{3:.2f}".format(sum(data[gm20_mergers_V][pmerger] < 0.2),sum(gm20_mergers_V),np.mean(data[gm20_mergers_V][pmerger]),np.std(data[gm20_mergers_V][pmerger]))
    
    # Plot the GZH mergers in G-M20 space

    fig,ax1 = plt.subplots(1,1,figsize=(6,6))

    hr4 = np.histogram2d(data[goodmorph]['M20_I'],data[goodmorph]['G_I'],bins=(bins_x,bins_y))
    fi = gaussian_filter(hr4[0].T,0.5)
    CS = ax1.contour(bins_x[1:],bins_y[1:],fi,np.linspace(3,200,20),colors='grey',linewidths=1.5)
    CS.collections[0].set_label('all galaxies ({0})'.format(sum(goodmorph)))

    ax1.scatter(data[goodmorph & gzh_mergers]['M20_I'],
                data[goodmorph & gzh_mergers]['G_I'],
                marker = 'o',
                color=green,
                label='GZH mergers ({0})'.format(sum(goodmorph & gzh_mergers)))

    ax1.set_xlim(-0.4,-3.2)
    ax1.set_ylim(0.29,0.75)
    ax1.legend(loc='upper left',fontsize=15)
    
    ax1.plot(xtrim,m1*xtrim + b1,linestyle='--',color='black')
    ax1.plot(xarr,m2*xarr + b2,linestyle='--',color='green')

    ax1.set_xlabel(r'$M_{20}$',fontsize=30)
    ax1.set_ylabel(r'$G$',fontsize=30)

    # What's the merger fraction as a function of redshift?

    # G-M20

    zvalues = (0.2,0.4,0.6,0.8,1.0)

    print ''
    for zlow in zvalues:
        zhi = zlow+0.2
        zbin = (data['Z_BEST'] >= zlow) & (data['Z_BEST'] < zhi)
        print "For {0:.1f}<z<{1:.1f}: GZH merger fraction is {2:.3f}".format(zlow,zhi,sum(zbin & gzh_mergers) * 1./sum(zbin))

    print ''
    for zlow in zvalues:
        zhi = zlow+0.2
        zbin = (data['Z_BEST'] >= zlow) & (data['Z_BEST'] < zhi)
        print "For {0:.1f}<z<{1:.1f}: G-M20 merger fraction (I-band) is {2:2d}%".format(zlow,zhi,int(sum(zbin & gm20_mergers_I) * 100./sum(zbin)))

    print ''
    for zlow in zvalues:
        zhi = zlow+0.2
        zbin = (data['Z_BEST'] >= zlow) & (data['Z_BEST'] < zhi)
        print "For {0:.1f}<z<{1:.1f}: G-M20 merger fraction (V-band) is {2:2d}%".format(zlow,zhi,int(sum(zbin & gm20_mergers_V) * 100./sum(zbin)))



    """

    # Skip for now - need external hard drive for images

    merger_imgs = plt.figure(figsize=(12,12))
    n = 10
    gs=gridspec.GridSpec(n,n)
    gs.update(wspace=0.01)
    gs.update(hspace=0.01)
    for i in range(n*n):
        ax=plt.subplot(gs[i/n,i%n])
        try:
            gal = data[gzh_mergers][i]
            plt.imshow(get_image(gal['survey_id']))
            plt.tick_params(labelbottom='off',labelleft='off')
            ax.annotate('$\mathrm{n = %s}$' % gal['n'],fontsize=15,xy=(0.02,.97),
                xycoords='axes fraction',verticalalignment='top',color='white')
    
            ax.annotate('$\mathrm{f_{merger} = %s}$'%round(gal[pmerger],2),fontsize=15,xy=(0.02,.02),
                xycoords='axes fraction',color='white')

    merger_imgs.text(0.5,0.92,'GEMS: GZH featured but GALFIT bulgy',fontsize=20,ha='center')
    
    """

    plt.show()

def goods():

    goods_n = fits.getdata('{0}/comparisons/GOODS/goods_n_gzh.fits'.format(gzh_path),1)
    goods_s = fits.getdata('{0}/comparisons/GOODS/goods_s_gzh.fits'.format(gzh_path),1)

    gzh_n_feat = (goods_n[pfeat] >= 0.8)
    gzh_n_int = (goods_n[pfeat] < 0.8) & (goods_n[pfeat] > 0.2)
    gzh_n_smooth = (goods_n[pfeat] <= 0.2)
    gzh_s_feat = (goods_s[pfeat] >= 0.8)
    gzh_s_int = (goods_s[pfeat] < 0.8) & (goods_s[pfeat] > 0.2)
    gzh_s_smooth = (goods_s[pfeat] <= 0.2)

    # Plot

    '''
    fig,(ax1,ax2) = plt.subplots(1,2,figsize=(10,6))

    ax1.scatter(goods_n['MORPH'],goods_n[pfeat],marker='.')
    ax2.scatter(goods_s['MORPH'],goods_s[pfeat],marker='.')

    for ax in (ax1,ax2):
        ax.set_xlim(-3,11)
        ax.set_ylim(-0.05,1.05)
        #ax.legend(loc='upper left',fontsize=15)

        ax.set_ylabel(r'$f_\mathrm{features}$',fontsize=20)
        ax.set_xlabel('MORPH',fontsize=20)

    ax1.set_title("GZH - GOODS-N")
    ax2.set_title("GZH - GOODS-S")

    fig.tight_layout()

    '''

    md = {-2:'star',
          -1:'compact',  
           0:'E',  
           1:'E/S0',  
           2:'S0',  
           3:'Sab',  
           4:'S',  
           5:'Scd',  
           6:'Irr',  
           7:'unclass',  
           8:'merger',  
           9:'fault'}

    md_gzh = {-2:{'value':pstar,'label':'star'},
          -1:{'value':pfeat,'label':'feat'},
           0:{'value':pfeat,'label':'feat'},
           1:{'value':pfeat,'label':'feat'},
           2:{'value':pfeat,'label':'feat'},
           3:{'value':pfeat,'label':'feat'},
           4:{'value':pspiral,'label':'spiral'},
           5:{'value':pspiral,'label':'spiral'},
           6:{'value':pirr,'label':'irr'},
           7:{'value':pstar,'label':'artifact'},
           8:{'value':pmerger,'label':'merger'},
           9:{'value':pstar,'label':'artifact'}}

    def goods_images(morphno=-2,ns='n'):

        if ns == 'n':
            data = goods_n
        if ns == 's':
            data = goods_s

        subsample = (data['MORPH'] == morphno)
        slen = sum(subsample)
        figimg = plt.figure(figsize=(12,12))
        n = 5
        gs=gridspec.GridSpec(n,n)
        gs.update(wspace=0.01)
        gs.update(hspace=0.01)

        for i in range(min(n*n,slen)):
            gal = data[subsample][i]
            try:
                img = get_image(gal['survey_id'],goods=True,ns=ns)
                foundit = True
            except:
                print "No image found for {0}".format(gal['survey_id'])
                foundit = False
            if foundit:
                ax=plt.subplot(gs[i/n,i%n])
                ax.imshow(img)
                ax.tick_params(labelbottom='off',labelleft='off')
                ax.annotate('{0}'.format(i+1),fontsize=15,xy=(0.02,.97),
                    xycoords='axes fraction',verticalalignment='top',color='white')
                label = md_gzh[morphno]['label']
                value = md_gzh[morphno]['value']
                ax.annotate('$\mathrm{f_{%s} = %.02f}$'%(label,round(gal[value],2)),fontsize=15,xy=(0.02,.02),
                    xycoords='axes fraction',color='white')

        figimg.text(0.5,0.92,'GZH - GOODS-{1}: MORPH={0}'.format(md[morphno],ns.upper()),fontsize=20,ha='center')

        return figimg

    '''
    for n in np.arange(-2,10):
        for ns in ('n','s'):
            figimg = goods_images(n,ns)
            figimg.savefig("{0}/comparisons/GOODS/gzh_goods_{2}_{1:d}.png".format(gzh_path,n,ns.upper()))
            plt.close(figimg)
    
    '''
    
    # Print the avg. f_features fraction for each morphological category

    print "                                         f_feat     f_smooth   f_artifact"
    for n in np.arange(-2,10):
        vf = goods_n[goods_n['MORPH'] == n][pfeat]
        vs = goods_n[goods_n['MORPH'] == n][psmooth]
        va = goods_n[goods_n['MORPH'] == n][pstar]
        print "Average GZH value for {2:3d} ({3:10}): {0:.2f}+-{1:.2f} {4:.2f}+-{5:.2f} {6:.2f}+-{7:.2f} ".format(np.mean(vf),np.std(vf),n,md[n],np.mean(vs),np.std(vs),np.mean(va),np.std(va))

    ef = goods_n[(goods_n['MORPH'] >= 0) & (goods_n['MORPH'] <= 2)][pfeat]
    es = goods_n[(goods_n['MORPH'] >= 0) & (goods_n['MORPH'] <= 2)][psmooth]
    ea = goods_n[(goods_n['MORPH'] >= 0) & (goods_n['MORPH'] <= 2)][pstar]
    sf = goods_n[(goods_n['MORPH'] >= 3) & (goods_n['MORPH'] <= 5)][pfeat]
    ss = goods_n[(goods_n['MORPH'] >= 3) & (goods_n['MORPH'] <= 5)][psmooth]
    sa = goods_n[(goods_n['MORPH'] >= 3) & (goods_n['MORPH'] <= 5)][pstar]
    pf = goods_n[(goods_n['MORPH'] >= 6) & (goods_n['MORPH'] <= 8)][pfeat]
    ps = goods_n[(goods_n['MORPH'] >= 6) & (goods_n['MORPH'] <= 8)][psmooth]
    pa = goods_n[(goods_n['MORPH'] >= 6) & (goods_n['MORPH'] <= 8)][pstar]
    print "\nAverage GZH value for all      E/S0: {0:.2f}+-{1:.2f} {2:.2f}+-{3:.2f} {4:.2f}+-{5:.2f} ".format(np.mean(ef),np.std(ef),np.mean(es),np.std(es),np.mean(ea),np.std(ea))
    print "Average GZH value for all   Spirals: {0:.2f}+-{1:.2f} {2:.2f}+-{3:.2f} {4:.2f}+-{5:.2f} ".format(np.mean(sf),np.std(vf),np.mean(ss),np.std(ss),np.mean(sa),np.std(sa))
    print "Average GZH value for all Peculiars: {0:.2f}+-{1:.2f} {2:.2f}+-{3:.2f} {4:.2f}+-{5:.2f} ".format(np.mean(pf),np.std(pf),np.mean(ps),np.std(ps),np.mean(pa),np.std(pa))

    print "                                         f_feat     f_smooth   f_artifact"
    for n in np.arange(-2,10):
        vf = goods_s[goods_s['MORPH'] == n][pfeat]
        vs = goods_s[goods_s['MORPH'] == n][psmooth]
        va = goods_s[goods_s['MORPH'] == n][pstar]
        print "Average GZH value for {2:3d} ({3:10}): {0:.2f}+-{1:.2f} {4:.2f}+-{5:.2f} {6:.2f}+-{7:.2f} ".format(np.mean(vf),np.std(vf),n,md[n],np.mean(vs),np.std(vs),np.mean(va),np.std(va))

    ef = goods_s[(goods_s['MORPH'] >= 0) & (goods_s['MORPH'] <= 2)][pfeat]
    es = goods_s[(goods_s['MORPH'] >= 0) & (goods_s['MORPH'] <= 2)][psmooth]
    ea = goods_s[(goods_s['MORPH'] >= 0) & (goods_s['MORPH'] <= 2)][pstar]
    sf = goods_s[(goods_s['MORPH'] >= 3) & (goods_s['MORPH'] <= 5)][pfeat]
    ss = goods_s[(goods_s['MORPH'] >= 3) & (goods_s['MORPH'] <= 5)][psmooth]
    sa = goods_s[(goods_s['MORPH'] >= 3) & (goods_s['MORPH'] <= 5)][pstar]
    pf = goods_s[(goods_s['MORPH'] >= 6) & (goods_s['MORPH'] <= 8)][pfeat]
    ps = goods_s[(goods_s['MORPH'] >= 6) & (goods_s['MORPH'] <= 8)][psmooth]
    pa = goods_s[(goods_s['MORPH'] >= 6) & (goods_s['MORPH'] <= 8)][pstar]
    print "\nAverage GZH value for all      E/S0: {0:.2f}+-{1:.2f} {2:.2f}+-{3:.2f} {4:.2f}+-{5:.2f} ".format(np.mean(ef),np.std(ef),np.mean(es),np.std(es),np.mean(ea),np.std(ea))
    print "Average GZH value for all   Spirals: {0:.2f}+-{1:.2f} {2:.2f}+-{3:.2f} {4:.2f}+-{5:.2f} ".format(np.mean(sf),np.std(vf),np.mean(ss),np.std(ss),np.mean(sa),np.std(sa))
    print "Average GZH value for all Peculiars: {0:.2f}+-{1:.2f} {2:.2f}+-{3:.2f} {4:.2f}+-{5:.2f} ".format(np.mean(pf),np.std(pf),np.mean(ps),np.std(ps),np.mean(pa),np.std(pa))
    plt.show()

def cosmos():

    # Fraction of total sample as function of f_features, f_odd

    tasca = fits.getdata('{0}/comparisons/COSMOS/cosmos_tasca_gzh.fits'.format(gzh_path),1)
    cassata = fits.getdata('{0}/comparisons/COSMOS/cosmos_cassata_gzh.fits'.format(gzh_path),1)
    zest = fits.getdata('{0}/comparisons/COSMOS/cosmos_zest_gzh.fits'.format(gzh_path),1)

    dlist = [
             {'data':cassata,
              'label':'Cassata+07',
              'classcol':'AUTOCLASS'},
             {'data':tasca,
              'label':'Tasca+11',
              'classcol':'CLASS_INT'},
             {'data':zest,
              'label':'Scarlata+07',
              'classcol':'TYPE'}]

    fig,axarr = plt.subplots(3,2,sharey=True,figsize=(8,12))
    bins = np.linspace(0,1,20)

    mlabels = ('early','spiral','irregular')
    colors = (red,blue,green)

    for i,(d,(ax1,ax2)) in enumerate(zip(dlist,axarr)):

        early = (d['data'][d['classcol']] == 1)
        spiral = (d['data'][d['classcol']] == 2)
        irr = (d['data'][d['classcol']] == 3)

        # Plot the fraction of each class for all galaxies as a function of f_features and f_odd

        for ax,frac in zip((ax1,ax2),(pfeat,podd)):
            for morph,color,label,ls in zip((early,spiral,irr),colors,mlabels,('-','--','-.')):
                h,edges = np.histogram(d['data'][morph][frac],bins=20)
                h2,edges2 = np.histogram(d['data'][frac],bins=20)

                e_frac = h * 1./h2
                e_frac_err = cameron_err(0.90,h,h2)

                ax.fill_between(bins,e_frac_err[0],e_frac_err[1],facecolor=color,alpha=0.5)
                ax.plot(bins,e_frac,label=label,lw=3,ls=ls,color=color)
        
                ax.set_xlim(-0.05,1.05)
                ax.set_ylim(0,1)
                ax.set_title(d['label'],fontsize=16)
                
                if i == 0:
                    ax.legend(loc='upper left',fontsize=11)

        ax1.set_ylabel('fraction',fontsize=16)
        if i == 2:
            ax1.set_xlabel(r'$f_\mathrm{features}$',fontsize=20)
            ax2.set_xlabel(r'$f_\mathrm{odd}$',fontsize=20)

    fig.tight_layout()

    # Plot features vs odd for each subsample

    bins_x = np.linspace(0,1,25)
    bins_y = np.linspace(0,1,25)
    extent = [bins_x[0],bins_x[-1],bins_y[0],bins_y[-1]]
    cmap = cm.gray_r
    cmap.set_bad('w')
    levels=10**np.linspace(1,3,10)
    
    fig2,axarr = plt.subplots(3,3,sharex=True,sharey=True,figsize=(10,10))
    
    for i,(d,axrow) in enumerate(zip(dlist,axarr.T)):

        early = (d['data'][d['classcol']] == 1)
        spiral = (d['data'][d['classcol']] == 2)
        irr = (d['data'][d['classcol']] == 3)

        for morph,color,label,ax in zip((early,spiral,irr),colors,mlabels,axrow):
            hall,xedges,yedges = np.histogram2d(d['data'][pfeat],d['data'][podd],bins=(bins_x,bins_y))
            hall_masked = ma.masked_array(hall,mask=(hall<10))
            im = ax.imshow(hall_masked.T, cmap=cmap, alpha=1.0, extent=extent, interpolation='nearest', origin='lower',norm=LogNorm())
            hmorph,xedges,yedges = np.histogram2d(d['data'][morph][pfeat],d['data'][morph][podd],bins=(bins_x,bins_y))
            fi = gaussian_filter(hmorph.T,0.5)
            CS = ax.contour(bins_x[:-1],bins_y[:-1],fi,levels,colors=color,linewidths=2)

            ax.set_xlabel(r'$f_\mathrm{features}$',fontsize=20)
            ax.set_ylabel(r'$f_\mathrm{odd}$',fontsize=20)
            ax.set_title("{0} ({1})".format(d['label'],label),fontsize=16)
 
            #cb = plt.colorbar(im,orientation='vertical')
            #cb.set_label('All COSMOS galaxies',fontsize=16)

    fig2.tight_layout()

    plt.show()

def volume_limit(sample,zmax,magcol,appmag):

    # Volume-limit the sample based on the detection limit and redshift
    
    appmag_lim = appmag * u.mag
    
    distmod = WMAP9.distmod(zmax)
    absmag_lim = appmag_lim - distmod

    absmags = sample[magcol] - WMAP9.distmod(sample['Z_BEST']).value

    #print sample['imaging'][0]
    #print 'absmag_lim',absmag_lim
    #print 'appmag_lim',appmag_lim
    #print 'distmod',distmod

    vl_ind = (absmags < absmag_lim.value) & (sample['Z_BEST'] > 0) & (sample['Z_BEST'] < zmax)
    
    return sample[vl_ind]


def frac_features(zint=0.0,zlow=0.0,volume_limited=False):

    # Fraction of total sample as function of f_features, f_odd

    tasca = fits.getdata('{0}/comparisons/COSMOS/cosmos_tasca_gzh.fits'.format(gzh_path),1)
    cassata = fits.getdata('{0}/comparisons/COSMOS/cosmos_cassata_gzh.fits'.format(gzh_path),1)
    zest = fits.getdata('{0}/comparisons/COSMOS/cosmos_zest_gzh.fits'.format(gzh_path),1)
    aegis = fits.getdata('{0}/comparisons/AEGIS/aegis_gzh.fits'.format(gzh_path),1)
    gems = fits.getdata('{0}/comparisons/GEMS/gems_gzh.fits'.format(gzh_path),1)
    goods = fits.getdata('{0}/comparisons/GOODS/goods_gzh.fits'.format(gzh_path),1)

    dlist = get_dlist(aegis,gems,goods,cassata,tasca,zest)

    if volume_limited:
        for s in dlist:
            zlim = 1.0
            maglim = 22.5

            d = s['data']
            newdata = volume_limit(d,zlim,s['magcol'],maglim)
            s['data'] = newdata

    fig,axarr = plt.subplots(2,3,sharey=True,figsize=(16,12))

    bnum = 10 if zint > 0 else 20
    bins = np.linspace(0,1,bnum)

    mlabels = ('early','spiral','irregular')
    colors = (red,blue,green)

    for i,(d,ax1) in enumerate(zip(dlist,axarr.ravel())):

        data = d['data']

        if zint > 0.:
            data = data[(data['Z_BEST'] >= zlow) & (data['Z_BEST'] < (zlow+zint))]

        if d['survey'] == 'COSMOS':

            # For ZEST galaxies, include tpye 2.0 (S0 disks) as ellipticals

            if d['citation'] == 'Scarlata+07':
                early = (data[d['classcol']] == 1) | ((data[d['classcol']] == 2) & (data['BULG'] == 0))
                spiral = (data[d['classcol']] == 2) & (data['BULG'] > 0.)
                irr = (data[d['classcol']] == 3)
            else:
                early = (data[d['classcol']] == 1)
                spiral = (data[d['classcol']] == 2)
                irr = (data[d['classcol']] == 3)

        if d['survey'] == 'AEGIS':
            m1,b1 = 0.14,0.80
            m2,b2 = -0.14,0.33

            good = (data['FLAG_I'] == 0) & (data['FLAG_V'] == 0) & (data['SN_V'] >= 3.) & (data['SN_I'] >= 3.)

            early = (data['G_I'] <= (m2*data['M20_I'] + b2)) & (data['G_I'] > (m1*data['M20_I'] + b1)) & good
            spiral = (data['G_I'] <= (m2*data['M20_I'] + b2)) & (data['G_I'] <= (m1*data['M20_I'] + b1)) & good
            irr = (data['G_I'] > (m2*data['M20_I'] + b2)) & good

        if d['survey'] == 'GEMS':
            spiral = (data['n'] < 2.5) & (data['flag_constraint'] == 1)
            early  = (data['n'] > 2.5) & (data['flag_constraint'] == 1)
            irr = (data['n'] < 0) & (data['flag_constraint'] == 1)

        if d['survey'] == 'GOODS':
            early = (data['MORPH'] >= 0) & (data['MORPH'] <= 2)
            spiral = (data['MORPH'] >= 3) & (data['MORPH'] <= 5)
            irr = (data['MORPH'] >= 6) & (data['MORPH'] <= 8)

        for morph,color,label,ls in zip((early,spiral,irr),colors,mlabels,('-','--','-.')):

            if morph.sum() > 0:
                h,edges = np.histogram(data[morph][pfeat],bins=bnum,range=(0,1))
                h2,edges2 = np.histogram(data[np.any([early,spiral,irr],axis=0)][pfeat],bins=bnum,range=(0,1))

                e_frac = h * 1./h2
                e_frac_err = cameron_err(0.90,h,h2)

                ax1.fill_between(bins,e_frac_err[0],e_frac_err[1],facecolor=color,alpha=0.5)
                ax1.plot(bins,e_frac,label="{0} ({1})".format(label,morph.sum()),lw=3,ls=ls,color=color)
        
            ax1.set_xlim(-0.05,1.05)
            ax1.set_ylim(0,1.1)
            ax1.set_title("{0} ({1})".format(d['survey'],d['citation']),fontsize=16)

            ax1.legend(loc='upper right',fontsize=12)
                
        if not (i % 3):
            ax1.set_ylabel('fraction',fontsize=16)
        if i > 2:
            ax1.set_xlabel(r'$f_\mathrm{features,best}$',fontsize=20)

    fig.tight_layout()

    '''
    plt.show()
    '''
    if zint > 0.:
        plt.savefig("{0}/comparisons/all_z{1:.1f}.pdf".format(gzh_path,zlow),dpi=100)
    else:
        plt.savefig("{0}/writeup/figures/comparisons_features.pdf".format(gzh_path),dpi=100)

    return None

def frac_odd(zint=0.0,zlow=0.0,volume_limited=False):

    # Fraction of total sample as function of f_features, f_odd

    tasca = fits.getdata('{0}/comparisons/COSMOS/cosmos_tasca_gzh.fits'.format(gzh_path),1)
    cassata = fits.getdata('{0}/comparisons/COSMOS/cosmos_cassata_gzh.fits'.format(gzh_path),1)
    zest = fits.getdata('{0}/comparisons/COSMOS/cosmos_zest_gzh.fits'.format(gzh_path),1)
    aegis = fits.getdata('{0}/comparisons/AEGIS/aegis_gzh.fits'.format(gzh_path),1)
    gems = fits.getdata('{0}/comparisons/GEMS/gems_gzh.fits'.format(gzh_path),1)
    goods = fits.getdata('{0}/comparisons/GOODS/goods_gzh.fits'.format(gzh_path),1)

    dlist = get_dlist(aegis,gems,goods,cassata,tasca,zest)

    if volume_limited:
        for s in dlist:
            zlim = 1.0
            maglim = 22.5

            d = s['data']
            newdata = volume_limit(d,zlim,s['magcol'],maglim)
            s['data'] = newdata

    fig,axarr = plt.subplots(2,3,sharey=True,figsize=(16,12))

    bnum = 10 if zint > 0 else 20
    bins = np.linspace(0,1,bnum)

    mlabels = ('early','spiral','irregular')
    colors = (red,blue,green)

    for i,(d,ax1) in enumerate(zip(dlist,axarr.ravel())):

        data = d['data']

        if zint > 0.:
            data = data[(data['Z_BEST'] >= zlow) & (data['Z_BEST'] < (zlow+zint))]

        if d['survey'] == 'COSMOS':

            # For ZEST galaxies, include tpye 2.0 (S0 disks) as ellipticals

            if d['citation'] == 'Scarlata+07':
                early = (data[d['classcol']] == 1) | ((data[d['classcol']] == 2) & (data['BULG'] == 0))
                spiral = (data[d['classcol']] == 2) & (data['BULG'] > 0.)
                irr = (data[d['classcol']] == 3)
            else:
                early = (data[d['classcol']] == 1)
                spiral = (data[d['classcol']] == 2)
                irr = (data[d['classcol']] == 3)

        if d['survey'] == 'AEGIS':
            m1,b1 = 0.14,0.80
            m2,b2 = -0.14,0.33

            good = (data['FLAG_I'] == 0) & (data['FLAG_V'] == 0) & (data['SN_V'] >= 3.) & (data['SN_I'] >= 3.)

            early = (data['G_I'] <= (m2*data['M20_I'] + b2)) & (data['G_I'] > (m1*data['M20_I'] + b1)) & good
            spiral = (data['G_I'] <= (m2*data['M20_I'] + b2)) & (data['G_I'] <= (m1*data['M20_I'] + b1)) & good
            irr = (data['G_I'] > (m2*data['M20_I'] + b2)) & good

        if d['survey'] == 'GEMS':
            spiral = (data['n'] < 2.5) & (data['flag_constraint'] == 1)
            early  = (data['n'] > 2.5) & (data['flag_constraint'] == 1)
            irr = (data['n'] < 0) & (data['flag_constraint'] == 1)

        if d['survey'] == 'GOODS':
            early = (data['MORPH'] >= 0) & (data['MORPH'] <= 2)
            spiral = (data['MORPH'] >= 3) & (data['MORPH'] <= 5)
            irr = (data['MORPH'] >= 6) & (data['MORPH'] <= 8)

        # Don't look at p_odd for galaxies dominated by dust lanes or lenses
        odd_ind = np.logical_not((data[podd] >= 0.5) & (data[pdustlane] + data[plens] > 0.5))

        hh,ee = np.histogram(data[odd_ind][podd],bins=bnum,range=(0,1))
        print sum(hh[10:])

        for morph,color,label,ls in zip((early,spiral,irr),colors,mlabels,('-','--','-.')):

            if morph.sum() > 0:
                h,edges = np.histogram(data[morph & odd_ind][podd],bins=bnum,range=(0,1))
                h2,edges2 = np.histogram(data[np.any([early,spiral,irr],axis=0) & odd_ind][podd],bins=bnum,range=(0,1))


                e_frac = h * 1./h2
                e_frac_err = cameron_err(0.90,h,h2)

                ax1.fill_between(bins,e_frac_err[0],e_frac_err[1],facecolor=color,alpha=0.5)
                ax1.plot(bins,e_frac,label="{0} ({1})".format(label,morph.sum()),lw=3,ls=ls,color=color)
        
            ax1.set_xlim(-0.05,1.05)
            ax1.set_ylim(0,1)
            ax1.set_title("{0} ({1})".format(d['survey'],d['citation']),fontsize=16)

            ax1.legend(loc='upper left',fontsize=12)
                
        if not (i % 3):
            ax1.set_ylabel('fraction',fontsize=16)
        if i > 2:
            ax1.set_xlabel(r'$f_\mathrm{odd}$',fontsize=20)

    fig.tight_layout()

    '''
    plt.show()
    '''
    if zint > 0.:
        plt.savefig("{0}/comparisons/all_z{1:.1f}.pdf".format(gzh_path,zlow),dpi=100)
    else:
        plt.savefig("{0}/writeup/figures/comparisons_odd.pdf".format(gzh_path),dpi=100)

    return None

if __name__ == "__main__":

    #zint = 0.2
    #for z in np.arange(0,1,zint):
    #    frac_all(zint,z)

    frac_odd(volume_limited=True)
