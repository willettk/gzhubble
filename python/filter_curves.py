import os
import urllib2

import numpy as np
from matplotlib import pyplot as plt
from matplotlib.patches import Arrow
from astropy.io import fits
import glob

def fetch_acs_filter(filt):

    url = 'http://www.stsci.edu/hst/acs/analysis/throughputs/tables/wfc_{0}.dat'.format(filt)
    
    if not os.path.exists('downloads'):
        os.makedirs('downloads')

    '''
    loc = os.path.join('downloads', '%s.dat' % filt)
    if not os.path.exists(loc):
        print "downloading from %s" % url
        F = urllib2.urlopen(url)
        open(loc, 'w').write(F.read())

    F = open(loc)
        
    data = np.loadtxt(F)
    '''

    # Picked filter curves instead from http://www.stsci.edu/hst/observatory/crds/SIfileInfo/pysynphottables/current_acs_throughput_html

    filenames = glob.glob("downloads/acs*{0}*wfc*.fits".format(filt))
    if len(filenames) > 0:
        data = fits.getdata(filenames[0],1)

    return data

def fetch_subaru_filter(filt):

    filename = '{0}_subaru'.format(filt)
    url = 'http://www.astro.caltech.edu/~capak/filters/{0}.res'.format(filename)
    
    if not os.path.exists('downloads'):
        os.makedirs('downloads')

    loc = os.path.join('downloads', '%s.dat' % filename)
    if not os.path.exists(loc):
        print "downloading from %s" % url
        F = urllib2.urlopen(url)
        open(loc, 'w').write(F.read())

    F = open(loc)
        
    data = np.loadtxt(F)
    return data


#----------------------------------------------------------------------
# Plot filters in color with a single spectrum
fig,axarr = plt.subplots(nrows=2,ncols=3,figsize=(16,10))

acs_filters = {'F435W':{'color':'b','labelpos':4300},
                'F606W':{'color':'g','labelpos':6000},
                'F775W':{'color':'r','labelpos':7800},
                'F814W':{'color':'m','labelpos':8200},
                'F850LP':{'color':'k','labelpos':9200}}

suprime_filters = {'B':{'label':r'$B_J$','labelpos':(3500,0.50)},
                    'g':{'label':r'$g^+$','labelpos':(5100,0.80)},
                    'r':{'label':r'$r^+$','labelpos':(6300,0.70)}}

gzh = {'AEGIS':('F606W','F814W'),'COSMOS':('F814W',),'GEMS, GOODS-S':('F606W','F850LP'),'GOODS-N':('F606W','F775W'),'deep GOODS-N,-S':('F435W','F606W','F775W','F814W')}
kwargs = dict(fontsize=14, ha='center', va='center')

for survey,ax in zip(gzh,axarr.ravel()):
    for acs_filter in gzh[survey]:
        """
        X = fetch_acs_filter(acs_filter)
        c = acs_filters[acs_filter]['color']
        ax.fill(X[:, 0], X[:, 1], ec=c, fc=c, alpha=0.4)
        ax.text(acs_filters[acs_filter]['labelpos'], 0.02, acs_filter, color=c, **kwargs)
        """

        X = fetch_acs_filter(acs_filter.lower())
        if acs_filter == 'F850LP':
            dlambda = 0.001
            X = np.append(X,np.array([(X['WAVELENGTH'][-1] + dlambda,0.0)],dtype=X.dtype))
        c = acs_filters[acs_filter]['color']
        ax.fill(X['WAVELENGTH'], X['THROUGHPUT'], ec=c, fc=c, alpha=0.4)
        pos = 0.12 if ('F775W' in gzh[survey] and 'F814W' in gzh[survey] and acs_filter == 'F814W') else 0.02
        ax.text(acs_filters[acs_filter]['labelpos'], pos, acs_filter, color=c, **kwargs)

    if survey == 'COSMOS':
        with open('downloads/suprime_QE.txt','r') as fn:
            QE = np.loadtxt(fn)
        with open('downloads/suprime_PFU_transmission.txt','r') as fn:
            PFU = np.loadtxt(fn)
        with open('downloads/suprime_reflectivity.txt','r') as fn:
            RF = np.loadtxt(fn)
        for sf in ('Bgr'):
            X = fetch_subaru_filter(sf)
            transmission = []
            for x,y in zip(X[:, 0], X[:, 1]):
                transmission.append(QE[(np.abs(QE[:,0]-x)).argmin(),1] * PFU[(np.abs(PFU[:,0]-x)).argmin(),1] * RF[(np.abs(RF[:,0]-x)).argmin(),1] * y)
            ax.plot(X[:, 0], X[:, 1], color='k')
            ax.text(suprime_filters[sf]['labelpos'][0],suprime_filters[sf]['labelpos'][1], suprime_filters[sf]['label'], color='k', **kwargs)

    ax.set_xlim(3000, 11000)
    ax.set_ylim(0,1)
    ax.set_xticks((3000,5000,7000,9000,11000))
    ax.set_title(survey,fontsize=20)
    ax.set_xlabel(r'wavelength [$\AA$]')
    ax.set_ylabel('filter transmission')

fig.delaxes(axarr.ravel()[-1])


plt.tight_layout()
plt.show()
