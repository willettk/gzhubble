import os
import urllib2
from astropy.io.votable import parse_single_table
from matplotlib import pyplot as plt

'''
Deprecated locations for filter curves:

ACS:
    http://www.stsci.edu/hst/acs/analysis/throughputs/tables/wfc_{0}.dat'.format(filt)
    http://www.stsci.edu/hst/observatory/crds/SIfileInfo/pysynphottables/current_acs_throughput_html

Suprime:

    http://www.astro.caltech.edu/~capak/filters/{0}.res'.format(filename)

'''
    
def fetch_filter(telescope,instrument,filt):

    # Downloads filter transmission curves as VOTables from Spanish Virtual Observatory

    # Use -81C for HST, operating temperature since 2006
    if instrument == 'ACS_WFC':
        filt += '_81'

    url = 'http://svo2.cab.inta-csic.es/svo/theory/fps3/fps.php?ID={0}/{1}.{2}'.format(telescope,instrument,filt)

    if not os.path.exists('downloads'):
        os.makedirs('downloads')

    loc = os.path.join('downloads', '{0}.{1}.{2}.vot'.format(telescope,instrument,filt))
    if not os.path.exists(loc):
        print "downloading from %s" % url
        F = urllib2.urlopen(url)
        open(loc, 'w').write(F.read())

    table = parse_single_table(loc)
    data = table.array
        
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
                    'SDSS_g':{'label':r'$g^+$','labelpos':(5100,0.80)},
                    'SDSS_r':{'label':r'$r^+$','labelpos':(6300,0.70)}}

gzh = [{'name':'AEGIS','filters':('F606W','F814W')},
        {'name':'COSMOS','filters':('F814W',)},
        {'name':'GEMS, GOODS-S','filters':('F606W','F850LP')},
        {'name':'GOODS-N','filters':('F606W','F775W')},
        {'name':'5-epoch GOODS-N,-S','filters':('F435W','F606W','F775W','F850LP')},
        ]
kwargs = dict(fontsize=14, ha='center', va='center')

efficiency = True
synphot = False

for survey,ax in zip(gzh,axarr.ravel()):
    
    for filt in survey['filters']:

        telescope = 'HST'
        instrument = 'ACS_WFC'
        X = fetch_filter(telescope,instrument,filt)
        c = acs_filters[filt]['color']
        ax.fill(X['Wavelength'], X['Transmission'], ec=c, fc=c, alpha=0.4)
        ax.text(acs_filters[filt]['labelpos'], 0.02, filt, color=c, **kwargs)

    if survey['name'] == 'COSMOS':
        telescope = 'Subaru'
        instrument = 'Suprime'
        for filt in suprime_filters:
            X = fetch_filter(telescope,instrument,filt)
            ax.plot(X['Wavelength'], X['Transmission'], color='k')
            ax.text(suprime_filters[filt]['labelpos'][0],suprime_filters[filt]['labelpos'][1], suprime_filters[filt]['label'], color='k', **kwargs)

    ax.set_xlim(3000, 11000)
    ax.set_ylim(0,1)
    ax.set_xticks((3000,5000,7000,9000,11000))
    ax.set_title(survey['name'],fontsize=20)
    ax.set_xlabel(r'wavelength [$\AA$]',fontsize=16)
    ax.set_ylabel('transmission',fontsize=16)

fig.delaxes(axarr.ravel()[-1])

plt.tight_layout()
plt.show()
