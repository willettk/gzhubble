# Compare the measured SNR in the original Hubble images to those in FERENGI
#
# Naively: plot the SNR as a function of galaxy (apparent) magnitude for both GZH and FERENGI

from astropy.stats import signal_to_noise_oir_ccd
from astropy.io import fits
import numpy as np
from matplotlib import pyplot as plt
import pandas as pd

def load_gzh_image():

    # Load a single GZH FITS image

    # First problem - don't have the original GZH FITS images anywhere
    
    img = fits.getdata(filename,0)

    return img

def load_ferengi_image():

    # Load a single FERENGI FITS image

    filename = "../588848899357999327_imout_acsf814_ugriz_z0.800_evo0.0.fits"
    img = fits.getdata(filename,0)

    return img

def example_snr():

    img = load_ferengi_image()

    snr = signal_to_noise_oir_ccd(t, source_eps, sky_eps, dark_eps, rd, npix, gain=1.0)

    return None

# Super naively. What if I just look at the rms of the background for all of the FERENGI
# images, both as a histogram and function of galaxy magnitude?
#
# Pick 25x25 patch in one corner of the image
#
def scrape_rms():

    path = '/Users/willettk/Astronomy/Research/GalaxyZoo/gzhubble'

    import glob

    files = glob.glob("{0}/ferengi_fits/*.fits".format(path))

    f606 = open("{0}/ferengi_rms_f606.csv".format(path),"w")
    f814 = open("{0}/ferengi_rms_f814.csv".format(path),"w")

    for filename in files:
        img = fits.getdata(filename,0)
        corner_patch = img[:25,:25]
        mean = np.mean(corner_patch)
        rms = np.std(corner_patch)

        fname = filename.split("/")[-1]
        objid = fname.split("_")[0]
        band = fname.split("_")[2]

        if band == 'acsf606':
            print >> f606,"{0},{2},{3}".format(objid,band,mean,rms)
        elif band == 'acsf814':
            print >> f814,"{0},{2},{3}".format(objid,band,mean,rms)

    return None

def plot_results():

    fig = plt.figure(1,figsize=(10,10))
    ax1 = fig.add_subplot(221)
    ax2 = fig.add_subplot(222)

    f606 = pd.read_csv("../ferengi_rms_f606.csv",names=('objid','mean','rms'))
    f814 = pd.read_csv("../ferengi_rms_f814.csv",names=('objid','mean','rms'))

    nbins = 50

    ax1.hist(f606['mean'],bins=nbins,range=(-0.002,0.02),histtype='stepfilled',label='f606')
    ax1.hist(f814['mean'],bins=nbins,range=(-0.002,0.02),histtype='step',label='f814')
    ax1.set_xlabel('mean')
    ax1.legend()

    ax2.hist(f606['rms'],bins=nbins,range=(0,0.02),histtype='stepfilled',label='f606')
    ax2.hist(f814['rms'],bins=nbins,range=(0,0.02),histtype='step',label='f814')
    ax2.set_xlabel('rms')
    ax2.legend()

    ax3 = fig.add_subplot(223)
    ax4 = fig.add_subplot(224)

    ax3.scatter(f606['mean'],f606['rms'],label='f606',color='#377EB8',marker='o',alpha=0.4)
    ax3.scatter(f814['mean'],f814['rms'],label='f814',color='#E41A1C',marker='+',alpha=0.4)
    ax3.set_title('f606')
    ax3.set_xlabel('mean')
    ax3.set_ylabel('rms')
    ax3.legend()

    plt.show()

    return None

if __name__ == "__main__":
    example_snr()
