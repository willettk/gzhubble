from astropy.io import fits
from scipy import optimize
from matplotlib import pyplot as plt
import numpy as np

def get_data():
    data = fits.getdata("../ferengi_data.fits",1)
    unique_galaxies = set(data['sdss_id'])

    return data,unique_galaxies

def fzeta_exp(p,x):
    
    y = p[0] * np.exp(-1 * (x-p[1])/p[2])

    return y
    
def fzeta_lin(p,x):
    
    y = p[0] + p[1] * x

    return y

def errfunc(p,x,y):
    
    err = y - fzeta_exp(p,x)

    return err

def common_labels(fig,xlabel=None,ylabel=None,xfontsize=12,yfontsize=12):

    # Set common labels
    cax = fig.add_subplot(111)    # The big subplot
    cax.set_axis_bgcolor('none')
    cax.spines['top'].set_color('none')
    cax.spines['bottom'].set_color('none')
    cax.spines['left'].set_color('none')
    cax.spines['right'].set_color('none')
    cax.tick_params(labelcolor='none', top='off', bottom='off', left='off', right='off')
    cax.set_xlabel(xlabel,fontsize=xfontsize)
    cax.set_ylabel(ylabel,fontsize=yfontsize)

    return cax

def plot_examples(data,unique_galaxies):

    p_guess = np.array([0.5,0.3,1.])
    # Set up plot
    fig,axarr = plt.subplots(nrows=2,ncols=2,sharex=True,sharey=True)
    bigax = common_labels(fig,'Redshift',r'$f_{features}$',14,20)
    zarr = np.linspace(0,1,50)
    
    for ax in axarr.ravel():
        
        slen = 0
        # Make sure there are enough points in the data to fit a function
        while slen < (len(p_guess)+1):
            ind = (data['sdss_id'] == unique_galaxies.pop())
            slen = sum(ind)
            
        galaxy1 = data[ind & (data['sim_evolution'] == 0)]
    
        z_gal = galaxy1['sim_redshift']
        f_gal = galaxy1['p_features']
    
        # Data must be explicitly cast as double-type precision for optimization to work. Incredibly frustrating.
        # Fix: http://stackoverflow.com/questions/12473406/scipy-optimize-leastsq-returns-best-guess-parameters-not-new-best-fit 
    
        p, cov, infodict, mesg, ier = optimize.leastsq(errfunc,p_guess,args=(z_gal.astype(np.float64),f_gal.astype(np.float64)),full_output=1)
        ax.plot(z_gal,f_gal,lw=2)
        ax.plot(zarr,fzeta_exp(p,zarr),'--',lw=1)
        ax.set_xlim(0,1)
        ax.set_ylim(0,1)
    
    
    plt.show()

    return None

def plot_zeta_sb(data):


if __name__ == "__main__":

    data,unique_galaxies = get_data()
    plot_examples(data,unique_galaxies)
