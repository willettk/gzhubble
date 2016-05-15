import numpy as np
import time
import warnings
from astropy.io import fits
from astropy.io.fits import Column

warnings.filterwarnings('ignore', category=RuntimeWarning, append=True)

'''
Trim the plotting functions to output only the new classification categories.
'''

def bins():

    data=fits.getdata('../data/ferengi_debiasable_data.fits',1)

    #surface brightness bins
    yedges=np.linspace(np.min(data['mu_max_i']),np.max(data['mu_max_i']),10)

    #simulated redshifts of FERENGI data:
    redshifts=[.3,.4,.5,.6,.7,.8,.9,1.]

    return yedges,redshifts
    
def polyfit_with_fixed_points(n, xs, ys, xf, yf) :
    
    mat = np.empty((n + 1 + len(xf),) * 2)
    vec = np.empty((n + 1 + len(xf),))
    x_n = xs**np.arange(2 * n + 1)[:, None]
    yx_n = np.sum(x_n[:n + 1] * ys, axis=1)
    x_n = np.sum(x_n, axis=1)
    idx = np.arange(n + 1) + np.arange(n + 1)[:, None]
    mat[:n + 1, :n + 1] = np.take(x_n, idx)
    xf_n = xf**np.arange(n + 1)[:, None]
    mat[:n + 1, n + 1:] = xf_n / 2
    mat[n + 1:, :n + 1] = xf_n.T
    mat[n + 1:, n + 1:] = 0
    vec[:n + 1] = yx_n
    vec[n + 1:] = yf
    params = np.linalg.solve(mat, vec)
    
    return params[:n + 1]
    
def calculate_hbins(deg,xs,ys):
    
    eps = 1e-3
    xf = (0+eps,1) # offset by small amount so matrix is non-singular
    yf = (0,1)
    p = polyfit_with_fixed_points(deg,xs,ys,xf,yf)
    poly = np.polynomial.Polynomial(p)
    
    return poly
    
def categorize_ferengi_data():

    '''
    Run zeta method on FERENGI data to determine which galaxies are correctable as function of SB, z, p_features
    '''
    
    data=fits.getdata('../data/ferengi_debiasable_data.fits',1)
    
    #Pick out unique galaxies
    galaxies = set(data['sdss_id'])
    
    yedges, redshifts = bins()

    #Defining lists of p_features at high and low (z=0.3) redshifts at given SB and redshifts. 
    scatter_dct={}
    for z in redshifts:
        for edge in yedges:
            scatter_dct[z,edge,'hi']=[]
            scatter_dct[z,edge,'lo']=[]
            scatter_dct[z,edge,'subj_id']=[]
    
    #In each SB range and simulated redshift, make a list of p_features at that redshift and p_features at z=0.3   
    for i,g in enumerate(galaxies):
        this_gal=(data['sdss_id']==g)
        evos = set(data[this_gal]['sim_evolution'])
        for e in evos:
            this_evo=(data[this_gal]['sim_evolution']==e)
            if len(set(data[this_gal][this_evo]['sim_redshift']))==8: #only want stuff where we have data down to 0.3
                p_at_3=(data[this_gal][this_evo]['sim_redshift']==.3)
                p_features_at_3 = data[this_gal][this_evo][p_at_3]['p_features'][0]#value of p_bar at redshift 0.3
                for row in data[this_gal][this_evo]:
                    for y in range(0,len(yedges)-1):
                        for hi_z in redshifts:
                            if round(row['sim_redshift'],2)==hi_z and row['mu_max_i'] > yedges[y] and row['mu_max_i'] < yedges[y+1]: #now look at high redshift data
                                scatter_dct[hi_z,yedges[y],'hi'].append(row['p_features']) # slap p_features in high list 
                                scatter_dct[hi_z,yedges[y],'lo'].append(p_features_at_3) # put z=0.3 value in low list 
                                scatter_dct[hi_z,yedges[y],'subj_id'].append(row['subject_id'])
                                
    '''
    Approach: fit data with low-order polynomial, 
    then evaluate polynomial near bin points. 
    
    Fit function must pass between (0,0) and (1,1); 
    uses method of Lagrange multipliers detailed
    here: http://stackoverflow.com/questions/15191088/how-to-do-a-polynomial-fit-with-fixed-points
    '''
    
    #track 'uncorrectable' regions of z/mu/p space (shaded)
    p_range_uncorrectable_dct={}
    #track 'correctable regions of z/mu/p space (unshaded *and* at least 5 points in bin) 
    p_range_correctable_dct={}
    #store range of spread of data in each bin
    interval_dct={}
    #redshift bins
    vbins = (0.2,0.5,0.8,1.0)
    
    yedge_int=[8,7,6,5,4,3,2,1,0]
    for z in redshifts:
        p_range_uncorrectable_dct[z] = {}
        p_range_correctable_dct[z] = {}
        interval_dct[z] = {}
        for y in yedge_int:
            xs=scatter_dct[z,yedges[y],'lo']
            ys=scatter_dct[z,yedges[y],'hi']
            p_range_uncorrectable_dct[z][yedges[y]]=[]
            p_range_correctable_dct[z][yedges[y]]=[]
            interval_dct[z][yedges[y]]=[]
            p_uncorrectable_list=[]
            p_correctable_list=[]
            
            if len(xs)>1:
                poly = calculate_hbins(3,xs,ys)
                if poly(vbins)[1] < poly(vbins)[0]:#bump in curve, bin edges out of order. No! 
                    poly = calculate_hbins(2,xs,ys)#redo with 2nd degree polynomial 
            this_dct={}  #stores p_features at z=0.3 value for each bin (xs); use to compute spread in each bin
                        #stored like: this_dct[lower_bin_edge,higher_bin_edge] = [distribution of xs in bin]
    
            bins_list=list(poly(vbins)) if len(xs)>1 else [0] #find appropriate bin edges
            bins_list.sort()
            bins_list.insert(0,0) #insert 0 for bottom bin edge
            #bins_list[0]=0 #bring lowest bin to zero
            for b in range(0,len(bins_list)-1):
                bin_bottom = round(bins_list[b],3)
                bin_top = round(bins_list[b+1],3)
                this_dct[bin_bottom,bin_top]=[]
                for l,val in enumerate(ys): #check p_features at z values
                    if val >= bin_bottom and val < bin_top: # if it falls inside bin in question:
                        this_dct[bin_bottom,bin_top].append(xs[l]) #then put corresponding p_features,z=0.3 value in list
            
            #Now compute the mean and spread of p_features,z=0.3 for each bin
            p_features_means=[]
            x_error_lo=[]
            x_error_hi=[]
            bin_centers=[]
    
            for b in range(0,len(bins_list)-1):
                bin_bottom = round(bins_list[b],3)
                bin_top = round(bins_list[b+1],3)
                p_features_means.append(np.median(this_dct[bin_bottom,bin_top]))
                try:
                    x_error_lo.append(p_features_means[b] - np.percentile(this_dct[bin_bottom,bin_top],25))
                except IndexError,ValueError:
                    x_error_lo.append(0)
                try:
                    x_error_hi.append(np.percentile(this_dct[bin_bottom,bin_top],75) - p_features_means[b])
                except IndexError,ValueError:
                    x_error_hi.append(0)
                bin_centers.append((bin_top-bin_bottom)/2.+bin_bottom)
                try:
                    interval_dct[z][yedges[y]].append({'bin_bottom':bin_bottom,'bin_top':bin_top,'low_limit':np.percentile(this_dct[bin_bottom,bin_top],25),'hi_limit':np.percentile(this_dct[bin_bottom,bin_top],75)})
                except IndexError,ValueError:
                    pass
            
           #now: check width of spread
            for b in range(0,len(bins_list)-1):
                bin_bottom = round(bins_list[b],3)
                bin_top = round(bins_list[b+1],3)
                bin_width = bin_top - bin_bottom
                spread = x_error_lo[b]+x_error_hi[b]
                points_in_bin=len(this_dct[bin_bottom,bin_top])
                try:
                    if spread > bin_width and points_in_bin>=5: #uncorrectable:
                        p_uncorrectable_list.append(bin_bottom)
                        p_uncorrectable_list.append(bin_top)
                    if spread < bin_width and points_in_bin >=5: #correctable! fo sho!
                        p_correctable_list.append(bin_bottom)
                        p_correctable_list.append(bin_top)
                except ValueError:
                    pass
            if len(p_uncorrectable_list)>0: #if square has any uncorrectable bins, store lowest and highest bin edge values 
                p_range_uncorrectable_dct[z][yedges[y]].append(np.min(p_uncorrectable_list))
                p_range_uncorrectable_dct[z][yedges[y]].append(np.max(p_uncorrectable_list))
            else: 
                p_range_uncorrectable_dct[z][yedges[y]].append(0)
                p_range_uncorrectable_dct[z][yedges[y]].append(0)
            if len(p_correctable_list)>0:
                p_range_correctable_dct[z][yedges[y]].append(np.min(p_correctable_list))
                p_range_correctable_dct[z][yedges[y]].append(np.max(p_correctable_list))
            else:
                p_range_correctable_dct[z][yedges[y]].append(0)
                p_range_correctable_dct[z][yedges[y]].append(0)
    
    #correctable range: white region 
    #uncorrectable range: blue region
    #NEI region: pink (or not on plot) 
    
    #get list of debiasable FERENGI galaxies 
    unique_galaxies = set(data['sdss_id'])
    z0ind = np.zeros(len(data),dtype=bool)
    for ug in unique_galaxies:
        ind = (data['sdss_id'] == ug)
        if data[ind]['sim_redshift'].min() < 0.301:
            z0ind[ind] = True
            
    data_z0 = data[z0ind]
    category_list_ferengi=[]
    for row in data_z0:
        if row['mu_max_i'] > yedges[0] and row['mu_max_i'] <= yedges[len(yedges)-1] and  row['sim_redshift'] > redshifts[0]-.05 and row['sim_redshift'] <= redshifts[len(redshifts)-1] + .05: # if within FERENGI space, check where it is. else, consider NEI or uncorrectable.
            for y in range(0,len(yedges)-1):
                if row['mu_max_i'] > yedges[y] and row['mu_max_i'] <= yedges[y+1]:  
                    for i,z in enumerate(redshifts): 
                        if row['sim_redshift'] > redshifts[i]-.05 and row['sim_redshift'] <= redshifts[i] + .05: # pick out where it is in SB/z and check color
                            if row['p_features'] > p_range_correctable_dct[z][yedges[y]][0] and row['p_features'] <= p_range_correctable_dct[z][yedges[y]][1]:# if it's in correctable range::
                                category_list_ferengi.append({'objid':row['subject_id'],'Correctable_Category':'correctable'})
                            elif row['p_features'] > p_range_uncorrectable_dct[z][yedges[y]][0] and row['p_features'] <= p_range_uncorrectable_dct[z][yedges[y]][1]:# if it's in uncorrectable range::
                                category_list_ferengi.append({'objid':row['subject_id'],'Correctable_Category':'uncorrectable'})
                            else: #not in correctable or uncorrectable range, so nei
                                category_list_ferengi.append({'objid':row['subject_id'],'Correctable_Category':'nei'})
        else: #galaxies outside FERENGI SB and z limits - still need to have meaasureable z and SB to possibly correct. 
            if row['sim_redshift'] > 0 and row['sim_redshift'] < 9 and row['mu_max_i'] >0: #these have measurements for z and SB, put in NEI
                category_list_ferengi.append({'objid':row['subject_id'],'Correctable_Category':'nei'})
            else: #these have nan or infinite values of z or mu, put in need_redshift_list
                category_list_ferengi.append({'objid':row['subject_id'],'Correctable_Category':'nei_needs_redshift'})
    
    #create fits file of galaxies with Correctable_Category Label for FERENGI
    strcolumn = np.array([' ']*len(category_list_ferengi),dtype='S24')
    c0 = Column(name='objid', format='A24', array=strcolumn) 
    c01 = Column(name='Correctable_Category',format='A24',array=strcolumn)
    cols = fits.TableHDU.from_columns([c0,c01])
    category_table_ferengi = fits.TableHDU.from_columns(cols.columns)           
    
    for i,gal in enumerate(category_list_ferengi):
        category_table_ferengi.data.field('objid')[i] = gal['objid']
        category_table_ferengi.data.field('Correctable_Category')[i] = gal['Correctable_Category']
    #category_table_ferengi.writeto('ferengi_debiasable_data.fits',clobber=True)

    return p_range_correctable_dct,p_range_uncorrectable_dct,interval_dct

def categorize_hubble_data(p_range_correctable_dct,p_range_uncorrectable_dct,interval_dct):

    '''
    Run zeta method on Hubble data to determine which galaxies are correctable as function of SB, z, p_features
    based on the FERENGI values
    '''
    hubble_data=fits.getdata('../data/Hubble_t01_data.fits',1)
    
    yedges, redshifts = bins()

    # Determine which parts of hubble sample are correctable (white), uncorrectable (green), or nei (purple) 
    
    category_list=[]

    '''
    Possibilities

    1. Galaxy was within FERENGI space, was correctable
    2. Galaxy was within FERENGI space, was uncorrectable
    3. Galaxy was within FERENGI space, had nei

    4. Galaxy was outside FERENGI space, had redshift z < 0.3
    5. Galaxy was outside FERENGI space, had no redshift
    6. Galaxy was outside FERENGI space, had nei

    Trying to run array-wise instead of looping everything (which takes 5-10 minutes). Haven't succeeded, probably not worth the effort.

    '''

    within_ferengi = (hubble_data['MU_HI'] > yedges[0]) & (hubble_data['MU_HI'] <= yedges[len(yedges)-1]) & (hubble_data['Z'] > redshifts[0]) & (hubble_data['Z'] <= redshifts[len(redshifts)-1] + .05)
    outside_ferengi = np.logical_not(within_ferengi)
    import itertools
    yz = list(itertools.product(*[yedges,redshifts]))
    for y,z in yz:
        con1 = [p_range_correctable_dct[entry[1]][entry[0]] for row in hubble_data[within_ferengi] if p_range_correctable_dct[entry[1]][entry[0]][0] <= row['p_features'] < p_range_correctable_dct[entry[1]][entry[0]][1]]



    
    for row in hubble_data: 
        # If galaxy is within FERENGI space, check where it is. else, consider NEI or uncorrectable.
        if row['MU_HI'] > yedges[0] and row['MU_HI'] <= yedges[len(yedges)-1] and  row['Z'] > redshifts[0] and row['Z'] <= redshifts[len(redshifts)-1] + .05:
            for y in range(0,len(yedges)-1):
                if row['MU_HI'] > yedges[y] and row['MU_HI'] <= yedges[y+1]:  
                    for i,z in enumerate(redshifts): 
                        if row['Z'] > redshifts[i]-.05 and row['Z'] <= redshifts[i] + .05: # pick out where it is in SB/z and check color
                            if row['p_features'] > p_range_correctable_dct[z,yedges[y]][0] and row['p_features'] <= p_range_correctable_dct[z,yedges[y]][1]:# if it's in correctable range::
                                category_list.append({'objid':row['OBJNO'],'Correctable_Category':'correctable','Imaging':row['Imaging']})
                            elif row['p_features'] > p_range_uncorrectable_dct[z,yedges[y]][0] and row['p_features'] <= p_range_uncorrectable_dct[z,yedges[y]][1]:# if it's in uncorrectable range::
                                category_list.append({'objid':row['OBJNO'],'Correctable_Category':'uncorrectable','Imaging':row['Imaging']})
                            else: #not in correctable or uncorrectable range, so nei
                                category_list.append({'objid':row['OBJNO'],'Correctable_Category':'nei','Imaging':row['Imaging']})
        else: #galaxies outside FERENGI SB and z limits - still need to have meaasureable z and SB to possibly correct. 
            if row['Z'] >=.3 and row['Z'] < 9 and row['MU_HI'] >0: #these have measurements for z and SB, put in NEI
                category_list.append({'objid':row['OBJNO'],'Correctable_Category':'nei','Imaging':row['Imaging']})
            elif row['Z'] > 0 and row['Z'] < .3: #don't need to be corrected, z < z0
                category_list.append({'objid':row['OBJNO'],'Correctable_Category':'z_lt_3','Imaging':row['Imaging']})
            else: #these have nan or infinite values of z or mu, put in need_redshift_list
                category_list.append({'objid':row['OBJNO'],'Correctable_Category':'nei_needs_redshift','Imaging':row['Imaging']})
    
    low_hi_limit_list=[]
    for row in hubble_data:
        if row['MU_HI'] > yedges[0] and row['MU_HI'] <= yedges[len(yedges)-1] and  row['Z'] > redshifts[0]-.05 and row['Z'] <= redshifts[len(redshifts)-1] + .05: 
            for y in range(0,len(yedges)-1):
                if row['MU_HI'] > yedges[y] and row['MU_HI'] <= yedges[y+1]:  
                    for i,z in enumerate(redshifts): 
                        if row['Z'] > redshifts[i]-.05 and row['Z'] <= redshifts[i] + .05: #now we have mu,z info:
                            for bin_range in interval_dct[z,yedges[y]]:
                                if row['p_features'] >= bin_range['bin_bottom'] and row['p_features'] < bin_range['bin_top']:
                                    low_hi_limit_list.append({'objid':row['OBJNO'],'low_limit':bin_range['low_limit'],'hi_limit':bin_range['hi_limit']})
                            
    imaging_list=set(hubble_data['Imaging'])
    total_correctable=0
    total_uncorrectable=0
    total_z_lt_3=0
    total_nei=0
    total_nr=0
    for survey in imaging_list:
        c=0
        u=0
        z_lt_3=0
        nei=0
        nr=0
        for row in category_list:
            if row['Correctable_Category']=='correctable' and row['Imaging']==survey:
                c+=1
            if row['Correctable_Category']=='uncorrectable' and row['Imaging']==survey:
                u+=1
            if row['Correctable_Category']=='z_lt_3' and row['Imaging']==survey:
                z_lt_3+=1
            if row['Correctable_Category']=='nei' and row['Imaging']==survey:
                nei+=1
            if row['Correctable_Category']=='nei_needs_redshift' and row['Imaging']==survey:
                nr+=1
        total_correctable+=c
        total_uncorrectable+=u
        total_z_lt_3+=z_lt_3
        total_nei+=nei
        total_nr+=nr
        print 'the number of correctable galaxies in %s is %i' %(survey,c)
        print 'the number of uncorrectable galaxies in %s is %i' %(survey,u)
        print 'the number of galaxies with z < 0.3 in %s is %i' %(survey,z_lt_3)
        print 'the number of NEI galaxies in %s is (due to not enough FERENGI galaxies in bin) is %i' %(survey,nei)
        print 'the number of NEI galaxies in %s is (due to needing redshift measurements) is %i' %(survey,nr)
    print 'total correctable: %i' %total_correctable
    print 'total uncorrectable: %i' %total_uncorrectable
    print 'total z less than .3: %i' %total_z_lt_3
    print 'total nei: %i' %total_nei
    print 'total nr: %i' %total_nr
    print 'total: %i' %len(category_list)
    
    #create fits file of galaxies with Correctable_Category Label 
    strcolumn = np.array([' ']*len(category_list),dtype='S24')
    c0 = Column(name='objid', format='A24', array=strcolumn) 
    c01 = Column(name='Correctable_Category',format='A24',array=strcolumn)
    cols = fits.TableHDU.from_columns([c0,c01])
    category_table = fits.TableHDU.from_columns(cols.columns)           
    
    for i,gal in enumerate(category_list):
        category_table.data.field('objid')[i] = gal['objid']
        category_table.data.field('Correctable_Category')[i] = gal['Correctable_Category']
    #category_table.writeto('category_table.fits',clobber=True)
    
    #create fits file of lower and upper limits for p_features
    floatcolumn = np.zeros(len(low_hi_limit_list),dtype=float)
    strcolumn = np.array([' ']*len(low_hi_limit_list),dtype='S24')
    c2 = Column(name='objid', format='A24', array=strcolumn) 
    c3 = Column(name='low_limit',format='D',array=floatcolumn)
    c4 = Column(name='high_limit',format='D',array=floatcolumn)
    
    cols2 = fits.TableHDU.from_columns([c2,c3,c4])
    limit_table = fits.TableHDU.from_columns(cols2.columns)        
    
    for i,gal in enumerate(low_hi_limit_list):
        limit_table.data.field('objid')[i] = gal['objid']
        limit_table.data.field('low_limit')[i] = gal['low_limit']
        limit_table.data.field('high_limit')[i] = gal['hi_limit']
        
    #limit_table.writeto('limits_table.fits',clobber=True)
    
    return None

if __name__ == "__main__":

    start = time.clock()
    p_range_correctable_dct,p_range_uncorrectable_dct,interval_dct = categorize_ferengi_data()
    categorize_hubble_data(p_range_correctable_dct,p_range_uncorrectable_dct,interval_dct)
    end = time.clock()
    print "{:.3f} seconds to run categorize_hubble_data.py".format(end-start)


