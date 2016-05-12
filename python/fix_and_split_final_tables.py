# Replace correctable_category with "correction_type", mapped to integers

from astropy.io import fits
import sys

assert len(sys.argv) == 3, \
    "python fix_and_split_final_tables.py infile.fits outfile.fits"

infile = sys.argv[1]
outfile = sys.argv[2]

# Load old file

data = fits.getdata(infile,1)

# Get a blank column of the correct length that has dtype = integers (to be overwritten)

newcol2 = data['t18_clumps_embedded_total_count'].copy()

# Mapping the correction type to an integer value

d = {'correctable':0,'uncorrectable':1,'z_lt_3':2,'nei_needs_redshift':4,'nei':3}

for k in d:
    newcol2[data['Correctable_Category'] == k] = d[k]

# Rearrange columns
newcols = fits.ColDefs([fits.Column(name='correction_type',format='I',array=newcol2)])
original_columns = data.columns
original_columns.del_col('Correctable_Category')
original_columns.change_name('OBJNO','survey_id')
original_columns.change_name('IMAGING','imaging')

# Write to a new file

hdu = fits.BinTableHDU.from_columns(original_columns[:3] + newcols + original_columns[3:])
hdu.writeto(outfile)
