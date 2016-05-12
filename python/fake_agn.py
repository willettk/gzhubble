# Fix the OBJNOs for the fake AGN

from astropy.io import ascii
from astropy.table.column import Column
import numpy as np

path = "/Users/willettk/Astronomy/Research/GalaxyZoo/csv"
data = ascii.read("{0}/gz_hst_lookup.csv".format(path))

ind = np.array([True if x[:4] == "AHZ7" else False for x in data['zooniverse_id']])

fake_agn_v1 = data[ind]['name'] < 1e9
fake_agn_v2 = data[ind]['name'] > 1e9

inds1 = np.arange(len(data))[ind][fake_agn_v1]
inds2 = np.arange(len(data))[ind][fake_agn_v2]

lrat = Column(data=np.zeros(len(data)),dtype=float,name='lrat')
agn_color = Column(data=np.zeros(len(data)),dtype=int,name='agn_color')
fake_bool = Column(data=np.zeros(len(data)),dtype=bool,name='fake_agn')
agn_type = Column(data=np.zeros(len(data)),dtype=int,name='fake_version')

ldict = {0:0.0,
         1:0.2,
         2:1.0,
         3:5.0,
         4:10.0,
         5:50.0}

for i in inds1[::15]:
    for c in range(3):
        for b in range(5):
            ind15 = i + c*5 + b
            data['name'][ind15] = data['name'][ind15]*100 + (c+1)*10 + (b+1)
            lrat[ind15] = ldict[b+1]
            agn_color[ind15] = c+1
            fake_bool[ind15] = True
            agn_type[ind15] = 1

for i in inds2[::16]:
    lrat[i] = ldict[0]
    agn_color[i] = 0
    fake_bool[i] = True
    agn_type[i] = 2
    for c in range(3):
        for b in range(5):
            ind15 = i + c*5 + b + 1
            lrat[ind15] = ldict[b+1]
            agn_color[ind15] = c+1
            fake_bool[ind15] = True
            agn_type[ind15] = 2

data.add_column(fake_bool)
data.add_column(agn_type)
data.add_column(lrat)
data.add_column(agn_color)

# Make some extra columns to ease matching later

data.write("{0}/gz_hst_lookup_fixed.csv".format(path))
