print "start"
import os
os.chdir('/Users/noah/Documents/3d_hst')
from astropy.io import ascii
import numpy as np
from astropy.io import fits
import matplotlib.pylab as plt
import pylab

data = ascii.read("aegis_3dhst.v4.1.cat")
use_phot = data["use_phot"]

idx, = np.where(use_phot == 1.0)

data_flagged = data[idx]

data_fast = ascii.read("aegis_3dhst.v4.1.fout")
data_fast_flagged = data_fast[idx]

lmass = data_fast_flagged["lmass"]
lsfr = data_fast_flagged["lsfr"]
idx2, = np.where(lsfr > -99)
lsfr_ = lsfr[idx2]
lmass_ = lmass[idx2]

print pylab.scatter(lmass_, lsfr_)
"""the z value goes from -1 to at least 5"""
pylab.title("aegis galaxy data")
pylab.xlabel("log mass")
pylab.ylabel("log star formation rate")
pylab.show()
print "end"