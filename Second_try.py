print "start"
import os
os.chdir("/Volumes/TOSHIBA EXT/3d_hst")
from astropy.io import ascii
import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
import pylab
from operator import itemgetter, attrgetter
import astropy.coordinates as coord

data = ascii.read("aegis_3dhst.v4.1.cat")
data_fast = ascii.read("aegis_3dhst.v4.1.fout")

idx, = np.where((data["use_phot"] == 1.0) & (data["star_flag"] == 0.0))
data_fast_flag = data_fast[idx]
data_flag = data[idx]

idx2, = np.where(data_fast_flag["lsfr"] != -99)
data_fast_flagged = data_fast_flag[idx2]
data_flagged = data_flag[idx2]

chunk1 = data_fast_flagged[(data_fast_flagged["z"] < 0.5)]
chunk2 = data_fast_flagged[(data_fast_flagged["z"] >= 0.5) & (data_fast_flagged["z"] < 1.5)]
chunk3 = data_fast_flagged[(data_fast_flagged["z"] >= 1.5) & (data_fast_flagged["z"] < 2.5)]
chunk4 = data_fast_flagged[(data_fast_flagged["z"] >= 2.5)]

newArray1 = [(x[6], x) for x in chunk1]
newArray1.sort(key = lambda x : x[0], reverse=True)
newArray2 = [(x[6], x) for x in chunk2]
newArray2.sort(key = lambda x : x[0], reverse=True)
newArray3 = [(x[6], x) for x in chunk3]
newArray3.sort(key = lambda x : x[0], reverse=True)
newArray4 = [(x[6], x) for x in chunk4]
newArray4.sort(key = lambda x : x[0], reverse=True)

massed1 = newArray1[0:10]
massed2 = newArray2[0:10]
massed3 = newArray3[0:10]
massed4 = newArray4[0:10]

p = data_flagged[(data_flagged['id'] == 37194)]

print p[0]



print "end"
