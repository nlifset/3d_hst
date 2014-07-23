print "start"
import os
os.chdir('/Volumes/TOSHIBA EXT/3d_hst')
from astropy.io import ascii
import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
import pylab

data = ascii.read("aegis_3dhst.v4.1.cat")
data_fast = ascii.read("aegis_3dhst.v4.1.fout")

idx, = np.where((data["use_phot"] == 1.0) & (data["star_flag"] == 0.0))
data_fast_flag = data_fast[idx]

idx2, = np.where(data_fast_flag["lsfr"] != -99)
data_fast_flagged = data_fast_flag[idx2]

chunk1 = data_fast_flagged[(data_fast_flagged["z"] < 0.5)]
chunk2 = data_fast_flagged[(data_fast_flagged["z"] >= 0.5) & (data_fast_flagged["z"] < 1.5)]
chunk3 = data_fast_flagged[(data_fast_flagged["z"] >= 1.5) & (data_fast_flagged["z"] < 2.5)]
chunk4 = data_fast_flagged[(data_fast_flagged["z"] >= 2.5)]

lmass1 = chunk1["lmass"]
lsfr1 = chunk1["lsfr"]
lmass2 = chunk2["lmass"]
lsfr2 = chunk2["lsfr"]
lmass3 = chunk3["lmass"]
lsfr3 = chunk3["lsfr"]
lmass4 = chunk4["lmass"]
lsfr4 = chunk4["lsfr"]

fig,pylab.axes = pylab.subplots(2, 2, sharex = True, sharey = True)

a1 = pylab.axes[0,0]
a2 = pylab.axes[0,1]
a3 = pylab.axes[1,0]
a4 = pylab.axes[1,1]

a1.plot(lmass1, lsfr1, alpha=0.3, marker="o", markeredgecolor = "none", linestyle="none")
a1.set_title("z < 0.5", fontsize=12)
a2.plot(lmass2, lsfr2, alpha=0.3, marker="o", markeredgecolor = "none", linestyle="none")
a2.set_title("0.5< z <1.5", fontsize=12)
a3.plot(lmass3, lsfr3, alpha=0.3, marker="o", markeredgecolor = "none", linestyle="none")
a3.set_xlabel("1.5< z <2.5")
a4.plot(lmass4, lsfr4, alpha=0.3, marker="o", markeredgecolor = "none", linestyle="none")
a4.set_xlabel("z > 2.5")

pylab.suptitle("aegis data", fontsize=20)
fig.text(0.44, 0.01, "log mass", fontsize=18)
fig.text(0.01, 0.67, "log star formation rate", rotation = "vertical", fontsize=18)
pylab.xlim([7,13])
pylab.ylim([-10,5])
fig.subplots_adjust(wspace=0.05, hspace=0.05)
pylab.ion()
pylab.show()
fig.savefig("test.jpg")
print "end"