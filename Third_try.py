print "start"
import os
os.chdir("/Volumes/TOSHIBA EXT/3d_hst")
from astropy.io import ascii
import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
import pylab
import operator 

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

newArray1 = sorted(chunk1, key = operator.itemgetter(6), reverse=True)
newArray2 = sorted(chunk2, key = operator.itemgetter(6), reverse=True)
newArray3 = sorted(chunk3, key = operator.itemgetter(6), reverse=True)
newArray4 = sorted(chunk4, key = operator.itemgetter(6), reverse=True)

massed1 = newArray1[:10]
massed2 = newArray2[:10]
massed3 = newArray3[:10]
massed4 = newArray4[:10]
unmassed1 = newArray1[10:]
unmassed2 = newArray2[10:]
unmassed3 = newArray3[10:]
unmassed4 = newArray4[10:]

lmass1 = []
for row in massed1:
	lmass1.append(row[6])
lunmass1 = []
for row in unmassed1:
	lunmass1.append(row[6])
lsfr1 = []
for row in massed1:
	lsfr1.append(row[7])
lunsfr1 =[]
for row in unmassed1:
	lunsfr1.append(row[7])
lmass2 = []
for row in massed2:
	lmass2.append(row[6])
lunmass2 = []
for row in unmassed2:
	lunmass2.append(row[6])
lsfr2 = []
for row in massed2:
	lsfr2.append(row[7])
lunsfr2 =[]
for row in unmassed2:
	lunsfr2.append(row[7])
lmass3 = []
for row in massed3:
	lmass3.append(row[6])
lunmass3 = []
for row in unmassed3:
	lunmass3.append(row[6])
lsfr3 = []
for row in massed3:
	lsfr3.append(row[7])
lunsfr3 =[]
for row in unmassed3:
	lunsfr3.append(row[7])
lmass4 = []
for row in massed4:
	lmass4.append(row[6])
lunmass4 = []
for row in unmassed4:
	lunmass4.append(row[6])
lsfr4 = []
for row in massed4:
	lsfr4.append(row[7])
lunsfr4 =[]
for row in unmassed4:
	lunsfr4.append(row[7])
	
fig,pylab.axes = pylab.subplots(2, 2, sharex = True, sharey = True)

a1 = pylab.axes[0,0]
a2 = pylab.axes[0,1]
a3 = pylab.axes[1,0]
a4 = pylab.axes[1,1]

a1.plot(lmass1, lsfr1, alpha=0.3, marker="o", markeredgecolor = "none", linestyle="none", color="r")
a1.plot(lunmass1, lunsfr1, alpha=0.3, marker="o", markeredgecolor = "none", linestyle="none", color="g")
a1.set_title("z < 0.5", fontsize=12)
a2.plot(lmass2, lsfr2, alpha=0.3, marker="o", markeredgecolor = "none", linestyle="none", color="r")
a2.plot(lunmass2, lunsfr2, alpha=0.3, marker="o", markeredgecolor = "none", linestyle="none", color="g")
a2.set_title("0.5< z <1.5", fontsize=12)
a3.plot(lmass3, lsfr3, alpha=0.3, marker="o", markeredgecolor = "none", linestyle="none", color="r")
a3.plot(lunmass3, lunsfr3, alpha=0.3, marker="o", markeredgecolor = "none", linestyle="none", color="g")
a3.set_xlabel("1.5< z <2.5")
a4.plot(lmass4, lsfr4, alpha=0.3, marker="o", markeredgecolor = "none", linestyle="none", color="r")
a4.plot(lunmass4, lunsfr4, alpha=0.3, marker="o", markeredgecolor = "none", linestyle="none", color="g")
a4.set_xlabel("z > 2.5")

pylab.suptitle("aegis data", fontsize=20)
fig.text(0.44, 0.01, "log mass", fontsize=18)
fig.text(0.01, 0.67, "log star formation rate", rotation = "vertical", fontsize=18)
pylab.xlim([7,13])
pylab.ylim([-10,5])
fig.subplots_adjust(wspace=0.05, hspace=0.05)
pylab.ion()
pylab.show()


print "end"