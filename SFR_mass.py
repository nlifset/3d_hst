print "start"
import os
os.chdir("/Volumes/TOSHIBA EXT/3d_hst")
from astropy.io import ascii
import numpy as np
from astropy.io import fits
import astropy.units as u
import astropy.constants
import matplotlib.pyplot as plt
import pylab
import operator 

data = ascii.read("goodsn_3dhst.v4.1.cat")
data_fast = ascii.read("goodsn_3dhst.v4.1.fout")
data_z = ascii.read("goodsn_3dhst.v4.0.sfr")

idx, = np.where((data["use_phot"] == 1.0) & (data["star_flag"] == 0.0))
data_fast_flag = data_fast[idx]
data_flag = data[idx]
data_z_flag = data_z[idx]

idx2, = np.where(data_fast_flag["lsfr"] != -99)
data_fast_flagged = data_fast_flag[idx2]
data_flagged = data_flag[idx2]
data_z_flagged = data_z_flag[idx2]

chunk1 = data_fast_flagged[(data_z_flagged["z"] < 0.5)]
chunk2 = data_fast_flagged[(data_z_flagged["z"] >= 0.5) & (data_z_flagged["z"] < 1.5)]
chunk3 = data_fast_flagged[(data_z_flagged["z"] >= 1.5) & (data_z_flagged["z"] < 2.5)]
chunk4 = data_fast_flagged[(data_z_flagged["z"] >= 2.5)]

newArray1 = sorted(chunk1, key = operator.itemgetter(6), reverse=True)
newArray2 = sorted(chunk2, key = operator.itemgetter(6), reverse=True)
newArray3 = sorted(chunk3, key = operator.itemgetter(6), reverse=True)
newArray4 = sorted(chunk4, key = operator.itemgetter(6), reverse=True)

massed1 = newArray1[:10]
massed2 = newArray2[:10]
massed3 = newArray3[:10]
massed4 = newArray4[:10]
massed1 = sorted(massed1, key = operator.itemgetter(0))
massed2 = sorted(massed2, key = operator.itemgetter(0))
massed3 = sorted(massed3, key = operator.itemgetter(0))
massed4 = sorted(massed4, key = operator.itemgetter(0))

lmass1 = []
for row in massed1:
	lmass1.append(row[6])
lmass2 = []
for row in massed2:
	lmass2.append(row[6])
lmass3 = []
for row in massed3:
	lmass3.append(row[6])
lmass4 = []
for row in massed4:
	lmass4.append(row[6])

	
z_1 = data_z_flagged[(data_z_flagged["z"] < 0.5)]
z_2 = data_z_flagged[(data_z_flagged["z"] >= 0.5) & (data_z_flagged["z"] < 1.5)]
z_3 = data_z_flagged[(data_z_flagged["z"] >= 1.5) & (data_z_flagged["z"] < 2.5)]
z_4 = data_z_flagged[(data_z_flagged["z"] >= 2.5)]

idx_mass1, = np.where((z_1["id"] == 37587) | (z_1["id"] == 7013) | (z_1["id"] == 21796) | (z_1["id"] == 21306) | (z_1["id"] == 27258) | (z_1["id"] == 24280) | (z_1["id"] == 21160) | (z_1["id"] == 34823) | (z_1["id"] == 7631) | (z_1["id"] == 28216))
z_1mass = z_1[idx_mass1]
idx_mass2, = np.where((z_2["id"] == 9623) | (z_2["id"] == 15872) | (z_2["id"] == 21700) | (z_2["id"] == 12343) | (z_2["id"] == 2302) | (z_2["id"] == 25340) | (z_2["id"] == 2151) | (z_2["id"] == 25813) | (z_2["id"] == 13530) | (z_2["id"] == 24107))
z_2mass = z_2[idx_mass2]
idx_mass3, = np.where((z_3["id"] == 10657) | (z_3["id"] == 19913) | (z_3["id"] == 5346) | (z_3["id"] == 32842) | (z_3["id"] == 338) | (z_3["id"] == 36988) | (z_3["id"] == 25942) | (z_3["id"] == 4117) | (z_3["id"] == 5371) | (z_3["id"] == 6789))
z_3mass = z_3[idx_mass3]
idx_mass4, = np.where((z_4["id"] == 7338) | (z_4["id"] == 679) | (z_4["id"] == 25670) | (z_4["id"] == 3042) | (z_4["id"] == 11444) | (z_4["id"] == 37143) | (z_4["id"] == 36086) | (z_4["id"] == 26436) | (z_4["id"] == 29400) | (z_4["id"] == 37632))
z_4mass = z_4[idx_mass4]

sfr1 = z_1mass["sfr"]
sfr2 = z_2mass["sfr"]
sfr3 = z_3mass["sfr"]
sfr4 = z_4mass["sfr"]


pylab.plot(sfr1, lmass1, marker="o", markeredgecolor = "k", linestyle="none", color="r", label="chunk1, z < 0.5")
pylab.plot(sfr2, lmass2, marker="o", markeredgecolor = "k", linestyle="none", color="g", label="chunk2, 0.5 < z < 1.5")
pylab.plot(sfr3, lmass3, marker="o", markeredgecolor = "k", linestyle="none", color="b", label="chunk3, 1.5 < z < 2.5")
pylab.plot(sfr4, lmass4, marker="o", markeredgecolor = "none", linestyle="none", color="k", label="chunk4, z > 2.5")

pylab.suptitle("10 most massive stars of 4 redshift chunks", fontsize=18)
pylab.xlabel("star formation rate", fontsize=16)
pylab.ylabel("log mass", fontsize=16)
pylab.legend(loc=2)
pylab.xscale("log")
pylab.ion()
pylab.show()
pylab.savefig("sfr_mass")

print "end"