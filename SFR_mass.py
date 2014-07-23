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

data = ascii.read("aegis_3dhst.v4.1.cat")
data_fast = ascii.read("aegis_3dhst.v4.1.fout")
data_z = ascii.read("aegis_3dhst.v4.0.sfr")

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

idx_mass1, = np.where((z_1["id"] == 18935) | (z_1["id"] == 3740) | (z_1["id"] == 37194) | (z_1["id"] == 33990) | (z_1["id"] == 7002) | (z_1["id"] == 11221) | (z_1["id"] == 25242) | (z_1["id"] == 10509) | (z_1["id"] == 33712) | (z_1["id"] == 8664))
z_1mass = z_1[idx_mass1]
idx_mass2, = np.where((z_2["id"] == 25354) | (z_2["id"] == 16238) | (z_2["id"] == 9849) | (z_2["id"] == 26850) | (z_2["id"] == 6424) | (z_2["id"] == 17527) | (z_2["id"] == 22126) | (z_2["id"] == 18045) | (z_2["id"] == 30421) | (z_2["id"] == 38065))
z_2mass = z_2[idx_mass2]
idx_mass3, = np.where((z_3["id"] == 2918) | (z_3["id"] == 38187) | (z_3["id"] == 24333) | (z_3["id"] == 18922) | (z_3["id"] == 33863) | (z_3["id"] == 31169) | (z_3["id"] == 23234) | (z_3["id"] == 15709) | (z_3["id"] == 28310) | (z_3["id"] == 15088))
z_3mass = z_3[idx_mass3]
idx_mass4, = np.where((z_4["id"] == 179) | (z_4["id"] == 8790) | (z_4["id"] == 41033) | (z_4["id"] == 21306) | (z_4["id"] == 3048) | (z_4["id"] == 1849) | (z_4["id"] == 27584) | (z_4["id"] == 19028) | (z_4["id"] == 38731) | (z_4["id"] == 37799))
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