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
unmassed = newArray1[10:] + newArray2[10:] + newArray3[10:] + newArray4[10:]
unmassed = sorted(unmassed, key = operator.itemgetter(0))

lmass1 = []
for row in massed1:
	lmass1.append(row[6])
lunmass = []
for row in unmassed:
	lunmass.append(row[6])
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

idx_mass1, = np.where((z_1["id"] == 18935) | (z_1["id"] == 3740) | (z_1["id"] == 37194) | (z_1["id"] == 33990) | (z_1["id"] == 7002) | (z_1["id"] == 11221) | (z_1["id"] == 25242) | (z_1["id"] == 10509) | (z_1["id"] == 33712) | (z_1["id"] == 6449))
z_1mass = z_1[idx_mass1]
zunmass = data_z_flagged[(data_z_flagged["id"] != 18935) & (data_z_flagged["id"] != 3740) & (data_z_flagged["id"] != 37194) & (data_z_flagged["id"] != 33990) & (data_z_flagged["id"] != 7002) & (data_z_flagged["id"] != 11221) & (data_z_flagged["id"] != 25242) & (data_z_flagged["id"] != 10509) & (data_z_flagged["id"] != 33712) & (data_z_flagged["id"] != 6449) & (data_z_flagged["id"] != 25354) & (data_z_flagged["id"] != 16238) & (data_z_flagged["id"] != 9849) & (data_z_flagged["id"] != 26850) & (data_z_flagged["id"] != 6424) & (data_z_flagged["id"] != 17527) & (data_z_flagged["id"] != 22126) & (data_z_flagged["id"] != 18045) & (data_z_flagged["id"] != 30421) & (data_z_flagged["id"] != 38065) & (data_z_flagged["id"] != 2918) & (data_z_flagged["id"] != 38187) & (data_z_flagged["id"] != 24333) & (data_z_flagged["id"] != 18922) & (data_z_flagged["id"] != 33863) & (data_z_flagged["id"] != 31169) & (data_z_flagged["id"] != 23234) & (data_z_flagged["id"] != 15709) & (data_z_flagged["id"] != 28310) & (data_z_flagged["id"] != 15088) & (data_z_flagged["id"] != 179) & (data_z_flagged["id"] != 8790) & (data_z_flagged["id"] != 41033) & (data_z_flagged["id"] != 21306) & (data_z_flagged["id"] != 3048) & (data_z_flagged["id"] != 1849) & (data_z_flagged["id"] != 27584) & (data_z_flagged["id"] != 19028) & (data_z_flagged["id"] != 38731) & (data_z_flagged["id"] != 37799)]
idx_mass2, = np.where((z_2["id"] == 25354) | (z_2["id"] == 16238) | (z_2["id"] == 9849) | (z_2["id"] == 26850) | (z_2["id"] == 6424) | (z_2["id"] == 17527) | (z_2["id"] == 22126) | (z_2["id"] == 18045) | (z_2["id"] == 30421) | (z_2["id"] == 38065))
z_2mass = z_2[idx_mass2]
idx_mass3, = np.where((z_3["id"] == 2918) | (z_3["id"] == 38187) | (z_3["id"] == 24333) | (z_3["id"] == 18922) | (z_3["id"] == 33863) | (z_3["id"] == 31169) | (z_3["id"] == 23234) | (z_3["id"] == 15709) | (z_3["id"] == 28310) | (z_3["id"] == 15088))
z_3mass = z_3[idx_mass3]
idx_mass4, = np.where((z_4["id"] == 179) | (z_4["id"] == 8790) | (z_4["id"] == 41033) | (z_4["id"] == 21306) | (z_4["id"] == 3048) | (z_4["id"] == 1849) | (z_4["id"] == 27584) | (z_4["id"] == 19028) | (z_4["id"] == 38731) | (z_4["id"] == 37799))
z_4mass = z_4[idx_mass4]

z_1massed = z_1mass["z"]
zunmassed = zunmass["z"]
z_2massed = z_2mass["z"]
z_3massed = z_3mass["z"]
z_4massed = z_4mass["z"]


pylab.scatter(z_1massed, lmass1, color="r", label="chunk1 points, z < 0.5", edgecolor="k")
pylab.scatter(z_2massed, lmass2, color="b", label="chunk2 points, 0.5 < z < 1.5", edgecolor="k")
pylab.scatter(z_3massed, lmass3, color="g", label="chunk3 points, 1.5 < z < 2.5", edgecolor="k")
pylab.scatter(z_4massed, lmass4, color="k", label="chunk4 points, z > 2.5")
pylab.scatter(zunmassed, lunmass, color="0.7", alpha=0.5, edgecolor="none", label="unused points")
pylab.xlabel("z")
pylab.ylabel("log mass")
pylab.title("Aegis data")
pylab.vlines(0.5, 3, 13, linestyle="dashed")
pylab.vlines(1.5, 3, 13, linestyle="dashed")
pylab.vlines(2.5, 3, 13, linestyle="dashed")
pylab.legend(loc=4)
pylab.xlim([-0.5,6])
pylab.ylim([3.5,13])
pylab.savefig("z_mass")
pylab.ion()
pylab.show()

print "end"