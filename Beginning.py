print "start"
import os
os.chdir("/Users/noah/Documents/3d_hst")
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
	
z_1 = data_z_flagged[(data_z_flagged["z"] < 0.5)]
z_2 = data_z_flagged[(data_z_flagged["z"] >= 0.5) & (data_z_flagged["z"] < 1.5)]
z_3 = data_z_flagged[(data_z_flagged["z"] >= 1.5) & (data_z_flagged["z"] < 2.5)]
z_4 = data_z_flagged[(data_z_flagged["z"] >= 2.5)]

idx_mass1, = np.where((z_1["id"] == 18935) | (z_1["id"] == 3740) | (z_1["id"] == 37194) | (z_1["id"] == 33990) | (z_1["id"] == 7002) | (z_1["id"] == 11221) | (z_1["id"] == 25242) | (z_1["id"] == 10509) | (z_1["id"] == 33712) | (z_1["id"] == 8664))
z_1mass = z_1[idx_mass1]
z_1unmass = z_1[(z_1["id"] != 18935) & (z_1["id"] != 3740) & (z_1["id"] != 37194) & (z_1["id"] != 33990) & (z_1["id"] != 7002) & (z_1["id"] != 11221) & (z_1["id"] != 25242) & (z_1["id"] != 10509) & (z_1["id"] != 33712) & (z_1["id"] != 8664)]
idx_mass2, = np.where((z_2["id"] == 25354) | (z_2["id"] == 16238) | (z_2["id"] == 9849) | (z_2["id"] == 26850) | (z_2["id"] == 6424) | (z_2["id"] == 17527) | (z_2["id"] == 22126) | (z_2["id"] == 18045) | (z_2["id"] == 30421) | (z_2["id"] == 38065))
z_2mass = z_2[idx_mass2]
z_2unmass = z_2[(z_2["id"] != 25354) & (z_2["id"] != 16238) & (z_2["id"] != 9849) & (z_2["id"] != 26850) & (z_2["id"] != 6424) & (z_2["id"] != 17527) & (z_2["id"] != 22126) & (z_2["id"] != 18045) & (z_2["id"] != 30421) & (z_2["id"] != 38065)]
idx_mass3, = np.where((z_3["id"] == 2918) | (z_3["id"] == 38187) | (z_3["id"] == 24333) | (z_3["id"] == 18922) | (z_3["id"] == 33863) | (z_3["id"] == 31169) | (z_3["id"] == 23234) | (z_3["id"] == 15709) | (z_3["id"] == 28310) | (z_3["id"] == 15088))
z_3mass = z_3[idx_mass3]
z_3unmass = z_3[(z_3["id"] != 2918) & (z_3["id"] != 38187) & (z_3["id"] != 24333) & (z_3["id"] != 18922) & (z_3["id"] != 33863) & (z_3["id"] != 31169) & (z_3["id"] != 23234) & (z_3["id"] != 15709) & (z_3["id"] != 28310) & (z_3["id"] != 15088)]
idx_mass4, = np.where((z_4["id"] == 179) | (z_4["id"] == 8790) | (z_4["id"] == 41033) | (z_4["id"] == 21306) | (z_4["id"] == 3048) | (z_4["id"] == 1849) | (z_4["id"] == 27584) | (z_4["id"] == 19028) | (z_4["id"] == 38731) | (z_4["id"] == 37799))
z_4mass = z_4[idx_mass4]
z_4unmass = z_4[(z_4["id"] != 179) & (z_4["id"] != 8790) & (z_4["id"] != 41033) & (z_4["id"] != 21306) & (z_4["id"] != 3048) & (z_4["id"] != 1849) & (z_4["id"] != 27584) & (z_4["id"] != 19028) & (z_4["id"] != 38731) & (z_4["id"] != 37799)]

z_1massed = z_1mass["z"]
z_1unmassed = z_1unmass["z"]
z_2massed = z_2mass["z"]
z_2unmassed = z_2unmass["z"]
z_3massed = z_3mass["z"]
z_3unmassed = z_3unmass["z"]
z_4massed = z_4mass["z"]
z_4unmassed = z_4unmass["z"]


pylab.scatter(z_1massed, lmass1, color="r", label="chunk1 points")
pylab.scatter(z_1unmassed, lunmass1, color="g", alpha=0.4, edgecolor="none", label="unused points")
pylab.scatter(z_2massed, lmass2, color="b", label="chunk2 points")
pylab.scatter(z_2unmassed, lunmass2, color="g", alpha=0.4, edgecolor="none", label="_nolegend_")
pylab.scatter(z_3massed, lmass3, color="c", label="chunk3 points")
pylab.scatter(z_3unmassed, lunmass3, color="g", alpha=0.4, edgecolor="none", label="_nolegend_")
pylab.scatter(z_4massed, lmass4, color="k", label="chunk4 points")
pylab.scatter(z_4unmassed, lunmass4, color="g", alpha=0.4, edgecolor="none", label="_nolegend_")
pylab.xlabel("z")
pylab.ylabel("log mass")
pylab.title("Aegis data")
pylab.vlines(0.5, 3, 13, linestyle="dashed")
pylab.vlines(1.5, 3, 13, linestyle="dashed")
pylab.vlines(2.5, 3, 13, linestyle="dashed")
pylab.legend(loc=4)

pylab.ion()
pylab.show()
pylab.fig.savefig("z vs mass")
print "end"