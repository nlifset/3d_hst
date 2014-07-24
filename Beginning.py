print "start"
import os
os.chdir("/Volumes/TOSHIBA EXT/3d_hst")
from astropy.io import ascii
import numpy as np
from astropy.io import fits
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

idx_mass1, = np.where((z_1["id"] == 37587) | (z_1["id"] == 7013) | (z_1["id"] == 21796) | (z_1["id"] == 21306) | (z_1["id"] == 27258) | (z_1["id"] == 24280) | (z_1["id"] == 21160) | (z_1["id"] == 34823) | (z_1["id"] == 7631) | (z_1["id"] == 28216))
z_1mass = z_1[idx_mass1]
zunmass = data_z_flagged[(data_z_flagged["id"] != 37587) & (data_z_flagged["id"] != 7013) & (data_z_flagged["id"] != 21796) & (data_z_flagged["id"] != 21306) & (data_z_flagged["id"] != 27258) & (data_z_flagged["id"] != 24280) & (data_z_flagged["id"] != 21160) & (data_z_flagged["id"] != 34823) & (data_z_flagged["id"] != 7631) & (data_z_flagged["id"] != 28216) & (data_z_flagged["id"] != 9623) & (data_z_flagged["id"] != 15872) & (data_z_flagged["id"] != 21700) & (data_z_flagged["id"] != 12343) & (data_z_flagged["id"] != 2302) & (data_z_flagged["id"] != 25340) & (data_z_flagged["id"] != 2151) & (data_z_flagged["id"] != 25813) & (data_z_flagged["id"] != 13530) & (data_z_flagged["id"] != 24107) & (data_z_flagged["id"] != 10657) & (data_z_flagged["id"] != 19913) & (data_z_flagged["id"] != 5346) & (data_z_flagged["id"] != 32842) & (data_z_flagged["id"] != 338) & (data_z_flagged["id"] != 36988) & (data_z_flagged["id"] != 25942) & (data_z_flagged["id"] != 4117) & (data_z_flagged["id"] != 5371) & (data_z_flagged["id"] != 6789) & (data_z_flagged["id"] != 7338) & (data_z_flagged["id"] != 679) & (data_z_flagged["id"] != 25670) & (data_z_flagged["id"] != 3042) & (data_z_flagged["id"] != 11444) & (data_z_flagged["id"] != 37143) & (data_z_flagged["id"] != 36086) & (data_z_flagged["id"] != 26436) & (data_z_flagged["id"] != 29400) & (data_z_flagged["id"] != 37632)]
idx_mass2, = np.where((z_2["id"] == 9623) | (z_2["id"] == 15872) | (z_2["id"] == 21700) | (z_2["id"] == 12343) | (z_2["id"] == 2302) | (z_2["id"] == 25340) | (z_2["id"] == 2151) | (z_2["id"] == 25813) | (z_2["id"] == 13530) | (z_2["id"] == 24107))
z_2mass = z_2[idx_mass2]
idx_mass3, = np.where((z_3["id"] == 10657) | (z_3["id"] == 19913) | (z_3["id"] == 5346) | (z_3["id"] == 32842) | (z_3["id"] == 338) | (z_3["id"] == 36988) | (z_3["id"] == 25942) | (z_3["id"] == 4117) | (z_3["id"] == 5371) | (z_3["id"] == 6789))
z_3mass = z_3[idx_mass3]
idx_mass4, = np.where((z_4["id"] == 7338) | (z_4["id"] == 679) | (z_4["id"] == 25670) | (z_4["id"] == 3042) | (z_4["id"] == 11444) | (z_4["id"] == 37143) | (z_4["id"] == 36086) | (z_4["id"] == 26436) | (z_4["id"] == 29400) | (z_4["id"] == 37632))
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
pylab.title("Goods-N data")
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