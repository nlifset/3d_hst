print "start"
import os
os.chdir("/Volumes/TOSHIBA EXT/3d_hst")
from astropy.io import ascii
import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
import pylab
import operator 
import math
np.set_printoptions(threshold="inf")
def z_mass():
	
	data1 = ascii.read("aegis_3dhst.v4.1.cat")
	data_fast1 = ascii.read("aegis_3dhst.v4.1.fout")
	data_z1 = ascii.read("aegis_3dhst.v4.0.sfr")
	
	idx1, = np.where((data1["use_phot"] == 1.0) & (data1["star_flag"] == 0.0))
	data_fast_flag1 = data_fast1[idx1]
	data_flag1 = data1[idx1]
	data_z_flag1 = data_z1[idx1]
	
	idx_1, = np.where(data_fast_flag1["lsfr"] != -99)
	data_fast_flagged1 = data_fast_flag1[idx_1]
	data_flagged1 = data_flag1[idx_1]
	data_z_flagged1 = data_z_flag1[idx_1]
	
	chunk1_a = data_fast_flagged1[(data_z_flagged1["z"] < 1)]
	chunk2_a = data_fast_flagged1[(data_z_flagged1["z"] >= 1) & (data_z_flagged1["z"] < 1.5)]
	chunk3_a = data_fast_flagged1[(data_z_flagged1["z"] >= 1.5) & (data_z_flagged1["z"] < 2)]
	chunk4_a = data_fast_flagged1[(data_z_flagged1["z"] >= 2) & (data_z_flagged1["z"] <2.5)]
	chunk5_a = data_fast_flagged1[(data_z_flagged1["z"] >= 2.5)]
	
	massed1_a = chunk1_a[(chunk1_a["lmass"] >= 11)]
	massed2_a = chunk2_a[(chunk2_a["lmass"] >= 11)]
	massed3_a = chunk3_a[(chunk3_a["lmass"] >= 11)]
	massed4_a = chunk4_a[(chunk4_a["lmass"] >= 11)]
	massed5_a = chunk5_a[(chunk5_a["lmass"] >= 11)]
	unmassed_a = data_fast_flagged1[(data_fast_flagged1["lmass"] < 11)]
	
	lmass1_a = massed1_a["lmass"]
	lmass2_a = massed2_a["lmass"]
	lmass3_a = massed3_a["lmass"]
	lmass4_a = massed4_a["lmass"]
	lmass5_a = massed5_a["lmass"]
	lunmass_a = unmassed_a["lmass"]
	
	chunk_1_a = data_z_flagged1[(data_z_flagged1["z"] < 1)]
	chunk_2_a = data_z_flagged1[(data_z_flagged1["z"] >= 1) & (data_z_flagged1["z"] < 1.5)]
	chunk_3_a = data_z_flagged1[(data_z_flagged1["z"] >= 1.5) & (data_z_flagged1["z"] < 2)]
	chunk_4_a = data_z_flagged1[(data_z_flagged1["z"] >= 2) & (data_z_flagged1["z"] <2.5)]
	chunk_5_a = data_z_flagged1[(data_z_flagged1["z"] >= 2.5)]
	
	massed_1_a = chunk_1_a[(chunk1_a["lmass"] >= 11)]
	massed_2_a = chunk_2_a[(chunk2_a["lmass"] >= 11)]
	massed_3_a = chunk_3_a[(chunk3_a["lmass"] >= 11)]
	massed_4_a = chunk_4_a[(chunk4_a["lmass"] >= 11)]
	massed_5_a = chunk_5_a[(chunk5_a["lmass"] >= 11)]
	unmassed_a_ = data_z_flagged1[(data_fast_flagged1["lmass"] < 11)]
	
	z1_a = massed_1_a["z"]
	z2_a = massed_2_a["z"]
	z3_a = massed_3_a["z"]
	z4_a = massed_4_a["z"]
	z5_a = massed_5_a["z"]
	unz_a = unmassed_a_["z"]
	
	
	
	data2 = ascii.read("cosmos_3dhst.v4.1.cat")
	data_fast2 = ascii.read("cosmos_3dhst.v4.1.fout")
	data_z2 = ascii.read("cosmos_3dhst.v4.0.sfr")
	
	idx2, = np.where((data2["use_phot"] == 1.0) & (data2["star_flag"] == 0.0))
	data_fast_flag2 = data_fast2[idx2]
	data_flag2 = data2[idx2]
	data_z_flag2 = data_z2[idx2]
	
	idx_2, = np.where(data_fast_flag2["lsfr"] != -99)
	data_fast_flagged2 = data_fast_flag2[idx_2]
	data_flagged2 = data_flag2[idx_2]
	data_z_flagged2 = data_z_flag2[idx_2]
	
	chunk1_c = data_fast_flagged2[(data_z_flagged2["z"] < 1)]
	chunk2_c = data_fast_flagged2[(data_z_flagged2["z"] >= 1) & (data_z_flagged2["z"] < 1.5)]
	chunk3_c = data_fast_flagged2[(data_z_flagged2["z"] >= 1.5) & (data_z_flagged2["z"] < 2)]
	chunk4_c = data_fast_flagged2[(data_z_flagged2["z"] >= 2) & (data_z_flagged2["z"] <2.5)]
	chunk5_c = data_fast_flagged2[(data_z_flagged2["z"] >= 2.5)]
	
	massed1_c = chunk1_c[(chunk1_c["lmass"] >= 11)]
	massed2_c = chunk2_c[(chunk2_c["lmass"] >= 11)]
	massed3_c = chunk3_c[(chunk3_c["lmass"] >= 11)]
	massed4_c = chunk4_c[(chunk4_c["lmass"] >= 11)]
	massed5_c = chunk5_c[(chunk5_c["lmass"] >= 11)]
	unmassed_c = data_fast_flagged2[(data_fast_flagged2["lmass"] < 11)]
	
	lmass1_c = massed1_c["lmass"]
	lmass2_c = massed2_c["lmass"]
	lmass3_c = massed3_c["lmass"]
	lmass4_c = massed4_c["lmass"]
	lmass5_c = massed5_c["lmass"]
	lunmass_c = unmassed_c["lmass"]
	
	chunk_1_c = data_z_flagged2[(data_z_flagged2["z"] < 1)]
	chunk_2_c = data_z_flagged2[(data_z_flagged2["z"] >= 1) & (data_z_flagged2["z"] < 1.5)]
	chunk_3_c = data_z_flagged2[(data_z_flagged2["z"] >= 1.5) & (data_z_flagged2["z"] < 2)]
	chunk_4_c = data_z_flagged2[(data_z_flagged2["z"] >= 2) & (data_z_flagged2["z"] <2.5)]
	chunk_5_c = data_z_flagged2[(data_z_flagged2["z"] >= 2.5)]
	
	massed_1_c = chunk_1_c[(chunk1_c["lmass"] >= 11)]
	massed_2_c = chunk_2_c[(chunk2_c["lmass"] >= 11)]
	massed_3_c = chunk_3_c[(chunk3_c["lmass"] >= 11)]
	massed_4_c = chunk_4_c[(chunk4_c["lmass"] >= 11)]
	massed_5_c = chunk_5_c[(chunk5_c["lmass"] >= 11)]
	unmassed_c_ = data_z_flagged2[(data_fast_flagged2["lmass"] < 11)]
	
	z1_c = massed_1_c["z"]
	z2_c = massed_2_c["z"]
	z3_c = massed_3_c["z"]
	z4_c = massed_4_c["z"]
	z5_c = massed_5_c["z"]
	unz_c = unmassed_c_["z"]
	
	
	
	
	data3 = ascii.read("goodsn_3dhst.v4.1.cat")
	data_fast3 = ascii.read("goodsn_3dhst.v4.1.fout")
	data_z3 = ascii.read("goodsn_3dhst.v4.0.sfr")
	
	idx3, = np.where((data3["use_phot"] == 1.0) & (data3["star_flag"] == 0.0))
	data_fast_flag3 = data_fast3[idx3]
	data_flag3 = data3[idx3]
	data_z_flag3 = data_z3[idx3]
	
	idx_3, = np.where(data_fast_flag3["lsfr"] != -99)
	data_fast_flagged3 = data_fast_flag3[idx_3]
	data_flagged3 = data_flag3[idx_3]
	data_z_flagged3 = data_z_flag3[idx_3]
	
	chunk1_n = data_fast_flagged3[(data_z_flagged3["z"] < 1)]
	chunk2_n = data_fast_flagged3[(data_z_flagged3["z"] >= 1) & (data_z_flagged3["z"] < 1.5)]
	chunk3_n = data_fast_flagged3[(data_z_flagged3["z"] >= 1.5) & (data_z_flagged3["z"] < 2)]
	chunk4_n = data_fast_flagged3[(data_z_flagged3["z"] >= 2) & (data_z_flagged3["z"] <2.5)]
	chunk5_n = data_fast_flagged3[(data_z_flagged3["z"] >= 2.5)]
	
	massed1_n = chunk1_n[(chunk1_n["lmass"] >= 11)]
	massed2_n = chunk2_n[(chunk2_n["lmass"] >= 11)]
	massed3_n = chunk3_n[(chunk3_n["lmass"] >= 11)]
	massed4_n = chunk4_n[(chunk4_n["lmass"] >= 11)]
	massed5_n = chunk5_n[(chunk5_n["lmass"] >= 11)]
	unmassed_n = data_fast_flagged3[(data_fast_flagged3["lmass"] < 11)]
	
	lmass1_n = massed1_n["lmass"]
	lmass2_n = massed2_n["lmass"]
	lmass3_n = massed3_n["lmass"]
	lmass4_n = massed4_n["lmass"]
	lmass5_n = massed5_n["lmass"]
	lunmass_n = unmassed_n["lmass"]
	
	chunk_1_n = data_z_flagged3[(data_z_flagged3["z"] < 1)]
	chunk_2_n = data_z_flagged3[(data_z_flagged3["z"] >= 1) & (data_z_flagged3["z"] < 1.5)]
	chunk_3_n = data_z_flagged3[(data_z_flagged3["z"] >= 1.5) & (data_z_flagged3["z"] < 2)]
	chunk_4_n = data_z_flagged3[(data_z_flagged3["z"] >= 2) & (data_z_flagged3["z"] <2.5)]
	chunk_5_n = data_z_flagged3[(data_z_flagged3["z"] >= 2.5)]
	
	massed_1_n = chunk_1_n[(chunk1_n["lmass"] >= 11)]
	massed_2_n = chunk_2_n[(chunk2_n["lmass"] >= 11)]
	massed_3_n = chunk_3_n[(chunk3_n["lmass"] >= 11)]
	massed_4_n = chunk_4_n[(chunk4_n["lmass"] >= 11)]
	massed_5_n = chunk_5_n[(chunk5_n["lmass"] >= 11)]
	unmassed_n_ = data_z_flagged3[(data_fast_flagged3["lmass"] < 11)]
	
	z1_n = massed_1_n["z"]
	z2_n = massed_2_n["z"]
	z3_n = massed_3_n["z"]
	z4_n = massed_4_n["z"]
	z5_n = massed_5_n["z"]
	unz_n = unmassed_n_["z"]
	
	
	
	
	data4 = ascii.read("goodss_3dhst.v4.1.cat")
	data_fast4 = ascii.read("goodss_3dhst.v4.1.fout")
	data_z4 = ascii.read("goodss_3dhst.v4.0.sfr")
	
	idx4, = np.where((data4["use_phot"] == 1.0) & (data4["star_flag"] == 0.0))
	data_fast_flag4 = data_fast4[idx4]
	data_flag4 = data4[idx4]
	data_z_flag4 = data_z4[idx4]
	
	idx_4, = np.where(data_fast_flag4["lsfr"] != -99)
	data_fast_flagged4 = data_fast_flag4[idx_4]
	data_flagged4 = data_flag4[idx_4]
	data_z_flagged4 = data_z_flag4[idx_4]
	
	chunk1_s = data_fast_flagged4[(data_z_flagged4["z"] < 1)]
	chunk2_s = data_fast_flagged4[(data_z_flagged4["z"] >= 1) & (data_z_flagged4["z"] < 1.5)]
	chunk3_s = data_fast_flagged4[(data_z_flagged4["z"] >= 1.5) & (data_z_flagged4["z"] < 2)]
	chunk4_s = data_fast_flagged4[(data_z_flagged4["z"] >= 2) & (data_z_flagged4["z"] <2.5)]
	chunk5_s = data_fast_flagged4[(data_z_flagged4["z"] >= 2.5)]
	
	massed1_s = chunk1_s[(chunk1_s["lmass"] >= 11)]
	massed2_s = chunk2_s[(chunk2_s["lmass"] >= 11)]
	massed3_s = chunk3_s[(chunk3_s["lmass"] >= 11)]
	massed4_s = chunk4_s[(chunk4_s["lmass"] >= 11)]
	massed5_s = chunk5_s[(chunk5_s["lmass"] >= 11)]
	unmassed_s = data_fast_flagged4[(data_fast_flagged4["lmass"] < 11)]
	
	lmass1_s = massed1_s["lmass"]
	lmass2_s = massed2_s["lmass"]
	lmass3_s = massed3_s["lmass"]
	lmass4_s = massed4_s["lmass"]
	lmass5_s = massed5_s["lmass"]
	lunmass_s = unmassed_s["lmass"]
	
	chunk_1_s = data_z_flagged4[(data_z_flagged4["z"] < 1)]
	chunk_2_s = data_z_flagged4[(data_z_flagged4["z"] >= 1) & (data_z_flagged4["z"] < 1.5)]
	chunk_3_s = data_z_flagged4[(data_z_flagged4["z"] >= 1.5) & (data_z_flagged4["z"] < 2)]
	chunk_4_s = data_z_flagged4[(data_z_flagged4["z"] >= 2) & (data_z_flagged4["z"] <2.5)]
	chunk_5_s = data_z_flagged4[(data_z_flagged4["z"] >= 2.5)]
	
	massed_1_s = chunk_1_s[(chunk1_s["lmass"] >= 11)]
	massed_2_s = chunk_2_s[(chunk2_s["lmass"] >= 11)]
	massed_3_s = chunk_3_s[(chunk3_s["lmass"] >= 11)]
	massed_4_s = chunk_4_s[(chunk4_s["lmass"] >= 11)]
	massed_5_s = chunk_5_s[(chunk5_s["lmass"] >= 11)]
	unmassed_s_ = data_z_flagged4[(data_fast_flagged4["lmass"] < 11)]
	
	z1_s = massed_1_s["z"]
	z2_s = massed_2_s["z"]
	z3_s = massed_3_s["z"]
	z4_s = massed_4_s["z"]
	z5_s = massed_5_s["z"]
	unz_s = unmassed_s_["z"]
	
	
	
	
	
	data5 = ascii.read("uds_3dhst.v4.1.cat")
	data_fast5 = ascii.read("uds_3dhst.v4.1.fout")
	data_z5 = ascii.read("uds_3dhst.v4.0.sfr")
	
	idx5, = np.where((data5["use_phot"] == 1.0) & (data5["star_flag"] == 0.0))
	data_fast_flag5 = data_fast5[idx5]
	data_flag5 = data5[idx5]
	data_z_flag5 = data_z5[idx5]
	
	idx_5, = np.where(data_fast_flag5["lsfr"] != -99)
	data_fast_flagged5 = data_fast_flag5[idx_5]
	data_flagged5 = data_flag5[idx_5]
	data_z_flagged5 = data_z_flag5[idx_5]
	
	chunk1_u = data_fast_flagged5[(data_z_flagged5["z"] < 1)]
	chunk2_u = data_fast_flagged5[(data_z_flagged5["z"] >= 1) & (data_z_flagged5["z"] < 1.5)]
	chunk3_u = data_fast_flagged5[(data_z_flagged5["z"] >= 1.5) & (data_z_flagged5["z"] < 2)]
	chunk4_u = data_fast_flagged5[(data_z_flagged5["z"] >= 2) & (data_z_flagged5["z"] <2.5)]
	chunk5_u = data_fast_flagged5[(data_z_flagged5["z"] >= 2.5)]
	
	massed1_u = chunk1_u[(chunk1_u["lmass"] >= 11)]
	massed2_u = chunk2_u[(chunk2_u["lmass"] >= 11)]
	massed3_u = chunk3_u[(chunk3_u["lmass"] >= 11)]
	massed4_u = chunk4_u[(chunk4_u["lmass"] >= 11)]
	massed5_u = chunk5_u[(chunk5_u["lmass"] >= 11)]
	unmassed_u = data_fast_flagged5[(data_fast_flagged5["lmass"] < 11)]
	
	lmass1_u = massed1_u["lmass"]
	lmass2_u = massed2_u["lmass"]
	lmass3_u = massed3_u["lmass"]
	lmass4_u = massed4_u["lmass"]
	lmass5_u = massed5_u["lmass"]
	lunmass_u = unmassed_u["lmass"]
	
	chunk_1_u = data_z_flagged5[(data_z_flagged5["z"] < 1)]
	chunk_2_u = data_z_flagged5[(data_z_flagged5["z"] >= 1) & (data_z_flagged5["z"] < 1.5)]
	chunk_3_u = data_z_flagged5[(data_z_flagged5["z"] >= 1.5) & (data_z_flagged5["z"] < 2)]
	chunk_4_u = data_z_flagged5[(data_z_flagged5["z"] >= 2) & (data_z_flagged5["z"] <2.5)]
	chunk_5_u = data_z_flagged5[(data_z_flagged5["z"] >= 2.5)]
	
	massed_1_u = chunk_1_u[(chunk1_u["lmass"] >= 11)]
	massed_2_u = chunk_2_u[(chunk2_u["lmass"] >= 11)]
	massed_3_u = chunk_3_u[(chunk3_u["lmass"] >= 11)]
	massed_4_u = chunk_4_u[(chunk4_u["lmass"] >= 11)]
	massed_5_u = chunk_5_u[(chunk5_u["lmass"] >= 11)]
	unmassed_u_ = data_z_flagged5[(data_fast_flagged5["lmass"] < 11)]
	
	z1_u = massed_1_u["z"]
	z2_u = massed_2_u["z"]
	z3_u = massed_3_u["z"]
	z4_u = massed_4_u["z"]
	z5_u = massed_5_u["z"]
	unz_u = unmassed_u_["z"]
	
	
	
	
	pylab.scatter(z1_a, lmass1_a, color="0.7", alpha=0.5, edgecolor="none")
	pylab.scatter(z2_a, lmass2_a, color="b", label="chunk1 points, 1 < z < 1.5", alpha=0.7, edgecolor="none")
	pylab.scatter(z3_a, lmass3_a, color="g", label="chunk2 points, 1.5 < z < 2", alpha=0.7, edgecolor="none")
	pylab.scatter(z4_a, lmass4_a, color="purple", label="chunk3 points, 2 < z < 2.5", alpha=0.7, edgecolor="none")
	pylab.scatter(z5_a, lmass5_a, color="0.7", alpha=0.5, edgecolor="none")
	pylab.scatter(unz_a, lunmass_a, color="0.7", alpha=0.5, edgecolor="none", label="unused points")
	
	pylab.scatter(z1_c, lmass1_c, color="0.7", alpha=0.5, edgecolor="none")
	pylab.scatter(z2_c, lmass2_c, color="b", alpha=0.7, edgecolor="none")
	pylab.scatter(z3_c, lmass3_c, color="g", alpha=0.7, edgecolor="none")
	pylab.scatter(z4_c, lmass4_c, color="purple", alpha=0.7, edgecolor="none")
	pylab.scatter(z5_c, lmass5_c, color="0.7", alpha=0.5, edgecolor="none")
	pylab.scatter(unz_c, lunmass_c, color="0.7", alpha=0.5, edgecolor="none")
	
	pylab.scatter(z1_n, lmass1_n, color="0.7", alpha=0.5, edgecolor="none")
	pylab.scatter(z2_n, lmass2_n, color="b", alpha=0.7, edgecolor="none")
	pylab.scatter(z3_n, lmass3_n, color="g", alpha=0.7, edgecolor="none")
	pylab.scatter(z4_n, lmass4_n, color="purple", alpha=0.7, edgecolor="none")
	pylab.scatter(z5_n, lmass5_n, color="0.7", alpha=0.5, edgecolor="none")
	pylab.scatter(unz_n, lunmass_n, color="0.7", alpha=0.5, edgecolor="none")
	
	pylab.scatter(z1_s, lmass1_s, color="0.7", alpha=0.5, edgecolor="none")
	pylab.scatter(z2_s, lmass2_s, color="b", alpha=0.7, edgecolor="none")
	pylab.scatter(z3_s, lmass3_s, color="g", alpha=0.7, edgecolor="none")
	pylab.scatter(z4_s, lmass4_s, color="purple", alpha=0.7, edgecolor="none")
	pylab.scatter(z5_s, lmass5_s, color="0.7", alpha=0.5, edgecolor="none")
	pylab.scatter(unz_s, lunmass_s, color="0.7", alpha=0.5, edgecolor="none")
	
	pylab.scatter(z1_u, lmass1_u, color="0.7", alpha=0.5, edgecolor="none")
	pylab.scatter(z2_u, lmass2_u, color="b", alpha=0.7, edgecolor="none")
	pylab.scatter(z3_u, lmass3_u, color="g", alpha=0.7, edgecolor="none")
	pylab.scatter(z4_u, lmass4_u, color="purple", alpha=0.7, edgecolor="none")
	pylab.scatter(z5_u, lmass5_u, color="0.7", alpha=0.5, edgecolor="none")
	pylab.scatter(unz_u, lunmass_u, color="0.7", alpha=0.5, edgecolor="none")
	
	pylab.xlabel("z")
	pylab.ylabel("log mass")
	pylab.title("Log Mass vs Z")
	pylab.vlines(1, 3, 13, linestyle="dashed")
	pylab.vlines(1.5, 3, 13, linestyle="dashed")
	pylab.vlines(2, 3, 13, linestyle="dashed")
	pylab.vlines(2.5, 3, 13, linestyle="dashed")
	pylab.hlines(11, -0.5, 6, linestyle="dashed")
	pylab.legend(loc=4)
	pylab.xlim([-0.5,6])
	pylab.ylim([3.5,13])
	
	pylab.ion()
	pylab.show()
	
	
def age_mass():
	
	data1 = ascii.read("aegis_3dhst.v4.1.cat")
	data_fast1 = ascii.read("aegis_3dhst.v4.1.fout")
	data_z1 = ascii.read("aegis_3dhst.v4.0.sfr")
	
	idx1, = np.where((data1["use_phot"] == 1.0) & (data1["star_flag"] == 0.0))
	data_fast_flag1 = data_fast1[idx1]
	data_flag1 = data1[idx1]
	data_z_flag1 = data_z1[idx1]
	
	idx_1, = np.where(data_fast_flag1["lsfr"] != -99)
	data_fast_flagged1 = data_fast_flag1[idx_1]
	data_flagged1 = data_flag1[idx_1]
	data_z_flagged1 = data_z_flag1[idx_1]
	
	chunk2_a = data_fast_flagged1[(data_z_flagged1["z"] >= 1) & (data_z_flagged1["z"] < 1.5)]
	chunk3_a = data_fast_flagged1[(data_z_flagged1["z"] >= 1.5) & (data_z_flagged1["z"] < 2)]
	chunk4_a = data_fast_flagged1[(data_z_flagged1["z"] >= 2) & (data_z_flagged1["z"] <2.5)]
	
	massed2_a = chunk2_a[(chunk2_a["lmass"] >= 11)]
	massed3_a = chunk3_a[(chunk3_a["lmass"] >= 11)]
	massed4_a = chunk4_a[(chunk4_a["lmass"] >= 11)]
	
	lmass2_a = massed2_a["lmass"]
	lmass3_a = massed3_a["lmass"]
	lmass4_a = massed4_a["lmass"]
	
	lage2_a = massed2_a["lage"]
	lage3_a = massed3_a["lage"]
	lage4_a = massed4_a["lage"]
	
	
	
	
	data2 = ascii.read("cosmos_3dhst.v4.1.cat")
	data_fast2 = ascii.read("cosmos_3dhst.v4.1.fout")
	data_z2 = ascii.read("cosmos_3dhst.v4.0.sfr")
	
	idx2, = np.where((data2["use_phot"] == 1.0) & (data2["star_flag"] == 0.0))
	data_fast_flag2 = data_fast2[idx2]
	data_flag2 = data2[idx2]
	data_z_flag2 = data_z2[idx2]
	
	idx_2, = np.where(data_fast_flag2["lsfr"] != -99)
	data_fast_flagged2 = data_fast_flag2[idx_2]
	data_flagged2 = data_flag2[idx_2]
	data_z_flagged2 = data_z_flag2[idx_2]
	
	chunk2_c = data_fast_flagged2[(data_z_flagged2["z"] >= 1) & (data_z_flagged2["z"] < 1.5)]
	chunk3_c = data_fast_flagged2[(data_z_flagged2["z"] >= 1.5) & (data_z_flagged2["z"] < 2)]
	chunk4_c = data_fast_flagged2[(data_z_flagged2["z"] >= 2) & (data_z_flagged2["z"] <2.5)]
	
	massed2_c = chunk2_c[(chunk2_c["lmass"] >= 11)]
	massed3_c = chunk3_c[(chunk3_c["lmass"] >= 11)]
	massed4_c = chunk4_c[(chunk4_c["lmass"] >= 11)]
	
	lmass2_c = massed2_c["lmass"]
	lmass3_c = massed3_c["lmass"]
	lmass4_c = massed4_c["lmass"]

	lage2_c = massed2_c["lage"]
	lage3_c = massed3_c["lage"]
	lage4_c = massed4_c["lage"]
	
	
	
	
	data3 = ascii.read("goodsn_3dhst.v4.1.cat")
	data_fast3 = ascii.read("goodsn_3dhst.v4.1.fout")
	data_z3 = ascii.read("goodsn_3dhst.v4.0.sfr")
	
	idx3, = np.where((data3["use_phot"] == 1.0) & (data3["star_flag"] == 0.0))
	data_fast_flag3 = data_fast3[idx3]
	data_flag3 = data3[idx3]
	data_z_flag3 = data_z3[idx3]
	
	idx_3, = np.where(data_fast_flag3["lsfr"] != -99)
	data_fast_flagged3 = data_fast_flag3[idx_3]
	data_flagged3 = data_flag3[idx_3]
	data_z_flagged3 = data_z_flag3[idx_3]
	
	chunk2_n = data_fast_flagged3[(data_z_flagged3["z"] >= 1) & (data_z_flagged3["z"] < 1.5)]
	chunk3_n = data_fast_flagged3[(data_z_flagged3["z"] >= 1.5) & (data_z_flagged3["z"] < 2)]
	chunk4_n = data_fast_flagged3[(data_z_flagged3["z"] >= 2) & (data_z_flagged3["z"] <2.5)]
	
	massed2_n = chunk2_n[(chunk2_n["lmass"] >= 11)]
	massed3_n = chunk3_n[(chunk3_n["lmass"] >= 11)]
	massed4_n = chunk4_n[(chunk4_n["lmass"] >= 11)]
	
	lmass2_n = massed2_n["lmass"]
	lmass3_n = massed3_n["lmass"]
	lmass4_n = massed4_n["lmass"]

	lage2_n = massed2_n["lage"]
	lage3_n = massed3_n["lage"]
	lage4_n = massed4_n["lage"]
	
	
	
	
	data4 = ascii.read("goodss_3dhst.v4.1.cat")
	data_fast4 = ascii.read("goodss_3dhst.v4.1.fout")
	data_z4 = ascii.read("goodss_3dhst.v4.0.sfr")
	
	idx4, = np.where((data4["use_phot"] == 1.0) & (data4["star_flag"] == 0.0))
	data_fast_flag4 = data_fast4[idx4]
	data_flag4 = data4[idx4]
	data_z_flag4 = data_z4[idx4]
	
	idx_4, = np.where(data_fast_flag4["lsfr"] != -99)
	data_fast_flagged4 = data_fast_flag4[idx_4]
	data_flagged4 = data_flag4[idx_4]
	data_z_flagged4 = data_z_flag4[idx_4]
	
	chunk2_s = data_fast_flagged4[(data_z_flagged4["z"] >= 1) & (data_z_flagged4["z"] < 1.5)]
	chunk3_s = data_fast_flagged4[(data_z_flagged4["z"] >= 1.5) & (data_z_flagged4["z"] < 2)]
	chunk4_s = data_fast_flagged4[(data_z_flagged4["z"] >= 2) & (data_z_flagged4["z"] <2.5)]
	
	massed2_s = chunk2_s[(chunk2_s["lmass"] >= 11)]
	massed3_s = chunk3_s[(chunk3_s["lmass"] >= 11)]
	massed4_s = chunk4_s[(chunk4_s["lmass"] >= 11)]
	
	lmass2_s = massed2_s["lmass"]
	lmass3_s = massed3_s["lmass"]
	lmass4_s = massed4_s["lmass"]

	lage2_s = massed2_s["lage"]
	lage3_s = massed3_s["lage"]
	lage4_s = massed4_s["lage"]
	
	
	
	
	data5 = ascii.read("uds_3dhst.v4.1.cat")
	data_fast5 = ascii.read("uds_3dhst.v4.1.fout")
	data_z5 = ascii.read("uds_3dhst.v4.0.sfr")
	
	idx5, = np.where((data5["use_phot"] == 1.0) & (data5["star_flag"] == 0.0))
	data_fast_flag5 = data_fast5[idx5]
	data_flag5 = data5[idx5]
	data_z_flag5 = data_z5[idx5]
	
	idx_5, = np.where(data_fast_flag5["lsfr"] != -99)
	data_fast_flagged5 = data_fast_flag5[idx_5]
	data_flagged5 = data_flag5[idx_5]
	data_z_flagged5 = data_z_flag5[idx_5]
	
	chunk2_u = data_fast_flagged5[(data_z_flagged5["z"] >= 1) & (data_z_flagged5["z"] < 1.5)]
	chunk3_u = data_fast_flagged5[(data_z_flagged5["z"] >= 1.5) & (data_z_flagged5["z"] < 2)]
	chunk4_u = data_fast_flagged5[(data_z_flagged5["z"] >= 2) & (data_z_flagged5["z"] <2.5)]
	
	massed2_u = chunk2_u[(chunk2_u["lmass"] >= 11)]
	massed3_u = chunk3_u[(chunk3_u["lmass"] >= 11)]
	massed4_u = chunk4_u[(chunk4_u["lmass"] >= 11)]
	
	lmass2_u = massed2_u["lmass"]
	lmass3_u = massed3_u["lmass"]
	lmass4_u = massed4_u["lmass"]

	lage2_u = massed2_u["lage"]
	lage3_u = massed3_u["lage"]
	lage4_u = massed4_u["lage"]
	
	
	fig,pylab.axes = pylab.subplots(3, 1)

	a1 = pylab.axes[0]
	a2 = pylab.axes[1]
	a3 = pylab.axes[2]
	
	a1.scatter(lage2_a, lmass2_a, color="b", label="chunk1 points, 1 < z < 1.5", alpha=0.7)
	a2.scatter(lage3_a, lmass3_a, color="g", label="chunk2 points, 1.5 < z < 2", alpha=0.7)
	a3.scatter(lage4_a, lmass4_a, color="purple", label="chunk3 points, 2 < z < 2.5", alpha=0.7)
	
	a1.scatter(lage2_c, lmass2_c, color="b", alpha=0.7)
	a2.scatter(lage3_c, lmass3_c, color="g", alpha=0.7)
	a3.scatter(lage4_c, lmass4_c, color="purple", alpha=0.7)
	
	a1.scatter(lage2_n, lmass2_n, color="b", alpha=0.7)
	a2.scatter(lage3_n, lmass3_n, color="g", alpha=0.7)
	a3.scatter(lage4_n, lmass4_n, color="purple", alpha=0.7)
	
	a1.scatter(lage2_s, lmass2_s, color="b", alpha=0.7)
	a2.scatter(lage3_s, lmass3_s, color="g", alpha=0.7)
	a3.scatter(lage4_s, lmass4_s, color="purple", alpha=0.7)
	
	a1.scatter(lage2_u, lmass2_u, color="b", alpha=0.7)
	a2.scatter(lage3_u, lmass3_u, color="g", alpha=0.7)
	a3.scatter(lage4_u, lmass4_u, color="purple", alpha=0.7)
	

	
	pylab.suptitle("Log Mass vs Age", fontsize=19)
	a1.legend(loc=2)
	a2.legend(loc=2)
	a3.legend(loc=2)
	a1.set_ylim([10.91,11.79])
	a2.set_ylim([10.91,11.79])
	a3.set_ylim([10.91,11.79])

	fig.text(0.44, 0.01, "Age", fontsize=18)
	fig.text(0.01, 0.5, "log mass", rotation = "vertical", fontsize=18)
	fig.subplots_adjust(hspace=0, wspace=0)
	
	pylab.ion()
	pylab.show()

def sersic_mass():
	
	data1 = ascii.read("aegis_3dhst.v4.1.cat")
	data_fast1 = ascii.read("aegis_3dhst.v4.1.fout")
	data_z1 = ascii.read("aegis_3dhst.v4.0.sfr")
	dataz1 = ascii.read("aegis_3dhst.v4.1_f160w.galfit")
	
	idx1, = np.where((data1["use_phot"] == 1.0) & (data1["star_flag"] == 0.0))
	data_fast_flag1 = data_fast1[idx1]
	data_flag1 = data1[idx1]
	data_z_flag1 = data_z1[idx1]
	dataz_flag1 = dataz1[idx1]
	
	idx_1, = np.where(data_fast_flag1["lsfr"] != -99)
	data_fast_flagged1 = data_fast_flag1[idx_1]
	data_flagged1 = data_flag1[idx_1]
	data_z_flagged1 = data_z_flag1[idx_1]
	dataz_flagged1 = dataz_flag1[idx_1]
	
	chunk2_a = data_fast_flagged1[(data_z_flagged1["z"] >= 1) & (data_z_flagged1["z"] < 1.5)]
	chunk3_a = data_fast_flagged1[(data_z_flagged1["z"] >= 1.5) & (data_z_flagged1["z"] < 2)]
	chunk4_a = data_fast_flagged1[(data_z_flagged1["z"] >= 2) & (data_z_flagged1["z"] <2.5)]
	
	massed2_a = chunk2_a[(chunk2_a["lmass"] >= 11)]
	massed3_a = chunk3_a[(chunk3_a["lmass"] >= 11)]
	massed4_a = chunk4_a[(chunk4_a["lmass"] >= 11)]
	
	lmass2_a = massed2_a["lmass"]
	lmass3_a = massed3_a["lmass"]
	lmass4_a = massed4_a["lmass"]
	
	chunk_2_a = dataz_flagged1[(data_z_flagged1["z"] >= 1) & (data_z_flagged1["z"] < 1.5)]
	chunk_3_a = dataz_flagged1[(data_z_flagged1["z"] >= 1.5) & (data_z_flagged1["z"] < 2)]
	chunk_4_a = dataz_flagged1[(data_z_flagged1["z"] >= 2) & (data_z_flagged1["z"] <2.5)]
	
	massed_2_a = chunk_2_a[(chunk2_a["lmass"] >= 11)]
	massed_3_a = chunk_3_a[(chunk3_a["lmass"] >= 11)]
	massed_4_a = chunk_4_a[(chunk4_a["lmass"] >= 11)]
	
	n2_a = massed_2_a["n"]
	n3_a = massed_3_a["n"]
	n4_a = massed_4_a["n"]
	
	
	
	
	data2 = ascii.read("cosmos_3dhst.v4.1.cat")
	data_fast2 = ascii.read("cosmos_3dhst.v4.1.fout")
	data_z2 = ascii.read("cosmos_3dhst.v4.0.sfr")
	dataz2 = ascii.read("cosmos_3dhst.v4.1_f160w.galfit")

	
	idx2, = np.where((data2["use_phot"] == 1.0) & (data2["star_flag"] == 0.0))
	data_fast_flag2 = data_fast2[idx2]
	data_flag2 = data2[idx2]
	data_z_flag2 = data_z2[idx2]
	dataz_flag2 = dataz2[idx2]
	
	idx_2, = np.where(data_fast_flag2["lsfr"] != -99)
	data_fast_flagged2 = data_fast_flag2[idx_2]
	data_flagged2 = data_flag2[idx_2]
	data_z_flagged2 = data_z_flag2[idx_2]
	dataz_flagged2 = dataz_flag2[idx_2]
	
	chunk2_c = data_fast_flagged2[(data_z_flagged2["z"] >= 1) & (data_z_flagged2["z"] < 1.5)]
	chunk3_c = data_fast_flagged2[(data_z_flagged2["z"] >= 1.5) & (data_z_flagged2["z"] < 2)]
	chunk4_c = data_fast_flagged2[(data_z_flagged2["z"] >= 2) & (data_z_flagged2["z"] <2.5)]
	
	massed2_c = chunk2_c[(chunk2_c["lmass"] >= 11)]
	massed3_c = chunk3_c[(chunk3_c["lmass"] >= 11)]
	massed4_c = chunk4_c[(chunk4_c["lmass"] >= 11)]
	
	lmass2_c = massed2_c["lmass"]
	lmass3_c = massed3_c["lmass"]
	lmass4_c = massed4_c["lmass"]

	chunk_2_c = dataz_flagged2[(data_z_flagged2["z"] >= 1) & (data_z_flagged2["z"] < 1.5)]
	chunk_3_c = dataz_flagged2[(data_z_flagged2["z"] >= 1.5) & (data_z_flagged2["z"] < 2)]
	chunk_4_c = dataz_flagged2[(data_z_flagged2["z"] >= 2) & (data_z_flagged2["z"] <2.5)]
	
	massed_2_c = chunk_2_c[(chunk2_c["lmass"] >= 11)]
	massed_3_c = chunk_3_c[(chunk3_c["lmass"] >= 11)]
	massed_4_c = chunk_4_c[(chunk4_c["lmass"] >= 11)]
	
	n2_c = massed_2_c["n"]
	n3_c = massed_3_c["n"]
	n4_c = massed_4_c["n"]
	
	
	
	
	data3 = ascii.read("goodsn_3dhst.v4.1.cat")
	data_fast3 = ascii.read("goodsn_3dhst.v4.1.fout")
	data_z3 = ascii.read("goodsn_3dhst.v4.0.sfr")
	dataz3 = ascii.read("goodsn_3dhst.v4.1_f160w.galfit")

	
	idx3, = np.where((data3["use_phot"] == 1.0) & (data3["star_flag"] == 0.0))
	data_fast_flag3 = data_fast3[idx3]
	data_flag3 = data3[idx3]
	data_z_flag3 = data_z3[idx3]
	dataz_flag3 = dataz3[idx3]
	
	idx_3, = np.where(data_fast_flag3["lsfr"] != -99)
	data_fast_flagged3 = data_fast_flag3[idx_3]
	data_flagged3 = data_flag3[idx_3]
	data_z_flagged3 = data_z_flag3[idx_3]
	dataz_flagged3 = dataz_flag3[idx_3]
	
	chunk2_n = data_fast_flagged3[(data_z_flagged3["z"] >= 1) & (data_z_flagged3["z"] < 1.5)]
	chunk3_n = data_fast_flagged3[(data_z_flagged3["z"] >= 1.5) & (data_z_flagged3["z"] < 2)]
	chunk4_n = data_fast_flagged3[(data_z_flagged3["z"] >= 2) & (data_z_flagged3["z"] <2.5)]
	
	massed2_n = chunk2_n[(chunk2_n["lmass"] >= 11)]
	massed3_n = chunk3_n[(chunk3_n["lmass"] >= 11)]
	massed4_n = chunk4_n[(chunk4_n["lmass"] >= 11)]
	
	lmass2_n = massed2_n["lmass"]
	lmass3_n = massed3_n["lmass"]
	lmass4_n = massed4_n["lmass"]

	chunk_2_n = dataz_flagged3[(data_z_flagged3["z"] >= 1) & (data_z_flagged3["z"] < 1.5)]
	chunk_3_n = dataz_flagged3[(data_z_flagged3["z"] >= 1.5) & (data_z_flagged3["z"] < 2)]
	chunk_4_n = dataz_flagged3[(data_z_flagged3["z"] >= 2) & (data_z_flagged3["z"] <2.5)]
	
	massed_2_n = chunk_2_n[(chunk2_n["lmass"] >= 11)]
	massed_3_n = chunk_3_n[(chunk3_n["lmass"] >= 11)]
	massed_4_n = chunk_4_n[(chunk4_n["lmass"] >= 11)]
	
	n2_n = massed_2_n["n"]
	n3_n = massed_3_n["n"]
	n4_n = massed_4_n["n"]
	
	
	
	
	data4 = ascii.read("goodss_3dhst.v4.1.cat")
	data_fast4 = ascii.read("goodss_3dhst.v4.1.fout")
	data_z4 = ascii.read("goodss_3dhst.v4.0.sfr")
	dataz4 = ascii.read("goodss_3dhst.v4.1_f160w.galfit")

	
	idx4, = np.where((data4["use_phot"] == 1.0) & (data4["star_flag"] == 0.0))
	data_fast_flag4 = data_fast4[idx4]
	data_flag4 = data4[idx4]
	data_z_flag4 = data_z4[idx4]
	dataz_flag4 = dataz4[idx4]
	
	idx_4, = np.where(data_fast_flag4["lsfr"] != -99)
	data_fast_flagged4 = data_fast_flag4[idx_4]
	data_flagged4 = data_flag4[idx_4]
	data_z_flagged4 = data_z_flag4[idx_4]
	dataz_flagged4 = dataz_flag4[idx_4]
	
	chunk2_s = data_fast_flagged4[(data_z_flagged4["z"] >= 1) & (data_z_flagged4["z"] < 1.5)]
	chunk3_s = data_fast_flagged4[(data_z_flagged4["z"] >= 1.5) & (data_z_flagged4["z"] < 2)]
	chunk4_s = data_fast_flagged4[(data_z_flagged4["z"] >= 2) & (data_z_flagged4["z"] <2.5)]
	
	massed2_s = chunk2_s[(chunk2_s["lmass"] >= 11)]
	massed3_s = chunk3_s[(chunk3_s["lmass"] >= 11)]
	massed4_s = chunk4_s[(chunk4_s["lmass"] >= 11)]
	
	lmass2_s = massed2_s["lmass"]
	lmass3_s = massed3_s["lmass"]
	lmass4_s = massed4_s["lmass"]

	chunk_2_s = dataz_flagged4[(data_z_flagged4["z"] >= 1) & (data_z_flagged4["z"] < 1.5)]
	chunk_3_s = dataz_flagged4[(data_z_flagged4["z"] >= 1.5) & (data_z_flagged4["z"] < 2)]
	chunk_4_s = dataz_flagged4[(data_z_flagged4["z"] >= 2) & (data_z_flagged4["z"] <2.5)]
	
	massed_2_s = chunk_2_s[(chunk2_s["lmass"] >= 11)]
	massed_3_s = chunk_3_s[(chunk3_s["lmass"] >= 11)]
	massed_4_s = chunk_4_s[(chunk4_s["lmass"] >= 11)]
	
	n2_s = massed_2_s["n"]
	n3_s = massed_3_s["n"]
	n4_s = massed_4_s["n"]
	
	
	
	
	data5 = ascii.read("uds_3dhst.v4.1.cat")
	data_fast5 = ascii.read("uds_3dhst.v4.1.fout")
	data_z5 = ascii.read("uds_3dhst.v4.0.sfr")
	dataz5 = ascii.read("uds_3dhst.v4.1_f160w.galfit")

	
	idx5, = np.where((data5["use_phot"] == 1.0) & (data5["star_flag"] == 0.0))
	data_fast_flag5 = data_fast5[idx5]
	data_flag5 = data5[idx5]
	data_z_flag5 = data_z5[idx5]
	dataz_flag5 = dataz5[idx5]
	
	idx_5, = np.where(data_fast_flag5["lsfr"] != -99)
	data_fast_flagged5 = data_fast_flag5[idx_5]
	data_flagged5 = data_flag5[idx_5]
	data_z_flagged5 = data_z_flag5[idx_5]
	dataz_flagged5 = dataz_flag5[idx_5]
	
	chunk2_u = data_fast_flagged5[(data_z_flagged5["z"] >= 1) & (data_z_flagged5["z"] < 1.5)]
	chunk3_u = data_fast_flagged5[(data_z_flagged5["z"] >= 1.5) & (data_z_flagged5["z"] < 2)]
	chunk4_u = data_fast_flagged5[(data_z_flagged5["z"] >= 2) & (data_z_flagged5["z"] <2.5)]
	
	massed2_u = chunk2_u[(chunk2_u["lmass"] >= 11)]
	massed3_u = chunk3_u[(chunk3_u["lmass"] >= 11)]
	massed4_u = chunk4_u[(chunk4_u["lmass"] >= 11)]
	
	lmass2_u = massed2_u["lmass"]
	lmass3_u = massed3_u["lmass"]
	lmass4_u = massed4_u["lmass"]

	chunk_2_u = dataz_flagged5[(data_z_flagged5["z"] >= 1) & (data_z_flagged5["z"] < 1.5)]
	chunk_3_u = dataz_flagged5[(data_z_flagged5["z"] >= 1.5) & (data_z_flagged5["z"] < 2)]
	chunk_4_u = dataz_flagged5[(data_z_flagged5["z"] >= 2) & (data_z_flagged5["z"] <2.5)]
	
	massed_2_u = chunk_2_u[(chunk2_u["lmass"] >= 11)]
	massed_3_u = chunk_3_u[(chunk3_u["lmass"] >= 11)]
	massed_4_u = chunk_4_u[(chunk4_u["lmass"] >= 11)]
	
	n2_u = massed_2_u["n"]
	n3_u = massed_3_u["n"]
	n4_u = massed_4_u["n"]
	
	
	fig,pylab.axes = pylab.subplots(3, 1, sharex=True)

	a1 = pylab.axes[0]
	a2 = pylab.axes[1]
	a3 = pylab.axes[2]
	
	a1.scatter(n2_a, lmass2_a, color="b", label="chunk1 points, 1 < z < 1.5", alpha=0.7)
	a2.scatter(n3_a, lmass3_a, color="g", label="chunk2 points, 1.5 < z < 2", alpha=0.7)
	a3.scatter(n4_a, lmass4_a, color="purple", label="chunk3 points, 2 < z < 2.5", alpha=0.7)
	
	a1.scatter(n2_c, lmass2_c, color="b", alpha=0.7)
	a2.scatter(n3_c, lmass3_c, color="g", alpha=0.7)
	a3.scatter(n4_c, lmass4_c, color="purple", alpha=0.7)
	
	a1.scatter(n2_n, lmass2_n, color="b", alpha=0.7)
	a2.scatter(n3_n, lmass3_n, color="g", alpha=0.7)
	a3.scatter(n4_n, lmass4_n, color="purple", alpha=0.7)
	
	a1.scatter(n2_s, lmass2_s, color="b", alpha=0.7)
	a2.scatter(n3_s, lmass3_s, color="g", alpha=0.7)
	a3.scatter(n4_s, lmass4_s, color="purple", alpha=0.7)
	
	a1.scatter(n2_u, lmass2_u, color="b", alpha=0.7)
	a2.scatter(n3_u, lmass3_u, color="g", alpha=0.7)
	a3.scatter(n4_u, lmass4_u, color="purple", alpha=0.7)
	

	
	pylab.suptitle("Log Mass vs Sersic Index", fontsize=19)
	a1.set_xlim([0,8.5])
	a1.set_ylim([10.91,11.79])
	a2.set_xlim([0,8.5])
	a2.set_ylim([10.91,11.79])
	a3.set_xlim([0,8.5])
	a3.set_ylim([10.91,11.79])
	a1.legend()
	a2.legend()
	a3.legend()
	fig.text(0.44, 0.01, "Sersic index", fontsize=18)
	fig.text(0.01, 0.5, "log mass", rotation = "vertical", fontsize=18)
	fig.subplots_adjust(hspace=0, wspace=0)

	
	pylab.ion()
	pylab.show()


def size_mass():
	
	data1 = ascii.read("aegis_3dhst.v4.1.cat")
	data_fast1 = ascii.read("aegis_3dhst.v4.1.fout")
	data_z1 = ascii.read("aegis_3dhst.v4.0.sfr")
	dataz1 = ascii.read("aegis_3dhst.v4.1_f160w.galfit")
	
	idx1, = np.where((data1["use_phot"] == 1.0) & (data1["star_flag"] == 0.0))
	data_fast_flag1 = data_fast1[idx1]
	data_flag1 = data1[idx1]
	data_z_flag1 = data_z1[idx1]
	dataz_flag1 = dataz1[idx1]
	
	idx_1, = np.where(data_fast_flag1["lsfr"] != -99)
	data_fast_flagged1 = data_fast_flag1[idx_1]
	data_flagged1 = data_flag1[idx_1]
	data_z_flagged1 = data_z_flag1[idx_1]
	dataz_flagged1 = dataz_flag1[idx_1]
	
	chunk2_a = data_fast_flagged1[(data_z_flagged1["z"] >= 1) & (data_z_flagged1["z"] < 1.5)]
	chunk3_a = data_fast_flagged1[(data_z_flagged1["z"] >= 1.5) & (data_z_flagged1["z"] < 2)]
	chunk4_a = data_fast_flagged1[(data_z_flagged1["z"] >= 2) & (data_z_flagged1["z"] <2.5)]
	
	massed2_a = chunk2_a[(chunk2_a["lmass"] >= 11)]
	massed3_a = chunk3_a[(chunk3_a["lmass"] >= 11)]
	massed4_a = chunk4_a[(chunk4_a["lmass"] >= 11)]
	
	lmass2_a = massed2_a["lmass"]
	lmass3_a = massed3_a["lmass"]
	lmass4_a = massed4_a["lmass"]
	
	chunk_2_a = dataz_flagged1[(data_z_flagged1["z"] >= 1) & (data_z_flagged1["z"] < 1.5)]
	chunk_3_a = dataz_flagged1[(data_z_flagged1["z"] >= 1.5) & (data_z_flagged1["z"] < 2)]
	chunk_4_a = dataz_flagged1[(data_z_flagged1["z"] >= 2) & (data_z_flagged1["z"] <2.5)]
	
	massed_2_a = chunk_2_a[(chunk2_a["lmass"] >= 11)]
	massed_3_a = chunk_3_a[(chunk3_a["lmass"] >= 11)]
	massed_4_a = chunk_4_a[(chunk4_a["lmass"] >= 11)]
	
	n2_a = massed_2_a["re"]
	n3_a = massed_3_a["re"]
	n4_a = massed_4_a["re"]
	
	
	
	
	data2 = ascii.read("cosmos_3dhst.v4.1.cat")
	data_fast2 = ascii.read("cosmos_3dhst.v4.1.fout")
	data_z2 = ascii.read("cosmos_3dhst.v4.0.sfr")
	dataz2 = ascii.read("cosmos_3dhst.v4.1_f160w.galfit")

	
	idx2, = np.where((data2["use_phot"] == 1.0) & (data2["star_flag"] == 0.0))
	data_fast_flag2 = data_fast2[idx2]
	data_flag2 = data2[idx2]
	data_z_flag2 = data_z2[idx2]
	dataz_flag2 = dataz2[idx2]
	
	idx_2, = np.where(data_fast_flag2["lsfr"] != -99)
	data_fast_flagged2 = data_fast_flag2[idx_2]
	data_flagged2 = data_flag2[idx_2]
	data_z_flagged2 = data_z_flag2[idx_2]
	dataz_flagged2 = dataz_flag2[idx_2]
	
	chunk2_c = data_fast_flagged2[(data_z_flagged2["z"] >= 1) & (data_z_flagged2["z"] < 1.5)]
	chunk3_c = data_fast_flagged2[(data_z_flagged2["z"] >= 1.5) & (data_z_flagged2["z"] < 2)]
	chunk4_c = data_fast_flagged2[(data_z_flagged2["z"] >= 2) & (data_z_flagged2["z"] <2.5)]
	
	massed2_c = chunk2_c[(chunk2_c["lmass"] >= 11)]
	massed3_c = chunk3_c[(chunk3_c["lmass"] >= 11)]
	massed4_c = chunk4_c[(chunk4_c["lmass"] >= 11)]
	
	lmass2_c = massed2_c["lmass"]
	lmass3_c = massed3_c["lmass"]
	lmass4_c = massed4_c["lmass"]

	chunk_2_c = dataz_flagged2[(data_z_flagged2["z"] >= 1) & (data_z_flagged2["z"] < 1.5)]
	chunk_3_c = dataz_flagged2[(data_z_flagged2["z"] >= 1.5) & (data_z_flagged2["z"] < 2)]
	chunk_4_c = dataz_flagged2[(data_z_flagged2["z"] >= 2) & (data_z_flagged2["z"] <2.5)]
	
	massed_2_c = chunk_2_c[(chunk2_c["lmass"] >= 11)]
	massed_3_c = chunk_3_c[(chunk3_c["lmass"] >= 11)]
	massed_4_c = chunk_4_c[(chunk4_c["lmass"] >= 11)]
	
	n2_c = massed_2_c["re"]
	n3_c = massed_3_c["re"]
	n4_c = massed_4_c["re"]
	
	
	
	
	data3 = ascii.read("goodsn_3dhst.v4.1.cat")
	data_fast3 = ascii.read("goodsn_3dhst.v4.1.fout")
	data_z3 = ascii.read("goodsn_3dhst.v4.0.sfr")
	dataz3 = ascii.read("goodsn_3dhst.v4.1_f160w.galfit")

	
	idx3, = np.where((data3["use_phot"] == 1.0) & (data3["star_flag"] == 0.0))
	data_fast_flag3 = data_fast3[idx3]
	data_flag3 = data3[idx3]
	data_z_flag3 = data_z3[idx3]
	dataz_flag3 = dataz3[idx3]
	
	idx_3, = np.where(data_fast_flag3["lsfr"] != -99)
	data_fast_flagged3 = data_fast_flag3[idx_3]
	data_flagged3 = data_flag3[idx_3]
	data_z_flagged3 = data_z_flag3[idx_3]
	dataz_flagged3 = dataz_flag3[idx_3]
	
	chunk2_n = data_fast_flagged3[(data_z_flagged3["z"] >= 1) & (data_z_flagged3["z"] < 1.5)]
	chunk3_n = data_fast_flagged3[(data_z_flagged3["z"] >= 1.5) & (data_z_flagged3["z"] < 2)]
	chunk4_n = data_fast_flagged3[(data_z_flagged3["z"] >= 2) & (data_z_flagged3["z"] <2.5)]
	
	massed2_n = chunk2_n[(chunk2_n["lmass"] >= 11)]
	massed3_n = chunk3_n[(chunk3_n["lmass"] >= 11)]
	massed4_n = chunk4_n[(chunk4_n["lmass"] >= 11)]
	
	lmass2_n = massed2_n["lmass"]
	lmass3_n = massed3_n["lmass"]
	lmass4_n = massed4_n["lmass"]

	chunk_2_n = dataz_flagged3[(data_z_flagged3["z"] >= 1) & (data_z_flagged3["z"] < 1.5)]
	chunk_3_n = dataz_flagged3[(data_z_flagged3["z"] >= 1.5) & (data_z_flagged3["z"] < 2)]
	chunk_4_n = dataz_flagged3[(data_z_flagged3["z"] >= 2) & (data_z_flagged3["z"] <2.5)]
	
	massed_2_n = chunk_2_n[(chunk2_n["lmass"] >= 11)]
	massed_3_n = chunk_3_n[(chunk3_n["lmass"] >= 11)]
	massed_4_n = chunk_4_n[(chunk4_n["lmass"] >= 11)]
	
	n2_n = massed_2_n["re"]
	n3_n = massed_3_n["re"]
	n4_n = massed_4_n["re"]
	
	
	
	
	data4 = ascii.read("goodss_3dhst.v4.1.cat")
	data_fast4 = ascii.read("goodss_3dhst.v4.1.fout")
	data_z4 = ascii.read("goodss_3dhst.v4.0.sfr")
	dataz4 = ascii.read("goodss_3dhst.v4.1_f160w.galfit")

	
	idx4, = np.where((data4["use_phot"] == 1.0) & (data4["star_flag"] == 0.0))
	data_fast_flag4 = data_fast4[idx4]
	data_flag4 = data4[idx4]
	data_z_flag4 = data_z4[idx4]
	dataz_flag4 = dataz4[idx4]
	
	idx_4, = np.where(data_fast_flag4["lsfr"] != -99)
	data_fast_flagged4 = data_fast_flag4[idx_4]
	data_flagged4 = data_flag4[idx_4]
	data_z_flagged4 = data_z_flag4[idx_4]
	dataz_flagged4 = dataz_flag4[idx_4]
	
	chunk2_s = data_fast_flagged4[(data_z_flagged4["z"] >= 1) & (data_z_flagged4["z"] < 1.5)]
	chunk3_s = data_fast_flagged4[(data_z_flagged4["z"] >= 1.5) & (data_z_flagged4["z"] < 2)]
	chunk4_s = data_fast_flagged4[(data_z_flagged4["z"] >= 2) & (data_z_flagged4["z"] <2.5)]
	
	massed2_s = chunk2_s[(chunk2_s["lmass"] >= 11)]
	massed3_s = chunk3_s[(chunk3_s["lmass"] >= 11)]
	massed4_s = chunk4_s[(chunk4_s["lmass"] >= 11)]
	
	lmass2_s = massed2_s["lmass"]
	lmass3_s = massed3_s["lmass"]
	lmass4_s = massed4_s["lmass"]

	chunk_2_s = dataz_flagged4[(data_z_flagged4["z"] >= 1) & (data_z_flagged4["z"] < 1.5)]
	chunk_3_s = dataz_flagged4[(data_z_flagged4["z"] >= 1.5) & (data_z_flagged4["z"] < 2)]
	chunk_4_s = dataz_flagged4[(data_z_flagged4["z"] >= 2) & (data_z_flagged4["z"] <2.5)]
	
	massed_2_s = chunk_2_s[(chunk2_s["lmass"] >= 11)]
	massed_3_s = chunk_3_s[(chunk3_s["lmass"] >= 11)]
	massed_4_s = chunk_4_s[(chunk4_s["lmass"] >= 11)]
	
	n2_s = massed_2_s["re"]
	n3_s = massed_3_s["re"]
	n4_s = massed_4_s["re"]
	
	
	
	
	data5 = ascii.read("uds_3dhst.v4.1.cat")
	data_fast5 = ascii.read("uds_3dhst.v4.1.fout")
	data_z5 = ascii.read("uds_3dhst.v4.0.sfr")
	dataz5 = ascii.read("uds_3dhst.v4.1_f160w.galfit")

	
	idx5, = np.where((data5["use_phot"] == 1.0) & (data5["star_flag"] == 0.0))
	data_fast_flag5 = data_fast5[idx5]
	data_flag5 = data5[idx5]
	data_z_flag5 = data_z5[idx5]
	dataz_flag5 = dataz5[idx5]
	
	idx_5, = np.where(data_fast_flag5["lsfr"] != -99)
	data_fast_flagged5 = data_fast_flag5[idx_5]
	data_flagged5 = data_flag5[idx_5]
	data_z_flagged5 = data_z_flag5[idx_5]
	dataz_flagged5 = dataz_flag5[idx_5]
	
	chunk2_u = data_fast_flagged5[(data_z_flagged5["z"] >= 1) & (data_z_flagged5["z"] < 1.5)]
	chunk3_u = data_fast_flagged5[(data_z_flagged5["z"] >= 1.5) & (data_z_flagged5["z"] < 2)]
	chunk4_u = data_fast_flagged5[(data_z_flagged5["z"] >= 2) & (data_z_flagged5["z"] <2.5)]
	
	massed2_u = chunk2_u[(chunk2_u["lmass"] >= 11)]
	massed3_u = chunk3_u[(chunk3_u["lmass"] >= 11)]
	massed4_u = chunk4_u[(chunk4_u["lmass"] >= 11)]
	
	lmass2_u = massed2_u["lmass"]
	lmass3_u = massed3_u["lmass"]
	lmass4_u = massed4_u["lmass"]

	chunk_2_u = dataz_flagged5[(data_z_flagged5["z"] >= 1) & (data_z_flagged5["z"] < 1.5)]
	chunk_3_u = dataz_flagged5[(data_z_flagged5["z"] >= 1.5) & (data_z_flagged5["z"] < 2)]
	chunk_4_u = dataz_flagged5[(data_z_flagged5["z"] >= 2) & (data_z_flagged5["z"] <2.5)]
	
	massed_2_u = chunk_2_u[(chunk2_u["lmass"] >= 11)]
	massed_3_u = chunk_3_u[(chunk3_u["lmass"] >= 11)]
	massed_4_u = chunk_4_u[(chunk4_u["lmass"] >= 11)]
	
	n2_u = massed_2_u["re"]
	n3_u = massed_3_u["re"]
	n4_u = massed_4_u["re"]
	
	
	
	fig,pylab.axes = pylab.subplots(3, 1, sharex=True)

	a1 = pylab.axes[0]
	a2 = pylab.axes[1]
	a3 = pylab.axes[2]
	
	
	a1.scatter(n2_a, lmass2_a, color="b", label="chunk1 points, 1 < z < 1.5", alpha=0.7)
	a2.scatter(n3_a, lmass3_a, color="g", label="chunk2 points, 1.5 < z < 2", alpha=0.7)
	a3.scatter(n4_a, lmass4_a, color="purple", label="chunk3 points, 2 < z < 2.5", alpha=0.7)
	
	a1.scatter(n2_c, lmass2_c, color="b", alpha=0.7)
	a2.scatter(n3_c, lmass3_c, color="g", alpha=0.7)
	a3.scatter(n4_c, lmass4_c, color="purple", alpha=0.7)
	
	a1.scatter(n2_n, lmass2_n, color="b", alpha=0.7)
	a2.scatter(n3_n, lmass3_n, color="g", alpha=0.7)
	a3.scatter(n4_n, lmass4_n, color="purple", alpha=0.7)
	
	a1.scatter(n2_s, lmass2_s, color="b", alpha=0.7)
	a2.scatter(n3_s, lmass3_s, color="g", alpha=0.7)
	a3.scatter(n4_s, lmass4_s, color="purple", alpha=0.7)
	
	a1.scatter(n2_u, lmass2_u, color="b", alpha=0.7)
	a2.scatter(n3_u, lmass3_u, color="g", alpha=0.7)
	a3.scatter(n4_u, lmass4_u, color="purple", alpha=0.7)
	

	
	pylab.suptitle("Log Mass vs Size", fontsize=19)
	a1.legend()
	a2.legend()
	a3.legend()
	a1.set_xlim([0,2.6])
	a1.set_ylim([10.91,11.79])
	a2.set_xlim([0,2.6])
	a2.set_ylim([10.91,11.79])
	a3.set_xlim([0,2.6])
	a3.set_ylim([10.91,11.79])
	fig.text(0.3, 0.01, "Size in arcseconds of semimajor axis", fontsize=18)
	fig.text(0.01, 0.5, "log mass", rotation = "vertical", fontsize=18)
	fig.subplots_adjust(hspace=0, wspace=0)
	
	pylab.ion()
	pylab.show()





def sfr_mass():
	
	data1 = ascii.read("aegis_3dhst.v4.1.cat")
	data_fast1 = ascii.read("aegis_3dhst.v4.1.fout")
	data_z1 = ascii.read("aegis_3dhst.v4.0.sfr")
	
	idx1, = np.where((data1["use_phot"] == 1.0) & (data1["star_flag"] == 0.0))
	data_fast_flag1 = data_fast1[idx1]
	data_flag1 = data1[idx1]
	data_z_flag1 = data_z1[idx1]
	
	idx_1, = np.where(data_fast_flag1["lsfr"] != -99)
	data_fast_flagged1 = data_fast_flag1[idx_1]
	data_flagged1 = data_flag1[idx_1]
	data_z_flagged1 = data_z_flag1[idx_1]
	
	chunk2_a = data_fast_flagged1[(data_z_flagged1["z"] >= 1) & (data_z_flagged1["z"] < 1.5)]
	chunk3_a = data_fast_flagged1[(data_z_flagged1["z"] >= 1.5) & (data_z_flagged1["z"] < 2)]
	chunk4_a = data_fast_flagged1[(data_z_flagged1["z"] >= 2) & (data_z_flagged1["z"] <2.5)]
	
	massed2_a = chunk2_a[(chunk2_a["lmass"] >= 11)]
	massed3_a = chunk3_a[(chunk3_a["lmass"] >= 11)]
	massed4_a = chunk4_a[(chunk4_a["lmass"] >= 11)]
	
	lmass2_a = massed2_a["lmass"]
	lmass3_a = massed3_a["lmass"]
	lmass4_a = massed4_a["lmass"]
	
	chunk_2_a = data_z_flagged1[(data_z_flagged1["z"] >= 1) & (data_z_flagged1["z"] < 1.5)]
	chunk_3_a = data_z_flagged1[(data_z_flagged1["z"] >= 1.5) & (data_z_flagged1["z"] < 2)]
	chunk_4_a = data_z_flagged1[(data_z_flagged1["z"] >= 2) & (data_z_flagged1["z"] <2.5)]
	
	massed_2_a = chunk_2_a[(chunk2_a["lmass"] >= 11)]
	massed_3_a = chunk_3_a[(chunk3_a["lmass"] >= 11)]
	massed_4_a = chunk_4_a[(chunk4_a["lmass"] >= 11)]
	
	lage2_a = massed_2_a["sfr"]
	lage3_a = massed_3_a["sfr"]
	lage4_a = massed_4_a["sfr"]
	
	
	
	
	data2 = ascii.read("cosmos_3dhst.v4.1.cat")
	data_fast2 = ascii.read("cosmos_3dhst.v4.1.fout")
	data_z2 = ascii.read("cosmos_3dhst.v4.0.sfr")
	
	idx2, = np.where((data2["use_phot"] == 1.0) & (data2["star_flag"] == 0.0))
	data_fast_flag2 = data_fast2[idx2]
	data_flag2 = data2[idx2]
	data_z_flag2 = data_z2[idx2]
	
	idx_2, = np.where(data_fast_flag2["lsfr"] != -99)
	data_fast_flagged2 = data_fast_flag2[idx_2]
	data_flagged2 = data_flag2[idx_2]
	data_z_flagged2 = data_z_flag2[idx_2]
	
	chunk2_c = data_fast_flagged2[(data_z_flagged2["z"] >= 1) & (data_z_flagged2["z"] < 1.5)]
	chunk3_c = data_fast_flagged2[(data_z_flagged2["z"] >= 1.5) & (data_z_flagged2["z"] < 2)]
	chunk4_c = data_fast_flagged2[(data_z_flagged2["z"] >= 2) & (data_z_flagged2["z"] <2.5)]
	
	massed2_c = chunk2_c[(chunk2_c["lmass"] >= 11)]
	massed3_c = chunk3_c[(chunk3_c["lmass"] >= 11)]
	massed4_c = chunk4_c[(chunk4_c["lmass"] >= 11)]
	
	lmass2_c = massed2_c["lmass"]
	lmass3_c = massed3_c["lmass"]
	lmass4_c = massed4_c["lmass"]
	
	chunk_2_c = data_z_flagged2[(data_z_flagged2["z"] >= 1) & (data_z_flagged2["z"] < 1.5)]
	chunk_3_c = data_z_flagged2[(data_z_flagged2["z"] >= 1.5) & (data_z_flagged2["z"] < 2)]
	chunk_4_c = data_z_flagged2[(data_z_flagged2["z"] >= 2) & (data_z_flagged2["z"] <2.5)]
	
	massed_2_c = chunk_2_c[(chunk2_c["lmass"] >= 11)]
	massed_3_c = chunk_3_c[(chunk3_c["lmass"] >= 11)]
	massed_4_c = chunk_4_c[(chunk4_c["lmass"] >= 11)]

	lage2_c = massed_2_c["sfr"]
	lage3_c = massed_3_c["sfr"]
	lage4_c = massed_4_c["sfr"]
	
	
	
	
	data3 = ascii.read("goodsn_3dhst.v4.1.cat")
	data_fast3 = ascii.read("goodsn_3dhst.v4.1.fout")
	data_z3 = ascii.read("goodsn_3dhst.v4.0.sfr")
	
	idx3, = np.where((data3["use_phot"] == 1.0) & (data3["star_flag"] == 0.0))
	data_fast_flag3 = data_fast3[idx3]
	data_flag3 = data3[idx3]
	data_z_flag3 = data_z3[idx3]
	
	idx_3, = np.where(data_fast_flag3["lsfr"] != -99)
	data_fast_flagged3 = data_fast_flag3[idx_3]
	data_flagged3 = data_flag3[idx_3]
	data_z_flagged3 = data_z_flag3[idx_3]
	
	chunk2_n = data_fast_flagged3[(data_z_flagged3["z"] >= 1) & (data_z_flagged3["z"] < 1.5)]
	chunk3_n = data_fast_flagged3[(data_z_flagged3["z"] >= 1.5) & (data_z_flagged3["z"] < 2)]
	chunk4_n = data_fast_flagged3[(data_z_flagged3["z"] >= 2) & (data_z_flagged3["z"] <2.5)]
	
	massed2_n = chunk2_n[(chunk2_n["lmass"] >= 11)]
	massed3_n = chunk3_n[(chunk3_n["lmass"] >= 11)]
	massed4_n = chunk4_n[(chunk4_n["lmass"] >= 11)]
	
	lmass2_n = massed2_n["lmass"]
	lmass3_n = massed3_n["lmass"]
	lmass4_n = massed4_n["lmass"]
	
	chunk_2_n = data_z_flagged3[(data_z_flagged3["z"] >= 1) & (data_z_flagged3["z"] < 1.5)]
	chunk_3_n = data_z_flagged3[(data_z_flagged3["z"] >= 1.5) & (data_z_flagged3["z"] < 2)]
	chunk_4_n = data_z_flagged3[(data_z_flagged3["z"] >= 2) & (data_z_flagged3["z"] <2.5)]
	
	massed_2_n = chunk_2_n[(chunk2_n["lmass"] >= 11)]
	massed_3_n = chunk_3_n[(chunk3_n["lmass"] >= 11)]
	massed_4_n = chunk_4_n[(chunk4_n["lmass"] >= 11)]

	lage2_n = massed_2_n["sfr"]
	lage3_n = massed_3_n["sfr"]
	lage4_n = massed_4_n["sfr"]
	
	
	
	
	data4 = ascii.read("goodss_3dhst.v4.1.cat")
	data_fast4 = ascii.read("goodss_3dhst.v4.1.fout")
	data_z4 = ascii.read("goodss_3dhst.v4.0.sfr")
	
	idx4, = np.where((data4["use_phot"] == 1.0) & (data4["star_flag"] == 0.0))
	data_fast_flag4 = data_fast4[idx4]
	data_flag4 = data4[idx4]
	data_z_flag4 = data_z4[idx4]
	
	idx_4, = np.where(data_fast_flag4["lsfr"] != -99)
	data_fast_flagged4 = data_fast_flag4[idx_4]
	data_flagged4 = data_flag4[idx_4]
	data_z_flagged4 = data_z_flag4[idx_4]
	
	chunk2_s = data_fast_flagged4[(data_z_flagged4["z"] >= 1) & (data_fast_flagged4["z"] < 1.5)]
	chunk3_s = data_fast_flagged4[(data_fast_flagged4["z"] >= 1.5) & (data_z_flagged4["z"] < 2)]
	chunk4_s = data_fast_flagged4[(data_z_flagged4["z"] >= 2) & (data_z_flagged4["z"] <2.5)]
	
	massed2_s = chunk2_s[(chunk2_s["lmass"] >= 11)]
	massed3_s = chunk3_s[(chunk3_s["lmass"] >= 11)]
	massed4_s = chunk4_s[(chunk4_s["lmass"] >= 11)]
	
	lmass2_s = massed2_s["lmass"]
	lmass3_s = massed3_s["lmass"]
	lmass4_s = massed4_s["lmass"]
	
	chunk_2_s = data_z_flagged4[(data_z_flagged4["z"] >= 1) & (data_fast_flagged4["z"] < 1.5)]
	chunk_3_s = data_z_flagged4[(data_fast_flagged4["z"] >= 1.5) & (data_z_flagged4["z"] < 2)]
	chunk_4_s = data_z_flagged4[(data_z_flagged4["z"] >= 2) & (data_z_flagged4["z"] <2.5)]
	
	massed_2_s = chunk_2_s[(chunk2_s["lmass"] >= 11)]
	massed_3_s = chunk_3_s[(chunk3_s["lmass"] >= 11)]
	massed_4_s = chunk_4_s[(chunk4_s["lmass"] >= 11)]

	lage2_s = massed_2_s["sfr"]
	lage3_s = massed_3_s["sfr"]
	lage4_s = massed_4_s["sfr"]
	
	
	
	
	data5 = ascii.read("uds_3dhst.v4.1.cat")
	data_fast5 = ascii.read("uds_3dhst.v4.1.fout")
	data_z5 = ascii.read("uds_3dhst.v4.0.sfr")
	
	idx5, = np.where((data5["use_phot"] == 1.0) & (data5["star_flag"] == 0.0))
	data_fast_flag5 = data_fast5[idx5]
	data_flag5 = data5[idx5]
	data_z_flag5 = data_z5[idx5]
	
	idx_5, = np.where(data_fast_flag5["lsfr"] != -99)
	data_fast_flagged5 = data_fast_flag5[idx_5]
	data_flagged5 = data_flag5[idx_5]
	data_z_flagged5 = data_z_flag5[idx_5]
	
	chunk2_u = data_fast_flagged5[(data_z_flagged5["z"] >= 1) & (data_z_flagged5["z"] < 1.5)]
	chunk3_u = data_fast_flagged5[(data_z_flagged5["z"] >= 1.5) & (data_z_flagged5["z"] < 2)]
	chunk4_u = data_fast_flagged5[(data_z_flagged5["z"] >= 2) & (data_z_flagged5["z"] <2.5)]
	
	massed2_u = chunk2_u[(chunk2_u["lmass"] >= 11)]
	massed3_u = chunk3_u[(chunk3_u["lmass"] >= 11)]
	massed4_u = chunk4_u[(chunk4_u["lmass"] >= 11)]
	
	lmass2_u = massed2_u["lmass"]
	lmass3_u = massed3_u["lmass"]
	lmass4_u = massed4_u["lmass"]
	
	chunk_2_u = data_z_flagged5[(data_z_flagged5["z"] >= 1) & (data_z_flagged5["z"] < 1.5)]
	chunk_3_u = data_z_flagged5[(data_z_flagged5["z"] >= 1.5) & (data_z_flagged5["z"] < 2)]
	chunk_4_u = data_z_flagged5[(data_z_flagged5["z"] >= 2) & (data_z_flagged5["z"] <2.5)]
	
	massed_2_u = chunk_2_u[(chunk2_u["lmass"] >= 11)]
	massed_3_u = chunk_3_u[(chunk3_u["lmass"] >= 11)]
	massed_4_u = chunk_4_u[(chunk4_u["lmass"] >= 11)]

	lage2_u = massed_2_u["sfr"]
	lage3_u = massed_3_u["sfr"]
	lage4_u = massed_4_u["sfr"]
	
	
	
	fig,pylab.axes = pylab.subplots(3, 1, sharex=True)

	a1 = pylab.axes[0]
	a2 = pylab.axes[1]
	a3 = pylab.axes[2]
	
	a1.scatter(lage2_a, lmass2_a, color="b", label="chunk1 points, 1 < z < 1.5", alpha=0.7)
	a2.scatter(lage3_a, lmass3_a, color="g", label="chunk2 points, 1.5 < z < 2", alpha=0.7)
	a3.scatter(lage4_a, lmass4_a, color="purple", label="chunk3 points, 2 < z < 2.5", alpha=0.7)
	
	a1.scatter(lage2_c, lmass2_c, color="b", alpha=0.7)
	a2.scatter(lage3_c, lmass3_c, color="g", alpha=0.7)
	a3.scatter(lage4_c, lmass4_c, color="purple", alpha=0.7)
	
	a1.scatter(lage2_n, lmass2_n, color="b", alpha=0.7)
	a2.scatter(lage3_n, lmass3_n, color="g", alpha=0.7)
	a3.scatter(lage4_n, lmass4_n, color="purple", alpha=0.7)
	
	a1.scatter(lage2_s, lmass2_s, color="b", alpha=0.7)
	a2.scatter(lage3_s, lmass3_s, color="g", alpha=0.7)
	a3.scatter(lage4_s, lmass4_s, color="purple", alpha=0.7)
	
	a1.scatter(lage2_u, lmass2_u, color="b", alpha=0.7)
	a2.scatter(lage3_u, lmass3_u, color="g", alpha=0.7)
	a3.scatter(lage4_u, lmass4_u, color="purple", alpha=0.7)
	

	
	
	pylab.suptitle("Log Mass vs SFR", fontsize=19)
	a1.legend(loc=2)
	a2.legend(loc=2)
	a3.legend(loc=2)
	a1.set_ylim([10.91,11.79])
	a2.set_ylim([10.91,11.79])
	a3.set_ylim([10.91,11.79])
	a1.set_xscale("log")
	a2.set_xscale("log")
	a3.set_xscale("log")
	fig.text(0.4, 0.01, "Star formation rate", fontsize=18)
	fig.text(0.01, 0.5, "log mass", rotation = "vertical", fontsize=18)
	fig.subplots_adjust(hspace=0, wspace=0)
	
	pylab.ion()
	pylab.show()




def size():
	
	
	data1 = ascii.read("aegis_3dhst.v4.1.cat")
	data_fast1 = ascii.read("aegis_3dhst.v4.1.fout")
	data_z1 = ascii.read("aegis_3dhst.v4.0.sfr")
	
	idx1, = np.where((data1["use_phot"] == 1.0) & (data1["star_flag"] == 0.0))
	data_fast_flag1 = data_fast1[idx1]
	data_flag1 = data1[idx1]
	data_z_flag1 = data_z1[idx1]
	
	idx_1, = np.where(data_fast_flag1["lsfr"] != -99)
	data_fast_flagged1 = data_fast_flag1[idx_1]
	data_flagged1 = data_flag1[idx_1]
	data_z_flagged1 = data_z_flag1[idx_1]
	
	chunk2_a = data_fast_flagged1[(data_z_flagged1["z"] >= 1) & (data_z_flagged1["z"] < 1.5)]
	chunk3_a = data_fast_flagged1[(data_z_flagged1["z"] >= 1.5) & (data_z_flagged1["z"] < 2)]
	chunk4_a = data_fast_flagged1[(data_z_flagged1["z"] >= 2) & (data_z_flagged1["z"] <2.5)]
	
	massed2_a = chunk2_a[(chunk2_a["lmass"] >= 11)]
	massed3_a = chunk3_a[(chunk3_a["lmass"] >= 11)]
	massed4_a = chunk4_a[(chunk4_a["lmass"] >= 11)]
	massed2_a_ = massed2_a["id"]
	massed3_a_ = massed3_a["id"]
	massed4_a_ = massed4_a["id"]
	
	mass2_a = []
	mass3_a = []
	mass4_a = []
	for i in massed2_a_:
		mass2_a.append(i)
	for i in massed3_a_:
		mass3_a.append(i)
	for i in massed4_a_:
		mass4_a.append(i)
		
	print mass2_a
	print "space"
	print mass3_a
	print "space"
	print mass4_a
	

	
	
	
	
	data2 = ascii.read("cosmos_3dhst.v4.1.cat")
	data_fast2 = ascii.read("cosmos_3dhst.v4.1.fout")
	data_z2 = ascii.read("cosmos_3dhst.v4.0.sfr")
	
	idx2, = np.where((data2["use_phot"] == 1.0) & (data2["star_flag"] == 0.0))
	data_fast_flag2 = data_fast2[idx2]
	data_flag2 = data2[idx2]
	data_z_flag2 = data_z2[idx2]
	
	idx_2, = np.where(data_fast_flag2["lsfr"] != -99)
	data_fast_flagged2 = data_fast_flag2[idx_2]
	data_flagged2 = data_flag2[idx_2]
	data_z_flagged2 = data_z_flag2[idx_2]
	
	chunk2_c = data_fast_flagged2[(data_z_flagged2["z"] >= 1) & (data_z_flagged2["z"] < 1.5)]
	chunk3_c = data_fast_flagged2[(data_z_flagged2["z"] >= 1.5) & (data_z_flagged2["z"] < 2)]
	chunk4_c = data_fast_flagged2[(data_z_flagged2["z"] >= 2) & (data_z_flagged2["z"] <2.5)]
	
	massed2_c = chunk2_c[(chunk2_c["lmass"] >= 11)]
	massed3_c = chunk3_c[(chunk3_c["lmass"] >= 11)]
	massed4_c = chunk4_c[(chunk4_c["lmass"] >= 11)]
	massed2_c_ = massed2_c["id"]
	massed3_c_ = massed3_c["id"]
	massed4_c_ = massed4_c["id"]
	
	mass2_c = []
	mass3_c = []
	mass4_c = []
	for i in massed2_c_:
		mass2_c.append(i)
	for i in massed3_c_:
		mass3_c.append(i)
	for i in massed4_c_:
		mass4_c.append(i)
		
	print mass2_c
	print "space"
	print mass3_c
	print "space"
	print mass4_c
	

	
	
	
	
	data3 = ascii.read("goodsn_3dhst.v4.1.cat")
	data_fast3 = ascii.read("goodsn_3dhst.v4.1.fout")
	data_z3 = ascii.read("goodsn_3dhst.v4.0.sfr")
	
	idx3, = np.where((data3["use_phot"] == 1.0) & (data3["star_flag"] == 0.0))
	data_fast_flag3 = data_fast3[idx3]
	data_flag3 = data3[idx3]
	data_z_flag3 = data_z3[idx3]
	
	idx_3, = np.where(data_fast_flag3["lsfr"] != -99)
	data_fast_flagged3 = data_fast_flag3[idx_3]
	data_flagged3 = data_flag3[idx_3]
	data_z_flagged3 = data_z_flag3[idx_3]
	
	chunk2_n = data_fast_flagged3[(data_z_flagged3["z"] >= 1) & (data_z_flagged3["z"] < 1.5)]
	chunk3_n = data_fast_flagged3[(data_z_flagged3["z"] >= 1.5) & (data_z_flagged3["z"] < 2)]
	chunk4_n = data_fast_flagged3[(data_z_flagged3["z"] >= 2) & (data_z_flagged3["z"] <2.5)]
	
	massed2_n = chunk2_n[(chunk2_n["lmass"] >= 11)]
	massed3_n = chunk3_n[(chunk3_n["lmass"] >= 11)]
	massed4_n = chunk4_n[(chunk4_n["lmass"] >= 11)]
	massed2_n_ = massed2_n["id"]
	massed3_n_ = massed3_n["id"]
	massed4_n_ = massed4_n["id"]
	
	mass2_n = []
	mass3_n = []
	mass4_n = []
	for i in massed2_n_:
		mass2_n.append(i)
	for i in massed3_n_:
		mass3_n.append(i)
	for i in massed4_n_:
		mass4_n.append(i)
		
	print mass2_n
	print "space"
	print mass3_n
	print "space"
	print mass4_n
	

	
	
	
	
	data4 = ascii.read("goodss_3dhst.v4.1.cat")
	data_fast4 = ascii.read("goodss_3dhst.v4.1.fout")
	data_z4 = ascii.read("goodss_3dhst.v4.0.sfr")
	
	idx4, = np.where((data4["use_phot"] == 1.0) & (data4["star_flag"] == 0.0))
	data_fast_flag4 = data_fast4[idx4]
	data_flag4 = data4[idx4]
	data_z_flag4 = data_z4[idx4]
	
	idx_4, = np.where(data_fast_flag4["lsfr"] != -99)
	data_fast_flagged4 = data_fast_flag4[idx_4]
	data_flagged4 = data_flag4[idx_4]
	data_z_flagged4 = data_z_flag4[idx_4]
	
	chunk2_s = data_fast_flagged4[(data_z_flagged4["z"] >= 1) & (data_z_flagged4["z"] < 1.5)]
	chunk3_s = data_fast_flagged4[(data_z_flagged4["z"] >= 1.5) & (data_z_flagged4["z"] < 2)]
	chunk4_s = data_fast_flagged4[(data_z_flagged4["z"] >= 2) & (data_z_flagged4["z"] <2.5)]
	
	massed2_s = chunk2_s[(chunk2_s["lmass"] >= 11)]
	massed3_s = chunk3_s[(chunk3_s["lmass"] >= 11)]
	massed4_s = chunk4_s[(chunk4_s["lmass"] >= 11)]
	massed2_s_ = massed2_s["id"]
	massed3_s_ = massed3_s["id"]
	massed4_s_ = massed4_s["id"]
	
	mass2_s = []
	mass3_s = []
	mass4_s = []
	for i in massed2_s_:
		mass2_s.append(i)
	for i in massed3_s_:
		mass3_s.append(i)
	for i in massed4_s_:
		mass4_s.append(i)
		
	print mass2_s
	print "space"
	print mass3_s
	print "space"
	print mass4_s
	

	
	
	
	data5 = ascii.read("uds_3dhst.v4.1.cat")
	data_fast5 = ascii.read("uds_3dhst.v4.1.fout")
	data_z5 = ascii.read("uds_3dhst.v4.0.sfr")
	
	idx5, = np.where((data5["use_phot"] == 1.0) & (data5["star_flag"] == 0.0))
	data_fast_flag5 = data_fast5[idx5]
	data_flag5 = data5[idx5]
	data_z_flag5 = data_z5[idx5]
	
	idx_5, = np.where(data_fast_flag5["lsfr"] != -99)
	data_fast_flagged5 = data_fast_flag5[idx_5]
	data_flagged5 = data_flag5[idx_5]
	data_z_flagged5 = data_z_flag5[idx_5]
	
	chunk2_u = data_fast_flagged5[(data_z_flagged5["z"] >= 1) & (data_z_flagged5["z"] < 1.5)]
	chunk3_u = data_fast_flagged5[(data_z_flagged5["z"] >= 1.5) & (data_z_flagged5["z"] < 2)]
	chunk4_u = data_fast_flagged5[(data_z_flagged5["z"] >= 2) & (data_z_flagged5["z"] <2.5)]
	
	massed2_u = chunk2_u[(chunk2_u["lmass"] >= 11)]
	massed3_u = chunk3_u[(chunk3_u["lmass"] >= 11)]
	massed4_u = chunk4_u[(chunk4_u["lmass"] >= 11)]
	massed2_u_ = massed2_u["id"]
	massed3_u_ = massed3_u["id"]
	massed4_u_ = massed4_u["id"]
	
	mass2_u = []
	mass3_u = []
	mass4_u = []
	for i in massed2_u_:
		mass2_u.append(i)
	for i in massed3_u_:
		mass3_u.append(i)
	for i in massed4_u_:
		mass4_u.append(i)
		
	print mass2_u
	print "space"
	print mass3_u
	print "space"
	print mass4_u
	



def U_V():
	
	data1 = ascii.read("aegis_3dhst.v4.1.cat")
	data_fast1 = ascii.read("aegis_3dhst.v4.1.fout")
	data_z1 = ascii.read("aegis_3dhst.v4.1.master.RF")
	dataz1 = ascii.read("aegis_3dhst.v4.0.sfr")
	
	idx1, = np.where((data1["use_phot"] == 1.0) & (data1["star_flag"] == 0.0))
	data_fast_flag1 = data_fast1[idx1]
	data_flag1 = data1[idx1]
	data_z_flag1 = data_z1[idx1]
	dataz_flag1 = dataz1[idx1]
	
	idx_1, = np.where(data_fast_flag1["lsfr"] != -99)
	data_fast_flagged1 = data_fast_flag1[idx_1]
	data_flagged1 = data_flag1[idx_1]
	data_z_flagged1 = data_z_flag1[idx_1]
	dataz_flagged1 = dataz_flag1[idx_1]
	
	chunk2_a = data_z_flagged1[(dataz_flagged1["z"] >= 1) & (dataz_flagged1["z"] < 1.5) & (data_fast_flagged1["lmass"] > 11)]
	chunk3_a = data_z_flagged1[(dataz_flagged1["z"] >= 1.5) & (dataz_flagged1["z"] < 2) & (data_fast_flagged1["lmass"] > 11)]
	chunk4_a = data_z_flagged1[(dataz_flagged1["z"] >= 2) & (dataz_flagged1["z"] < 2.5) & (data_fast_flagged1["lmass"] > 11)]
	
	U2_a = chunk2_a["L153"]
	U3_a = chunk3_a["L153"]
	U4_a = chunk4_a["L153"]
	V2_a = chunk2_a["L155"]
	V3_a = chunk3_a["L155"]
	V4_a = chunk4_a["L155"]
	J2_a = chunk2_a["L161"]
	J3_a = chunk3_a["L161"]
	J4_a = chunk4_a["L161"]
	
	UV2_a = np.array([-2.5*math.log10(x) for x in U2_a]) - np.array([-2.5*math.log10(x) for x in V2_a])
	UV3_a = np.array([-2.5*math.log10(x) for x in U3_a]) - np.array([-2.5*math.log10(x) for x in V3_a])
	UV4_a = np.array([-2.5*math.log10(x) for x in U4_a]) - np.array([-2.5*math.log10(x) for x in V4_a])
	
	VJ2_a = np.array([-2.5*math.log10(x) for x in V2_a]) - np.array([-2.5*math.log10(x) for x in J2_a])
	VJ3_a = np.array([-2.5*math.log10(x) for x in V3_a]) - np.array([-2.5*math.log10(x) for x in J3_a])
	VJ4_a = np.array([-2.5*math.log10(x) for x in V4_a]) - np.array([-2.5*math.log10(x) for x in J4_a])
	
	
	
	
	data2 = ascii.read("cosmos_3dhst.v4.1.cat")
	data_fast2 = ascii.read("cosmos_3dhst.v4.1.fout")
	data_z2 = ascii.read("cosmos_3dhst.v4.1.master.RF")
	dataz2 = ascii.read("cosmos_3dhst.v4.0.sfr")
	
	idx2, = np.where((data2["use_phot"] == 1.0) & (data2["star_flag"] == 0.0))
	data_fast_flag2 = data_fast2[idx2]
	data_flag2 = data2[idx2]
	data_z_flag2 = data_z2[idx2]
	dataz_flag2 = dataz2[idx2]
	
	idx_2, = np.where(data_fast_flag2["lsfr"] != -99)
	data_fast_flagged2 = data_fast_flag2[idx_2]
	data_flagged2 = data_flag2[idx_2]
	data_z_flagged2 = data_z_flag2[idx_2]
	dataz_flagged2 = dataz_flag2[idx_2]
	
	chunk2_c = data_z_flagged2[(dataz_flagged2["z"] >= 1) & (dataz_flagged2["z"] < 1.5) & (data_fast_flagged2["lmass"] > 11)]
	chunk3_c = data_z_flagged2[(dataz_flagged2["z"] >= 1.5) & (dataz_flagged2["z"] < 2) & (data_fast_flagged2["lmass"] > 11)]
	chunk4_c = data_z_flagged2[(dataz_flagged2["z"] >= 2) & (dataz_flagged2["z"] < 2.5) & (data_fast_flagged2["lmass"] > 11)]
	
	U2_c = chunk2_c["L153"]
	U3_c = chunk3_c["L153"]
	U4_c = chunk4_c["L153"]
	V2_c = chunk2_c["L155"]
	V3_c = chunk3_c["L155"]
	V4_c = chunk4_c["L155"]
	J2_c = chunk2_c["L161"]
	J3_c = chunk3_c["L161"]
	J4_c = chunk4_c["L161"]
	
	UV2_c = np.array([-2.5*math.log10(x) for x in U2_c]) - np.array([-2.5*math.log10(x) for x in V2_c])
	UV3_c = np.array([-2.5*math.log10(x) for x in U3_c]) - np.array([-2.5*math.log10(x) for x in V3_c])
	UV4_c = np.array([-2.5*math.log10(x) for x in U4_c]) - np.array([-2.5*math.log10(x) for x in V4_c])
	
	VJ2_c = np.array([-2.5*math.log10(x) for x in V2_c]) - np.array([-2.5*math.log10(x) for x in J2_c])
	VJ3_c = np.array([-2.5*math.log10(x) for x in V3_c]) - np.array([-2.5*math.log10(x) for x in J3_c])
	VJ4_c = np.array([-2.5*math.log10(x) for x in V4_c]) - np.array([-2.5*math.log10(x) for x in J4_c])
	
	
	
	
	data3 = ascii.read("goodsn_3dhst.v4.1.cat")
	data_fast3 = ascii.read("goodsn_3dhst.v4.1.fout")
	data_z3 = ascii.read("goodsn_3dhst.v4.1.master.RF")
	dataz3 = ascii.read("goodsn_3dhst.v4.0.sfr")
	
	idx3, = np.where((data3["use_phot"] == 1.0) & (data3["star_flag"] == 0.0))
	data_fast_flag3 = data_fast3[idx3]
	data_flag3 = data3[idx3]
	data_z_flag3 = data_z3[idx3]
	dataz_flag3 = dataz3[idx3]
	
	idx_3, = np.where(data_fast_flag3["lsfr"] != -99)
	data_fast_flagged3 = data_fast_flag3[idx_3]
	data_flagged3 = data_flag3[idx_3]
	data_z_flagged3 = data_z_flag3[idx_3]
	dataz_flagged3 = dataz_flag3[idx_3]
	
	chunk2_n = data_z_flagged3[(dataz_flagged3["z"] >= 1) & (dataz_flagged3["z"] < 1.5) & (data_fast_flagged3["lmass"] > 11)]
	chunk3_n = data_z_flagged3[(dataz_flagged3["z"] >= 1.5) & (dataz_flagged3["z"] < 2) & (data_fast_flagged3["lmass"] > 11)]
	chunk4_n = data_z_flagged3[(dataz_flagged3["z"] >= 2) & (dataz_flagged3["z"] < 2.5) & (data_fast_flagged3["lmass"] > 11)]
	
	U2_n = chunk2_n["L153"]
	U3_n = chunk3_n["L153"]
	U4_n = chunk4_n["L153"]
	V2_n = chunk2_n["L155"]
	V3_n = chunk3_n["L155"]
	V4_n = chunk4_n["L155"]
	J2_n = chunk2_n["L161"]
	J3_n = chunk3_n["L161"]
	J4_n = chunk4_n["L161"]
	
	UV2_n = np.array([-2.5*math.log10(x) for x in U2_n]) - np.array([-2.5*math.log10(x) for x in V2_n])
	UV3_n = np.array([-2.5*math.log10(x) for x in U3_n]) - np.array([-2.5*math.log10(x) for x in V3_n])
	UV4_n = np.array([-2.5*math.log10(x) for x in U4_n]) - np.array([-2.5*math.log10(x) for x in V4_n])
	
	VJ2_n = np.array([-2.5*math.log10(x) for x in V2_n]) - np.array([-2.5*math.log10(x) for x in J2_n])
	VJ3_n = np.array([-2.5*math.log10(x) for x in V3_n]) - np.array([-2.5*math.log10(x) for x in J3_n])
	VJ4_n = np.array([-2.5*math.log10(x) for x in V4_n]) - np.array([-2.5*math.log10(x) for x in J4_n])
	
	
	
	
	data4 = ascii.read("goodss_3dhst.v4.1.cat")
	data_fast4 = ascii.read("goodss_3dhst.v4.1.fout")
	data_z4 = ascii.read("goodss_3dhst.v4.1.master.RF")
	dataz4 = ascii.read("goodss_3dhst.v4.0.sfr")
	
	idx4, = np.where((data4["use_phot"] == 1.0) & (data4["star_flag"] == 0.0))
	data_fast_flag4 = data_fast4[idx4]
	data_flag4 = data4[idx4]
	data_z_flag4 = data_z4[idx4]
	dataz_flag4 = dataz4[idx4]
	
	idx_4, = np.where(data_fast_flag4["lsfr"] != -99)
	data_fast_flagged4 = data_fast_flag4[idx_4]
	data_flagged4 = data_flag4[idx_4]
	data_z_flagged4 = data_z_flag4[idx_4]
	dataz_flagged4 = dataz_flag4[idx_4]
	
	chunk2_s = data_z_flagged4[(dataz_flagged4["z"] >= 1) & (dataz_flagged4["z"] < 1.5) & (data_fast_flagged4["lmass"] > 11)]
	chunk3_s = data_z_flagged4[(dataz_flagged4["z"] >= 1.5) & (dataz_flagged4["z"] < 2) & (data_fast_flagged4["lmass"] > 11)]
	chunk4_s = data_z_flagged4[(dataz_flagged4["z"] >= 2) & (dataz_flagged4["z"] < 2.5) & (data_fast_flagged4["lmass"] > 11)]
	
	U2_s = chunk2_s["L153"]
	U3_s = chunk3_s["L153"]
	U4_s = chunk4_s["L153"]
	V2_s = chunk2_s["L155"]
	V3_s = chunk3_s["L155"]
	V4_s = chunk4_s["L155"]
	J2_s = chunk2_s["L161"]
	J3_s = chunk3_s["L161"]
	J4_s = chunk4_s["L161"]
	
	UV2_s = np.array([-2.5*math.log10(x) for x in U2_s]) - np.array([-2.5*math.log10(x) for x in V2_s])
	UV3_s = np.array([-2.5*math.log10(x) for x in U3_s]) - np.array([-2.5*math.log10(x) for x in V3_s])
	UV4_s = np.array([-2.5*math.log10(x) for x in U4_s]) - np.array([-2.5*math.log10(x) for x in V4_s])
	
	VJ2_s = np.array([-2.5*math.log10(x) for x in V2_s]) - np.array([-2.5*math.log10(x) for x in J2_s])
	VJ3_s = np.array([-2.5*math.log10(x) for x in V3_s]) - np.array([-2.5*math.log10(x) for x in J3_s])
	VJ4_s = np.array([-2.5*math.log10(x) for x in V4_s]) - np.array([-2.5*math.log10(x) for x in J4_s])
	
	
	
	
	data5 = ascii.read("uds_3dhst.v4.1.cat")
	data_fast5 = ascii.read("uds_3dhst.v4.1.fout")
	data_z5 = ascii.read("uds_3dhst.v4.1.master.RF")
	dataz5 = ascii.read("uds_3dhst.v4.0.sfr")
	
	idx5, = np.where((data5["use_phot"] == 1.0) & (data5["star_flag"] == 0.0))
	data_fast_flag5 = data_fast5[idx5]
	data_flag5 = data5[idx5]
	data_z_flag5 = data_z5[idx5]
	dataz_flag5 = dataz5[idx5]
	
	idx_5, = np.where(data_fast_flag5["lsfr"] != -99)
	data_fast_flagged5 = data_fast_flag5[idx_5]
	data_flagged5 = data_flag5[idx_5]
	data_z_flagged5 = data_z_flag5[idx_5]
	dataz_flagged5 = dataz_flag5[idx_5]
	
	chunk2_u = data_z_flagged5[(dataz_flagged5["z"] >= 1) & (dataz_flagged5["z"] < 1.5) & (data_fast_flagged5["lmass"] > 11)]
	chunk3_u = data_z_flagged5[(dataz_flagged5["z"] >= 1.5) & (dataz_flagged5["z"] < 2) & (data_fast_flagged5["lmass"] > 11)]
	chunk4_u = data_z_flagged5[(dataz_flagged5["z"] >= 2) & (dataz_flagged5["z"] < 2.5) & (data_fast_flagged5["lmass"] > 11)]
	
	U2_u = chunk2_u["L153"]
	U3_u = chunk3_u["L153"]
	U4_u = chunk4_u["L153"]
	V2_u = chunk2_u["L155"]
	V3_u = chunk3_u["L155"]
	V4_u = chunk4_u["L155"]
	J2_u = chunk2_u["L161"]
	J3_u = chunk3_u["L161"]
	J4_u = chunk4_u["L161"]
	
	UV2_u = np.array([-2.5*math.log10(x) for x in U2_u]) - np.array([-2.5*math.log10(x) for x in V2_u])
	UV3_u = np.array([-2.5*math.log10(x) for x in U3_u]) - np.array([-2.5*math.log10(x) for x in V3_u])
	UV4_u = np.array([-2.5*math.log10(x) for x in U4_u]) - np.array([-2.5*math.log10(x) for x in V4_u])
	
	VJ2_u = np.array([-2.5*math.log10(x) for x in V2_u]) - np.array([-2.5*math.log10(x) for x in J2_u])
	VJ3_u = np.array([-2.5*math.log10(x) for x in V3_u]) - np.array([-2.5*math.log10(x) for x in J3_u])
	VJ4_u = np.array([-2.5*math.log10(x) for x in V4_u]) - np.array([-2.5*math.log10(x) for x in J4_u])
	
	
	fig,pylab.axes = pylab.subplots(3, 1, sharex=True)

	a1 = pylab.axes[0]
	a2 = pylab.axes[1]
	a3 = pylab.axes[2]
	
	a1.scatter(VJ2_a, UV2_a, color="g", edgecolor="none", label="1 < z < 1.5")
	a2.scatter(VJ3_a, UV3_a, color="b", edgecolor="none", label="1.5 < z < 2")
	a3.scatter(VJ4_a, UV4_a, color="purple", edgecolor="none", label="2 < z < 2.5")
	a1.scatter(VJ2_c, UV2_c, color="g", edgecolor="none")
	a2.scatter(VJ3_c, UV3_c, color="b", edgecolor="none")
	a3.scatter(VJ4_c, UV4_c, color="purple", edgecolor="none")
	a1.scatter(VJ2_n, UV2_n, color="g", edgecolor="none")
	a2.scatter(VJ3_n, UV3_n, color="b", edgecolor="none")
	a3.scatter(VJ4_n, UV4_n, color="purple", edgecolor="none")
	a1.scatter(VJ2_s, UV2_s, color="g", edgecolor="none")
	a2.scatter(VJ3_s, UV3_s, color="b", edgecolor="none")
	a3.scatter(VJ4_s, UV4_s, color="purple", edgecolor="none")
	a1.scatter(VJ2_u, UV2_u, color="g", edgecolor="none")
	a2.scatter(VJ3_u, UV3_u, color="b", edgecolor="none")
	a3.scatter(VJ4_u, UV4_u, color="purple", edgecolor="none")
	pylab.suptitle("V-J vs U-V", fontsize=19)
	a1.set_ylim([0.9,2.3])
	a2.set_ylim([0.9,2.3])
	a3.set_ylim([0.9,2.3])
	a1.set_xlim([0,2.5])
	a2.set_xlim([0,2.5])
	a3.set_xlim([0,2.5])
	a1.legend(loc=2)
	a2.legend(loc=2)
	a3.legend(loc=2)
	a1.vlines(1.6, 1.9, 2.3, color="r", lw=2)
	a2.vlines(1.6, 1.9, 2.3, color="r", lw=2)
	a3.vlines(1.6, 1.9, 2.3, color="r", lw=2)
	a1.hlines(1.3, 0, 0.9, color="r", lw=2)
	a2.hlines(1.3, 0, 0.9, color="r", lw=2)
	a3.hlines(1.3, 0, 0.9, color="r", lw=2)
	a1.plot([1.6, 0.9], [1.9, 1.3], color="r", lw=2)
	a2.plot([1.6, 0.9], [1.9, 1.3], color="r", lw=2)
	a3.plot([1.6, 0.9], [1.9, 1.3], color="r", lw=2)
	fig.text(0.5, 0.01, "V-J", fontsize=18)
	fig.text(0.01, 0.5, "U-V", rotation = "vertical", fontsize=18)
	fig.subplots_adjust(hspace=0, wspace=0)
	
	
	
	pylab.ion()
	pylab.show()



def wavelength():

	factor = 3.0*(10.0**5.56)
	
	os.chdir("/Volumes/TOSHIBA EXT/3d_hst/noah_massive/aegis_massive")
	
	data1_a1 = ascii.read("2250..temp_sed")
	data1_a2 = ascii.read("2481..temp_sed")
	data1_a3 = ascii.read("6115..temp_sed")
	data1_a4 = ascii.read("6234..temp_sed")
	data1_a5 = ascii.read("6424..temp_sed")
	data1_a6 = ascii.read("6691..temp_sed")
	data1_a7 = ascii.read("7771..temp_sed")
	data1_a8 = ascii.read("8838..temp_sed")
	data1_a9 = ascii.read("11022..temp_sed")
	data1_a10 = ascii.read("13524..temp_sed")
	data1_a11 = ascii.read("14012..temp_sed")
	data1_a12 = ascii.read("14056..temp_sed")
	data1_a13 = ascii.read("14499..temp_sed")
	data1_a14 = ascii.read("14609..temp_sed")
	data1_a15 = ascii.read("16238..temp_sed")
	data1_a16 = ascii.read("17297..temp_sed")
	data1_a17 = ascii.read("17570..temp_sed")
	data1_a18 = ascii.read("17775..temp_sed")
	data1_a19 = ascii.read("18045..temp_sed")
	data1_a20 = ascii.read("18257..temp_sed")
	data1_a21 = ascii.read("19970..temp_sed")
	data1_a22 = ascii.read("20106..temp_sed")
	data1_a23 = ascii.read("20250..temp_sed")
	data1_a24 = ascii.read("20794..temp_sed")
	data1_a25 = ascii.read("21028..temp_sed")
	data1_a26 = ascii.read("21357..temp_sed")
	data1_a27 = ascii.read("22126..temp_sed")
	data1_a28 = ascii.read("22678..temp_sed")
	data1_a29 = ascii.read("23011..temp_sed")
	data1_a30 = ascii.read("23089..temp_sed")
	data1_a31 = ascii.read("24456..temp_sed")
	data1_a32 = ascii.read("26850..temp_sed")
	data1_a33 = ascii.read("26884..temp_sed")
	data1_a34 = ascii.read("27177..temp_sed")
	data1_a35 = ascii.read("28328..temp_sed")
	data1_a36 = ascii.read("29863..temp_sed")
	data1_a37 = ascii.read("30393..temp_sed")
	data1_a38 = ascii.read("30421..temp_sed")
	data1_a39 = ascii.read("30735..temp_sed")
	data1_a40 = ascii.read("30920..temp_sed")
	data1_a41 = ascii.read("31326..temp_sed")
	data1_a42 = ascii.read("32114..temp_sed")
	data1_a43 = ascii.read("32425..temp_sed")
	data1_a44 = ascii.read("33028..temp_sed")
	data1_a45 = ascii.read("33158..temp_sed")
	data1_a46 = ascii.read("34141..temp_sed")
	data1_a47 = ascii.read("34254..temp_sed")
	data1_a48 = ascii.read("34722..temp_sed")
	data1_a49 = ascii.read("35604..temp_sed")
	data1_a50 = ascii.read("37919..temp_sed")
	data1_a51 = ascii.read("38065..temp_sed")
	data1_a52 = ascii.read("38130..temp_sed")
	data1_a53 = ascii.read("38167..temp_sed")

	
	lam1_a1 = data1_a1["lambda"]
	lam1_a2 = data1_a2["lambda"]
	lam1_a3 = data1_a3["lambda"]
	lam1_a4 = data1_a4["lambda"]
	lam1_a5 = data1_a5["lambda"]
	lam1_a6 = data1_a6["lambda"]
	lam1_a7 = data1_a7["lambda"]
	lam1_a8 = data1_a8["lambda"]
	lam1_a9 = data1_a9["lambda"]
	lam1_a10 = data1_a10["lambda"]
	lam1_a11 = data1_a11["lambda"]
	lam1_a12 = data1_a12["lambda"]
	lam1_a13 = data1_a13["lambda"]
	lam1_a14 = data1_a14["lambda"]
	lam1_a15 = data1_a15["lambda"]
	lam1_a16 = data1_a16["lambda"]
	lam1_a17 = data1_a17["lambda"]
	lam1_a18 = data1_a18["lambda"]
	lam1_a19 = data1_a19["lambda"]
	lam1_a20 = data1_a20["lambda"]
	lam1_a21 = data1_a21["lambda"]
	lam1_a22 = data1_a22["lambda"]
	lam1_a23 = data1_a23["lambda"]
	lam1_a24 = data1_a24["lambda"]
	lam1_a25 = data1_a25["lambda"]
	lam1_a26 = data1_a26["lambda"]
	lam1_a27 = data1_a27["lambda"]
	lam1_a28 = data1_a28["lambda"]
	lam1_a29 = data1_a29["lambda"]
	lam1_a30 = data1_a30["lambda"]
	lam1_a31 = data1_a31["lambda"]
	lam1_a32 = data1_a32["lambda"]
	lam1_a33 = data1_a33["lambda"]
	lam1_a34 = data1_a34["lambda"]
	lam1_a35 = data1_a35["lambda"]
	lam1_a36 = data1_a36["lambda"]
	lam1_a37 = data1_a37["lambda"]
	lam1_a38 = data1_a38["lambda"]
	lam1_a39 = data1_a39["lambda"]
	lam1_a40 = data1_a40["lambda"]
	lam1_a41 = data1_a41["lambda"]
	lam1_a42 = data1_a42["lambda"]
	lam1_a43 = data1_a43["lambda"]
	lam1_a44 = data1_a44["lambda"]
	lam1_a45 = data1_a45["lambda"]
	lam1_a46 = data1_a46["lambda"]
	lam1_a47 = data1_a47["lambda"]
	lam1_a48 = data1_a48["lambda"]
	lam1_a49 = data1_a49["lambda"]
	lam1_a50 = data1_a50["lambda"]
	lam1_a51 = data1_a51["lambda"]
	lam1_a52 = data1_a52["lambda"]
	lam1_a53 = data1_a53["lambda"]

	
	flux_1_a1 = data1_a1["tempflux"]
	flux_1_a2 = data1_a2["tempflux"]
	flux_1_a3 = data1_a3["tempflux"]
	flux_1_a4 = data1_a4["tempflux"]
	flux_1_a5 = data1_a5["tempflux"]
	flux_1_a6 = data1_a6["tempflux"]
	flux_1_a7 = data1_a7["tempflux"]
	flux_1_a8 = data1_a8["tempflux"]
	flux_1_a9 = data1_a9["tempflux"]
	flux_1_a10 = data1_a10["tempflux"]
	flux_1_a11 = data1_a11["tempflux"]
	flux_1_a12 = data1_a12["tempflux"]
	flux_1_a13 = data1_a13["tempflux"]
	flux_1_a14 = data1_a14["tempflux"]
	flux_1_a15 = data1_a15["tempflux"]
	flux_1_a16 = data1_a16["tempflux"]
	flux_1_a17 = data1_a17["tempflux"]
	flux_1_a18 = data1_a18["tempflux"]
	flux_1_a19 = data1_a19["tempflux"]
	flux_1_a20 = data1_a20["tempflux"]
	flux_1_a21 = data1_a21["tempflux"]
	flux_1_a22 = data1_a22["tempflux"]
	flux_1_a23 = data1_a23["tempflux"]
	flux_1_a24 = data1_a24["tempflux"]
	flux_1_a25 = data1_a25["tempflux"]
	flux_1_a26 = data1_a26["tempflux"]
	flux_1_a27 = data1_a27["tempflux"]
	flux_1_a28 = data1_a28["tempflux"]
	flux_1_a29 = data1_a29["tempflux"]
	flux_1_a30 = data1_a30["tempflux"]
	flux_1_a31 = data1_a31["tempflux"]
	flux_1_a32 = data1_a32["tempflux"]
	flux_1_a33 = data1_a33["tempflux"]
	flux_1_a34 = data1_a34["tempflux"]
	flux_1_a35 = data1_a35["tempflux"]
	flux_1_a36 = data1_a36["tempflux"]
	flux_1_a37 = data1_a37["tempflux"]
	flux_1_a38 = data1_a38["tempflux"]
	flux_1_a39 = data1_a39["tempflux"]
	flux_1_a40 = data1_a40["tempflux"]
	flux_1_a41 = data1_a41["tempflux"]
	flux_1_a42 = data1_a42["tempflux"]
	flux_1_a43 = data1_a43["tempflux"]
	flux_1_a44 = data1_a44["tempflux"]
	flux_1_a45 = data1_a45["tempflux"]
	flux_1_a46 = data1_a46["tempflux"]
	flux_1_a47 = data1_a47["tempflux"]
	flux_1_a48 = data1_a48["tempflux"]
	flux_1_a49 = data1_a49["tempflux"]
	flux_1_a50 = data1_a50["tempflux"]
	flux_1_a51 = data1_a51["tempflux"]
	flux_1_a52 = data1_a52["tempflux"]
	flux_1_a53 = data1_a53["tempflux"]

	
	
	flux1_a1 = (flux_1_a1*(lam1_a1**-2.0))*factor
	flux1_a2 = (flux_1_a2*(lam1_a2**-2.0))*factor
	flux1_a3 = (flux_1_a3*(lam1_a3**-2.0))*factor
	flux1_a4 = (flux_1_a4*(lam1_a4**-2.0))*factor
	flux1_a5 = (flux_1_a5*(lam1_a5**-2.0))*factor
	flux1_a6 = (flux_1_a6*(lam1_a6**-2.0))*factor
	flux1_a7 = (flux_1_a7*(lam1_a7**-2.0))*factor
	flux1_a8 = (flux_1_a8*(lam1_a8**-2.0))*factor
	flux1_a9 = (flux_1_a9*(lam1_a9**-2.0))*factor
	flux1_a10 = (flux_1_a10*(lam1_a10**-2.0))*factor
	flux1_a11 = (flux_1_a11*(lam1_a11**-2.0))*factor
	flux1_a12 = (flux_1_a12*(lam1_a12**-2.0))*factor
	flux1_a13 = (flux_1_a13*(lam1_a13**-2.0))*factor
	flux1_a14 = (flux_1_a14*(lam1_a14**-2.0))*factor
	flux1_a15 = (flux_1_a15*(lam1_a15**-2.0))*factor
	flux1_a16 = (flux_1_a16*(lam1_a16**-2.0))*factor
	flux1_a17 = (flux_1_a17*(lam1_a17**-2.0))*factor
	flux1_a18 = (flux_1_a18*(lam1_a18**-2.0))*factor
	flux1_a19 = (flux_1_a19*(lam1_a19**-2.0))*factor
	flux1_a20 = (flux_1_a20*(lam1_a20**-2.0))*factor
	flux1_a21 = (flux_1_a21*(lam1_a21**-2.0))*factor
	flux1_a22 = (flux_1_a22*(lam1_a22**-2.0))*factor
	flux1_a23 = (flux_1_a23*(lam1_a23**-2.0))*factor
	flux1_a24 = (flux_1_a24*(lam1_a24**-2.0))*factor
	flux1_a25 = (flux_1_a25*(lam1_a25**-2.0))*factor
	flux1_a26 = (flux_1_a26*(lam1_a26**-2.0))*factor
	flux1_a27 = (flux_1_a27*(lam1_a27**-2.0))*factor
	flux1_a28 = (flux_1_a28*(lam1_a28**-2.0))*factor
	flux1_a29 = (flux_1_a29*(lam1_a29**-2.0))*factor
	flux1_a30 = (flux_1_a30*(lam1_a30**-2.0))*factor
	flux1_a31 = (flux_1_a31*(lam1_a31**-2.0))*factor
	flux1_a32 = (flux_1_a32*(lam1_a32**-2.0))*factor
	flux1_a33 = (flux_1_a33*(lam1_a33**-2.0))*factor
	flux1_a34 = (flux_1_a34*(lam1_a34**-2.0))*factor
	flux1_a35 = (flux_1_a35*(lam1_a35**-2.0))*factor
	flux1_a36 = (flux_1_a36*(lam1_a36**-2.0))*factor
	flux1_a37 = (flux_1_a37*(lam1_a37**-2.0))*factor
	flux1_a38 = (flux_1_a38*(lam1_a38**-2.0))*factor
	flux1_a39 = (flux_1_a39*(lam1_a39**-2.0))*factor
	flux1_a40 = (flux_1_a40*(lam1_a40**-2.0))*factor
	flux1_a41 = (flux_1_a41*(lam1_a41**-2.0))*factor
	flux1_a42 = (flux_1_a42*(lam1_a42**-2.0))*factor
	flux1_a43 = (flux_1_a43*(lam1_a43**-2.0))*factor
	flux1_a44 = (flux_1_a44*(lam1_a44**-2.0))*factor
	flux1_a45 = (flux_1_a45*(lam1_a45**-2.0))*factor
	flux1_a46 = (flux_1_a46*(lam1_a46**-2.0))*factor
	flux1_a47 = (flux_1_a47*(lam1_a47**-2.0))*factor
	flux1_a48 = (flux_1_a48*(lam1_a48**-2.0))*factor
	flux1_a49 = (flux_1_a49*(lam1_a49**-2.0))*factor
	flux1_a50 = (flux_1_a50*(lam1_a50**-2.0))*factor
	flux1_a51 = (flux_1_a51*(lam1_a51**-2.0))*factor
	flux1_a52 = (flux_1_a52*(lam1_a52**-2.0))*factor
	flux1_a53 = (flux_1_a53*(lam1_a53**-2.0))*factor

	
	data2_a1 = ascii.read("195..temp_sed")
	data2_a2 = ascii.read("766..temp_sed")
	data2_a3 = ascii.read("1420..temp_sed")
	data2_a4 = ascii.read("1821..temp_sed")
	data2_a5 = ascii.read("1925..temp_sed")
	data2_a6 = ascii.read("2016..temp_sed")
	data2_a7 = ascii.read("2289..temp_sed")
	data2_a8 = ascii.read("2427..temp_sed")
	data2_a9 = ascii.read("3106..temp_sed")
	data2_a10 = ascii.read("3311..temp_sed")
	data2_a11 = ascii.read("4826..temp_sed")
	data2_a12 = ascii.read("5016..temp_sed")
	data2_a13 = ascii.read("6262..temp_sed")
	data2_a14 = ascii.read("6310..temp_sed")
	data2_a15 = ascii.read("6880..temp_sed")
	data2_a16 = ascii.read("7601..temp_sed")
	data2_a17 = ascii.read("8646..temp_sed")
	data2_a18 = ascii.read("10858..temp_sed")
	data2_a19 = ascii.read("11465..temp_sed")
	data2_a20 = ascii.read("11730..temp_sed")
	data2_a21 = ascii.read("12049..temp_sed")
	data2_a22 = ascii.read("12219..temp_sed")
	data2_a23 = ascii.read("14481..temp_sed")
	data2_a24 = ascii.read("14495..temp_sed")
	data2_a25 = ascii.read("15069..temp_sed")
	data2_a26 = ascii.read("15088..temp_sed")
	data2_a27 = ascii.read("16491..temp_sed")
	data2_a28 = ascii.read("17691..temp_sed")
	data2_a29 = ascii.read("18280..temp_sed")
	data2_a30 = ascii.read("18922..temp_sed")
	data2_a31 = ascii.read("19743..temp_sed")
	data2_a32 = ascii.read("21156..temp_sed")
	data2_a33 = ascii.read("22423..temp_sed")
	data2_a34 = ascii.read("22713..temp_sed")
	data2_a35 = ascii.read("22741..temp_sed")
	data2_a36 = ascii.read("22887..temp_sed")
	data2_a37 = ascii.read("23027..temp_sed")
	data2_a38 = ascii.read("23045..temp_sed")
	data2_a39 = ascii.read("23234..temp_sed")
	data2_a40 = ascii.read("24059..temp_sed")
	data2_a41 = ascii.read("24333..temp_sed")
	data2_a42 = ascii.read("24448..temp_sed")
	data2_a43 = ascii.read("25346..temp_sed")
	data2_a44 = ascii.read("25969..temp_sed")
	data2_a45 = ascii.read("26079..temp_sed")
	data2_a46 = ascii.read("26649..temp_sed")
	data2_a47 = ascii.read("27315..temp_sed")
	data2_a48 = ascii.read("28310..temp_sed")
	data2_a49 = ascii.read("28843..temp_sed")
	data2_a50 = ascii.read("28864..temp_sed")
	data2_a51 = ascii.read("29144..temp_sed")
	data2_a52 = ascii.read("29329..temp_sed")
	data2_a53 = ascii.read("29399..temp_sed")
	data2_a54 = ascii.read("31169..temp_sed")
	data2_a55 = ascii.read("31680..temp_sed")
	data2_a56 = ascii.read("31839..temp_sed")
	data2_a57 = ascii.read("33196..temp_sed")
	data2_a58 = ascii.read("33551..temp_sed")
	data2_a59 = ascii.read("33770..temp_sed")
	data2_a60 = ascii.read("34185..temp_sed")
	data2_a61 = ascii.read("35002..temp_sed")
	data2_a62 = ascii.read("35677..temp_sed")
	data2_a63 = ascii.read("36347..temp_sed")
	data2_a64 = ascii.read("37556..temp_sed")
	data2_a65 = ascii.read("37592..temp_sed")
	data2_a66 = ascii.read("38174..temp_sed")
	data2_a67 = ascii.read("38290..temp_sed")
	data2_a68 = ascii.read("39028..temp_sed")
	data2_a69 = ascii.read("39239..temp_sed")
	data2_a70 = ascii.read("40470..temp_sed")
	
	lam2_a1 = data2_a1["lambda"]
	lam2_a2 = data2_a2["lambda"]
	lam2_a3 = data2_a3["lambda"]
	lam2_a4 = data2_a4["lambda"]
	lam2_a5 = data2_a5["lambda"]
	lam2_a6 = data2_a6["lambda"]
	lam2_a7 = data2_a7["lambda"]
	lam2_a8 = data2_a8["lambda"]
	lam2_a9 = data2_a9["lambda"]
	lam2_a10 = data2_a10["lambda"]
	lam2_a11 = data2_a11["lambda"]
	lam2_a12 = data2_a12["lambda"]
	lam2_a13 = data2_a13["lambda"]
	lam2_a14 = data2_a14["lambda"]
	lam2_a15 = data2_a15["lambda"]
	lam2_a16 = data2_a16["lambda"]
	lam2_a17 = data2_a17["lambda"]
	lam2_a18 = data2_a18["lambda"]
	lam2_a19 = data2_a19["lambda"]
	lam2_a20 = data2_a20["lambda"]
	lam2_a21 = data2_a21["lambda"]
	lam2_a22 = data2_a22["lambda"]
	lam2_a23 = data2_a23["lambda"]
	lam2_a24 = data2_a24["lambda"]
	lam2_a25 = data2_a25["lambda"]
	lam2_a26 = data2_a26["lambda"]
	lam2_a27 = data2_a27["lambda"]
	lam2_a28 = data2_a28["lambda"]
	lam2_a29 = data2_a29["lambda"]
	lam2_a30 = data2_a30["lambda"]
	lam2_a31 = data2_a31["lambda"]
	lam2_a32 = data2_a32["lambda"]
	lam2_a33 = data2_a33["lambda"]
	lam2_a34 = data2_a34["lambda"]
	lam2_a35 = data2_a35["lambda"]
	lam2_a36 = data2_a36["lambda"]
	lam2_a37 = data2_a37["lambda"]
	lam2_a38 = data2_a38["lambda"]
	lam2_a39 = data2_a39["lambda"]
	lam2_a40 = data2_a40["lambda"]
	lam2_a41 = data2_a41["lambda"]
	lam2_a42 = data2_a42["lambda"]
	lam2_a43 = data2_a43["lambda"]
	lam2_a44 = data2_a44["lambda"]
	lam2_a45 = data2_a45["lambda"]
	lam2_a46 = data2_a46["lambda"]
	lam2_a47 = data2_a47["lambda"]
	lam2_a48 = data2_a48["lambda"]
	lam2_a49 = data2_a49["lambda"]
	lam2_a50 = data2_a50["lambda"]
	lam2_a51 = data2_a51["lambda"]
	lam2_a52 = data2_a52["lambda"]
	lam2_a53 = data2_a53["lambda"]
	lam2_a54 = data2_a54["lambda"]
	lam2_a55 = data2_a55["lambda"]
	lam2_a56 = data2_a56["lambda"]
	lam2_a57 = data2_a57["lambda"]
	lam2_a58 = data2_a58["lambda"]
	lam2_a59 = data2_a59["lambda"]
	lam2_a60 = data2_a60["lambda"]
	lam2_a61 = data2_a51["lambda"]
	lam2_a62 = data2_a52["lambda"]
	lam2_a63 = data2_a53["lambda"]
	lam2_a64 = data2_a54["lambda"]
	lam2_a65 = data2_a55["lambda"]
	lam2_a66 = data2_a56["lambda"]
	lam2_a67 = data2_a57["lambda"]
	lam2_a68 = data2_a58["lambda"]
	lam2_a69 = data2_a59["lambda"]
	lam2_a70 = data2_a60["lambda"]
	
	flux_2_a1 = data2_a1["tempflux"]
	flux_2_a2 = data2_a2["tempflux"]
	flux_2_a3 = data2_a3["tempflux"]
	flux_2_a4 = data2_a4["tempflux"]
	flux_2_a5 = data2_a5["tempflux"]
	flux_2_a6 = data2_a6["tempflux"]
	flux_2_a7 = data2_a7["tempflux"]
	flux_2_a8 = data2_a8["tempflux"]
	flux_2_a9 = data2_a9["tempflux"]
	flux_2_a10 = data2_a10["tempflux"]
	flux_2_a11 = data2_a11["tempflux"]
	flux_2_a12 = data2_a12["tempflux"]
	flux_2_a13 = data2_a13["tempflux"]
	flux_2_a14 = data2_a14["tempflux"]
	flux_2_a15 = data2_a15["tempflux"]
	flux_2_a16 = data2_a16["tempflux"]
	flux_2_a17 = data2_a17["tempflux"]
	flux_2_a18 = data2_a18["tempflux"]
	flux_2_a19 = data2_a19["tempflux"]
	flux_2_a20 = data2_a20["tempflux"]
	flux_2_a21 = data2_a21["tempflux"]
	flux_2_a22 = data2_a22["tempflux"]
	flux_2_a23 = data2_a23["tempflux"]
	flux_2_a24 = data2_a24["tempflux"]
	flux_2_a25 = data2_a25["tempflux"]
	flux_2_a26 = data2_a26["tempflux"]
	flux_2_a27 = data2_a27["tempflux"]
	flux_2_a28 = data2_a28["tempflux"]
	flux_2_a29 = data2_a29["tempflux"]
	flux_2_a30 = data2_a30["tempflux"]
	flux_2_a31 = data2_a31["tempflux"]
	flux_2_a32 = data2_a32["tempflux"]
	flux_2_a33 = data2_a33["tempflux"]
	flux_2_a34 = data2_a34["tempflux"]
	flux_2_a35 = data2_a35["tempflux"]
	flux_2_a36 = data2_a36["tempflux"]
	flux_2_a37 = data2_a37["tempflux"]
	flux_2_a38 = data2_a38["tempflux"]
	flux_2_a39 = data2_a39["tempflux"]
	flux_2_a40 = data2_a40["tempflux"]
	flux_2_a41 = data2_a41["tempflux"]
	flux_2_a42 = data2_a42["tempflux"]
	flux_2_a43 = data2_a43["tempflux"]
	flux_2_a44 = data2_a44["tempflux"]
	flux_2_a45 = data2_a45["tempflux"]
	flux_2_a46 = data2_a46["tempflux"]
	flux_2_a47 = data2_a47["tempflux"]
	flux_2_a48 = data2_a48["tempflux"]
	flux_2_a49 = data2_a49["tempflux"]
	flux_2_a50 = data2_a50["tempflux"]
	flux_2_a51 = data2_a51["tempflux"]
	flux_2_a52 = data2_a52["tempflux"]
	flux_2_a53 = data2_a53["tempflux"]
	flux_2_a54 = data2_a54["tempflux"]
	flux_2_a55 = data2_a55["tempflux"]
	flux_2_a56 = data2_a56["tempflux"]
	flux_2_a57 = data2_a57["tempflux"]
	flux_2_a58 = data2_a58["tempflux"]
	flux_2_a59 = data2_a59["tempflux"]
	flux_2_a60 = data2_a60["tempflux"]
	flux_2_a61 = data2_a51["tempflux"]
	flux_2_a62 = data2_a52["tempflux"]
	flux_2_a63 = data2_a53["tempflux"]
	flux_2_a64 = data2_a54["tempflux"]
	flux_2_a65 = data2_a55["tempflux"]
	flux_2_a66 = data2_a56["tempflux"]
	flux_2_a67 = data2_a57["tempflux"]
	flux_2_a68 = data2_a58["tempflux"]
	flux_2_a69 = data2_a59["tempflux"]
	flux_2_a70 = data2_a60["tempflux"]
	
	flux2_a1 = (flux_2_a1*(lam2_a1**-2.0))*factor
	flux2_a2 = (flux_2_a2*(lam2_a2**-2.0))*factor
	flux2_a3 = (flux_2_a3*(lam2_a3**-2.0))*factor
	flux2_a4 = (flux_2_a4*(lam2_a4**-2.0))*factor
	flux2_a5 = (flux_2_a5*(lam2_a5**-2.0))*factor
	flux2_a6 = (flux_2_a6*(lam2_a6**-2.0))*factor
	flux2_a7 = (flux_2_a7*(lam2_a7**-2.0))*factor
	flux2_a8 = (flux_2_a8*(lam2_a8**-2.0))*factor
	flux2_a9 = (flux_2_a9*(lam2_a9**-2.0))*factor
	flux2_a10 = (flux_2_a10*(lam2_a10**-2.0))*factor
	flux2_a11 = (flux_2_a11*(lam2_a11**-2.0))*factor
	flux2_a12 = (flux_2_a12*(lam2_a12**-2.0))*factor
	flux2_a13 = (flux_2_a13*(lam2_a13**-2.0))*factor
	flux2_a14 = (flux_2_a14*(lam2_a14**-2.0))*factor
	flux2_a15 = (flux_2_a15*(lam2_a15**-2.0))*factor
	flux2_a16 = (flux_2_a16*(lam2_a16**-2.0))*factor
	flux2_a17 = (flux_2_a17*(lam2_a17**-2.0))*factor
	flux2_a18 = (flux_2_a18*(lam2_a18**-2.0))*factor
	flux2_a19 = (flux_2_a19*(lam2_a19**-2.0))*factor
	flux2_a20 = (flux_2_a20*(lam2_a20**-2.0))*factor
	flux2_a21 = (flux_2_a21*(lam2_a21**-2.0))*factor
	flux2_a22 = (flux_2_a22*(lam2_a22**-2.0))*factor
	flux2_a23 = (flux_2_a23*(lam2_a23**-2.0))*factor
	flux2_a24 = (flux_2_a24*(lam2_a24**-2.0))*factor
	flux2_a25 = (flux_2_a25*(lam2_a25**-2.0))*factor
	flux2_a26 = (flux_2_a26*(lam2_a26**-2.0))*factor
	flux2_a27 = (flux_2_a27*(lam2_a27**-2.0))*factor
	flux2_a28 = (flux_2_a28*(lam2_a28**-2.0))*factor
	flux2_a29 = (flux_2_a29*(lam2_a29**-2.0))*factor
	flux2_a30 = (flux_2_a30*(lam2_a30**-2.0))*factor
	flux2_a31 = (flux_2_a31*(lam2_a31**-2.0))*factor
	flux2_a32 = (flux_2_a32*(lam2_a32**-2.0))*factor
	flux2_a33 = (flux_2_a33*(lam2_a33**-2.0))*factor
	flux2_a34 = (flux_2_a34*(lam2_a34**-2.0))*factor
	flux2_a35 = (flux_2_a35*(lam2_a35**-2.0))*factor
	flux2_a36 = (flux_2_a36*(lam2_a36**-2.0))*factor
	flux2_a37 = (flux_2_a37*(lam2_a37**-2.0))*factor
	flux2_a38 = (flux_2_a38*(lam2_a38**-2.0))*factor
	flux2_a39 = (flux_2_a39*(lam2_a39**-2.0))*factor
	flux2_a40 = (flux_2_a40*(lam2_a40**-2.0))*factor
	flux2_a41 = (flux_2_a41*(lam2_a41**-2.0))*factor
	flux2_a42 = (flux_2_a42*(lam2_a42**-2.0))*factor
	flux2_a43 = (flux_2_a43*(lam2_a43**-2.0))*factor
	flux2_a44 = (flux_2_a44*(lam2_a44**-2.0))*factor
	flux2_a45 = (flux_2_a45*(lam2_a45**-2.0))*factor
	flux2_a46 = (flux_2_a46*(lam2_a46**-2.0))*factor
	flux2_a47 = (flux_2_a47*(lam2_a47**-2.0))*factor
	flux2_a48 = (flux_2_a48*(lam2_a48**-2.0))*factor
	flux2_a49 = (flux_2_a49*(lam2_a49**-2.0))*factor
	flux2_a50 = (flux_2_a50*(lam2_a50**-2.0))*factor
	flux2_a51 = (flux_2_a51*(lam2_a51**-2.0))*factor
	flux2_a52 = (flux_2_a52*(lam2_a52**-2.0))*factor
	flux2_a53 = (flux_2_a53*(lam2_a53**-2.0))*factor
	flux2_a54 = (flux_2_a54*(lam2_a54**-2.0))*factor
	flux2_a55 = (flux_2_a55*(lam2_a55**-2.0))*factor
	flux2_a56 = (flux_2_a56*(lam2_a56**-2.0))*factor
	flux2_a57 = (flux_2_a57*(lam2_a57**-2.0))*factor
	flux2_a58 = (flux_2_a58*(lam2_a58**-2.0))*factor
	flux2_a59 = (flux_2_a59*(lam2_a59**-2.0))*factor
	flux2_a60 = (flux_2_a60*(lam2_a60**-2.0))*factor
	flux2_a61 = (flux_2_a51*(lam2_a51**-2.0))*factor
	flux2_a62 = (flux_2_a52*(lam2_a52**-2.0))*factor
	flux2_a63 = (flux_2_a53*(lam2_a53**-2.0))*factor
	flux2_a64 = (flux_2_a54*(lam2_a54**-2.0))*factor
	flux2_a65 = (flux_2_a55*(lam2_a55**-2.0))*factor
	flux2_a66 = (flux_2_a56*(lam2_a56**-2.0))*factor
	flux2_a67 = (flux_2_a57*(lam2_a57**-2.0))*factor
	flux2_a68 = (flux_2_a58*(lam2_a58**-2.0))*factor
	flux2_a69 = (flux_2_a59*(lam2_a59**-2.0))*factor
	flux2_a70 = (flux_2_a60*(lam2_a60**-2.0))*factor
	
	data3_a1 = ascii.read("531..temp_sed")
	data3_a2 = ascii.read("1606..temp_sed")
	data3_a3 = ascii.read("2578..temp_sed")
	data3_a4 = ascii.read("2918..temp_sed")
	data3_a5 = ascii.read("2957..temp_sed")
	data3_a6 = ascii.read("8635..temp_sed")
	data3_a7 = ascii.read("9128..temp_sed")
	data3_a8 = ascii.read("9870..temp_sed")
	data3_a9 = ascii.read("10755..temp_sed")
	data3_a10 = ascii.read("10893..temp_sed")
	data3_a11 = ascii.read("11416..temp_sed")
	data3_a12 = ascii.read("11773..temp_sed")
	data3_a13 = ascii.read("12227..temp_sed")
	data3_a14 = ascii.read("12479..temp_sed")
	data3_a15 = ascii.read("15709..temp_sed")
	data3_a16 = ascii.read("15871..temp_sed")
	data3_a17 = ascii.read("16065..temp_sed")
	data3_a18 = ascii.read("17754..temp_sed")
	data3_a19 = ascii.read("23040..temp_sed")
	data3_a20 = ascii.read("23645..temp_sed")
	data3_a21 = ascii.read("25300..temp_sed")
	data3_a22 = ascii.read("26508..temp_sed")
	data3_a23 = ascii.read("26952..temp_sed")
	data3_a24 = ascii.read("27802..temp_sed")
	data3_a25 = ascii.read("29106..temp_sed")
	data3_a26 = ascii.read("29178..temp_sed")
	data3_a27 = ascii.read("29861..temp_sed")
	data3_a28 = ascii.read("29987..temp_sed")
	data3_a29 = ascii.read("30967..temp_sed")
	data3_a30 = ascii.read("31353..temp_sed")
	data3_a31 = ascii.read("32014..temp_sed")
	data3_a32 = ascii.read("32686..temp_sed")
	data3_a33 = ascii.read("33863..temp_sed")
	data3_a34 = ascii.read("33925..temp_sed")
	data3_a35 = ascii.read("34685..temp_sed")
	data3_a36 = ascii.read("34918..temp_sed")
	data3_a37 = ascii.read("36104..temp_sed")
	data3_a38 = ascii.read("36574..temp_sed")
	data3_a39 = ascii.read("37853..temp_sed")
	data3_a40 = ascii.read("38187..temp_sed")
	
	
	lam3_a1 = data3_a1["lambda"]
	lam3_a2 = data3_a2["lambda"]
	lam3_a3 = data3_a3["lambda"]
	lam3_a4 = data3_a4["lambda"]
	lam3_a5 = data3_a5["lambda"]
	lam3_a6 = data3_a6["lambda"]
	lam3_a7 = data3_a7["lambda"]
	lam3_a8 = data3_a8["lambda"]
	lam3_a9 = data3_a9["lambda"]
	lam3_a10 = data3_a10["lambda"]
	lam3_a11 = data3_a11["lambda"]
	lam3_a12 = data3_a12["lambda"]
	lam3_a13 = data3_a13["lambda"]
	lam3_a14 = data3_a14["lambda"]
	lam3_a15 = data3_a15["lambda"]
	lam3_a16 = data3_a16["lambda"]
	lam3_a17 = data3_a17["lambda"]
	lam3_a18 = data3_a18["lambda"]
	lam3_a19 = data3_a19["lambda"]
	lam3_a20 = data3_a20["lambda"]
	lam3_a21 = data3_a21["lambda"]
	lam3_a22 = data3_a22["lambda"]
	lam3_a23 = data3_a23["lambda"]
	lam3_a24 = data3_a24["lambda"]
	lam3_a25 = data3_a25["lambda"]
	lam3_a26 = data3_a26["lambda"]
	lam3_a27 = data3_a27["lambda"]
	lam3_a28 = data3_a28["lambda"]
	lam3_a29 = data3_a29["lambda"]
	lam3_a30 = data3_a30["lambda"]
	lam3_a31 = data3_a31["lambda"]
	lam3_a32 = data3_a32["lambda"]
	lam3_a33 = data3_a33["lambda"]
	lam3_a34 = data3_a34["lambda"]
	lam3_a35 = data3_a35["lambda"]
	lam3_a36 = data3_a36["lambda"]
	lam3_a37 = data3_a37["lambda"]
	lam3_a38 = data3_a38["lambda"]
	lam3_a39 = data3_a39["lambda"]
	lam3_a40 = data3_a40["lambda"]
	
	
	flux_3_a1 = data3_a1["tempflux"]
	flux_3_a2 = data3_a2["tempflux"]
	flux_3_a3 = data3_a3["tempflux"]
	flux_3_a4 = data3_a4["tempflux"]
	flux_3_a5 = data3_a5["tempflux"]
	flux_3_a6 = data3_a6["tempflux"]
	flux_3_a7 = data3_a7["tempflux"]
	flux_3_a8 = data3_a8["tempflux"]
	flux_3_a9 = data3_a9["tempflux"]
	flux_3_a10 = data3_a10["tempflux"]
	flux_3_a11 = data3_a11["tempflux"]
	flux_3_a12 = data3_a12["tempflux"]
	flux_3_a13 = data3_a13["tempflux"]
	flux_3_a14 = data3_a14["tempflux"]
	flux_3_a15 = data3_a15["tempflux"]
	flux_3_a16 = data3_a16["tempflux"]
	flux_3_a17 = data3_a17["tempflux"]
	flux_3_a18 = data3_a18["tempflux"]
	flux_3_a19 = data3_a19["tempflux"]
	flux_3_a20 = data3_a20["tempflux"]
	flux_3_a21 = data3_a21["tempflux"]
	flux_3_a22 = data3_a22["tempflux"]
	flux_3_a23 = data3_a23["tempflux"]
	flux_3_a24 = data3_a24["tempflux"]
	flux_3_a25 = data3_a25["tempflux"]
	flux_3_a26 = data3_a26["tempflux"]
	flux_3_a27 = data3_a27["tempflux"]
	flux_3_a28 = data3_a28["tempflux"]
	flux_3_a29 = data3_a29["tempflux"]
	flux_3_a30 = data3_a30["tempflux"]
	flux_3_a31 = data3_a31["tempflux"]
	flux_3_a32 = data3_a32["tempflux"]
	flux_3_a33 = data3_a33["tempflux"]
	flux_3_a34 = data3_a34["tempflux"]
	flux_3_a35 = data3_a35["tempflux"]
	flux_3_a36 = data3_a36["tempflux"]
	flux_3_a37 = data3_a37["tempflux"]
	flux_3_a38 = data3_a38["tempflux"]
	flux_3_a39 = data3_a39["tempflux"]
	flux_3_a40 = data3_a40["tempflux"]
	
	
	flux3_a1 = (flux_3_a1*(lam3_a1**-2.0))*factor
	flux3_a2 = (flux_3_a2*(lam3_a2**-2.0))*factor
	flux3_a3 = (flux_3_a3*(lam3_a3**-2.0))*factor
	flux3_a4 = (flux_3_a4*(lam3_a4**-2.0))*factor
	flux3_a5 = (flux_3_a5*(lam3_a5**-2.0))*factor
	flux3_a6 = (flux_3_a6*(lam3_a6**-2.0))*factor
	flux3_a7 = (flux_3_a7*(lam3_a7**-2.0))*factor
	flux3_a8 = (flux_3_a8*(lam3_a8**-2.0))*factor
	flux3_a9 = (flux_3_a9*(lam3_a9**-2.0))*factor
	flux3_a10 = (flux_3_a10*(lam3_a10**-2.0))*factor
	flux3_a11 = (flux_3_a11*(lam3_a11**-2.0))*factor
	flux3_a12 = (flux_3_a12*(lam3_a12**-2.0))*factor
	flux3_a13 = (flux_3_a13*(lam3_a13**-2.0))*factor
	flux3_a14 = (flux_3_a14*(lam3_a14**-2.0))*factor
	flux3_a15 = (flux_3_a15*(lam3_a15**-2.0))*factor
	flux3_a16 = (flux_3_a16*(lam3_a16**-2.0))*factor
	flux3_a17 = (flux_3_a17*(lam3_a17**-2.0))*factor
	flux3_a18 = (flux_3_a18*(lam3_a18**-2.0))*factor
	flux3_a19 = (flux_3_a19*(lam3_a19**-2.0))*factor
	flux3_a20 = (flux_3_a20*(lam3_a20**-2.0))*factor
	flux3_a21 = (flux_3_a21*(lam3_a21**-2.0))*factor
	flux3_a22 = (flux_3_a22*(lam3_a22**-2.0))*factor
	flux3_a23 = (flux_3_a23*(lam3_a23**-2.0))*factor
	flux3_a24 = (flux_3_a24*(lam3_a24**-2.0))*factor
	flux3_a25 = (flux_3_a25*(lam3_a25**-2.0))*factor
	flux3_a26 = (flux_3_a26*(lam3_a26**-2.0))*factor
	flux3_a27 = (flux_3_a27*(lam3_a27**-2.0))*factor
	flux3_a28 = (flux_3_a28*(lam3_a28**-2.0))*factor
	flux3_a29 = (flux_3_a29*(lam3_a29**-2.0))*factor
	flux3_a30 = (flux_3_a30*(lam3_a30**-2.0))*factor
	flux3_a31 = (flux_3_a31*(lam3_a31**-2.0))*factor
	flux3_a32 = (flux_3_a32*(lam3_a32**-2.0))*factor
	flux3_a33 = (flux_3_a33*(lam3_a33**-2.0))*factor
	flux3_a34 = (flux_3_a34*(lam3_a34**-2.0))*factor
	flux3_a35 = (flux_3_a35*(lam3_a35**-2.0))*factor
	flux3_a36 = (flux_3_a36*(lam3_a36**-2.0))*factor
	flux3_a37 = (flux_3_a37*(lam3_a37**-2.0))*factor
	flux3_a38 = (flux_3_a38*(lam3_a38**-2.0))*factor
	flux3_a39 = (flux_3_a39*(lam3_a39**-2.0))*factor
	flux3_a40 = (flux_3_a40*(lam3_a40**-2.0))*factor
	
	
	
	os.chdir("/Volumes/TOSHIBA EXT/3d_hst/noah_massive/cosmos_massive")
	
	data1_c1 = ascii.read("796..temp_sed")
	data1_c2 = ascii.read("2348..temp_sed")
	data1_c3 = ascii.read("5238..temp_sed")
	data1_c4 = ascii.read("9111..temp_sed")
	data1_c5 = ascii.read("10703..temp_sed")
	data1_c6 = ascii.read("11783..temp_sed")
	data1_c7 = ascii.read("11871..temp_sed")
	data1_c8 = ascii.read("12699..temp_sed")
	data1_c9 = ascii.read("12767..temp_sed")
	data1_c10 = ascii.read("13206..temp_sed")
	data1_c11 = ascii.read("13890..temp_sed")
	data1_c12 = ascii.read("15066..temp_sed")
	data1_c13 = ascii.read("17263..temp_sed")
	data1_c14 = ascii.read("18575..temp_sed")
	data1_c15 = ascii.read("25627..temp_sed")
	data1_c16 = ascii.read("27769..temp_sed")
	data1_c17 = ascii.read("31555..temp_sed")
	data1_c18 = ascii.read("32549..temp_sed")
	
	
	lam1_c1 = data1_c1["lambda"]
	lam1_c2 = data1_c2["lambda"]
	lam1_c3 = data1_c3["lambda"]
	lam1_c4 = data1_c4["lambda"]
	lam1_c5 = data1_c5["lambda"]
	lam1_c6 = data1_c6["lambda"]
	lam1_c7 = data1_c7["lambda"]
	lam1_c8 = data1_c8["lambda"]
	lam1_c9 = data1_c9["lambda"]
	lam1_c10 = data1_c10["lambda"]
	lam1_c11 = data1_c11["lambda"]
	lam1_c12 = data1_c12["lambda"]
	lam1_c13 = data1_c13["lambda"]
	lam1_c14 = data1_c14["lambda"]
	lam1_c15 = data1_c15["lambda"]
	lam1_c16 = data1_c16["lambda"]
	lam1_c17 = data1_c17["lambda"]
	lam1_c18 = data1_c18["lambda"]
	
	
	flux_1_c1 = data1_c1["tempflux"]
	flux_1_c2 = data1_c2["tempflux"]
	flux_1_c3 = data1_c3["tempflux"]
	flux_1_c4 = data1_c4["tempflux"]
	flux_1_c5 = data1_c5["tempflux"]
	flux_1_c6 = data1_c6["tempflux"]
	flux_1_c7 = data1_c7["tempflux"]
	flux_1_c8 = data1_c8["tempflux"]
	flux_1_c9 = data1_c9["tempflux"]
	flux_1_c10 = data1_c10["tempflux"]
	flux_1_c11 = data1_c11["tempflux"]
	flux_1_c12 = data1_c12["tempflux"]
	flux_1_c13 = data1_c13["tempflux"]
	flux_1_c14 = data1_c14["tempflux"]
	flux_1_c15 = data1_c15["tempflux"]
	flux_1_c16 = data1_c16["tempflux"]
	flux_1_c17 = data1_c17["tempflux"]
	flux_1_c18 = data1_c18["tempflux"]
	
	
	flux1_c1 = (flux_1_c1*(lam1_c1**-2.0))*factor
	flux1_c2 = (flux_1_c2*(lam1_c2**-2.0))*factor
	flux1_c3 = (flux_1_c3*(lam1_c3**-2.0))*factor
	flux1_c4 = (flux_1_c4*(lam1_c4**-2.0))*factor
	flux1_c5 = (flux_1_c5*(lam1_c5**-2.0))*factor
	flux1_c6 = (flux_1_c6*(lam1_c6**-2.0))*factor
	flux1_c7 = (flux_1_c7*(lam1_c7**-2.0))*factor
	flux1_c8 = (flux_1_c8*(lam1_c8**-2.0))*factor
	flux1_c9 = (flux_1_c9*(lam1_c9**-2.0))*factor
	flux1_c10 = (flux_1_c10*(lam1_c10**-2.0))*factor
	flux1_c11 = (flux_1_c11*(lam1_c11**-2.0))*factor
	flux1_c12 = (flux_1_c12*(lam1_c12**-2.0))*factor
	flux1_c13 = (flux_1_c13*(lam1_c13**-2.0))*factor
	flux1_c14 = (flux_1_c14*(lam1_c14**-2.0))*factor
	flux1_c15 = (flux_1_c15*(lam1_c15**-2.0))*factor
	flux1_c16 = (flux_1_c16*(lam1_c16**-2.0))*factor
	flux1_c17 = (flux_1_c17*(lam1_c17**-2.0))*factor
	flux1_c18 = (flux_1_c18*(lam1_c18**-2.0))*factor
	

	data2_c1 = ascii.read("312..temp_sed")
	data2_c2 = ascii.read("363..temp_sed")
	data2_c3 = ascii.read("728..temp_sed")
	data2_c4 = ascii.read("2616..temp_sed")
	data2_c5 = ascii.read("2816..temp_sed")
	data2_c6 = ascii.read("3200..temp_sed")
	data2_c7 = ascii.read("4536..temp_sed")
	data2_c8 = ascii.read("7216..temp_sed")
	data2_c9 = ascii.read("7411..temp_sed")
	data2_c10 = ascii.read("9667..temp_sed")
	data2_c11 = ascii.read("10128..temp_sed")
	data2_c12 = ascii.read("10592..temp_sed")
	data2_c13 = ascii.read("10989..temp_sed")
	data2_c14 = ascii.read("11973..temp_sed")
	data2_c15 = ascii.read("17089..temp_sed")
	data2_c16 = ascii.read("17406..temp_sed")
	data2_c17 = ascii.read("18688..temp_sed")
	data2_c18 = ascii.read("20983..temp_sed")
	data2_c19 = ascii.read("21723..temp_sed")
	data2_c20 = ascii.read("24462..temp_sed")
	data2_c21 = ascii.read("25534..temp_sed")
	data2_c22 = ascii.read("25581..temp_sed")
	data2_c23 = ascii.read("28492..temp_sed")
	data2_c24 = ascii.read("28523..temp_sed")
	data2_c25 = ascii.read("29222..temp_sed")
	data2_c26 = ascii.read("31090..temp_sed")
	
	
	lam2_c1 = data2_c1["lambda"]
	lam2_c2 = data2_c2["lambda"]
	lam2_c3 = data2_c3["lambda"]
	lam2_c4 = data2_c4["lambda"]
	lam2_c5 = data2_c5["lambda"]
	lam2_c6 = data2_c6["lambda"]
	lam2_c7 = data2_c7["lambda"]
	lam2_c8 = data2_c8["lambda"]
	lam2_c9 = data2_c9["lambda"]
	lam2_c10 = data2_c10["lambda"]
	lam2_c11 = data2_c11["lambda"]
	lam2_c12 = data2_c12["lambda"]
	lam2_c13 = data2_c13["lambda"]
	lam2_c14 = data2_c14["lambda"]
	lam2_c15 = data2_c15["lambda"]
	lam2_c16 = data2_c16["lambda"]
	lam2_c17 = data2_c17["lambda"]
	lam2_c18 = data2_c18["lambda"]
	lam2_c19 = data2_c19["lambda"]
	lam2_c20 = data2_c20["lambda"]
	lam2_c21 = data2_c21["lambda"]
	lam2_c22 = data2_c22["lambda"]
	lam2_c23 = data2_c23["lambda"]
	lam2_c24 = data2_c24["lambda"]
	lam2_c25 = data2_c25["lambda"]
	lam2_c26 = data2_c26["lambda"]
	
	
	flux_2_c1 = data2_c1["tempflux"]
	flux_2_c2 = data2_c2["tempflux"]
	flux_2_c3 = data2_c3["tempflux"]
	flux_2_c4 = data2_c4["tempflux"]
	flux_2_c5 = data2_c5["tempflux"]
	flux_2_c6 = data2_c6["tempflux"]
	flux_2_c7 = data2_c7["tempflux"]
	flux_2_c8 = data2_c8["tempflux"]
	flux_2_c9 = data2_c9["tempflux"]
	flux_2_c10 = data2_c10["tempflux"]
	flux_2_c11 = data2_c11["tempflux"]
	flux_2_c12 = data2_c12["tempflux"]
	flux_2_c13 = data2_c13["tempflux"]
	flux_2_c14 = data2_c14["tempflux"]
	flux_2_c15 = data2_c15["tempflux"]
	flux_2_c16 = data2_c16["tempflux"]
	flux_2_c17 = data2_c17["tempflux"]
	flux_2_c18 = data2_c18["tempflux"]
	flux_2_c19 = data2_c19["tempflux"]
	flux_2_c20 = data2_c20["tempflux"]
	flux_2_c21 = data2_c21["tempflux"]
	flux_2_c22 = data2_c22["tempflux"]
	flux_2_c23 = data2_c23["tempflux"]
	flux_2_c24 = data2_c24["tempflux"]
	flux_2_c25 = data2_c25["tempflux"]
	flux_2_c26 = data2_c26["tempflux"]
	
	
	flux2_c1 = (flux_2_c1*(lam2_c1**-2.0))*factor
	flux2_c2 = (flux_2_c2*(lam2_c2**-2.0))*factor
	flux2_c3 = (flux_2_c3*(lam2_c3**-2.0))*factor
	flux2_c4 = (flux_2_c4*(lam2_c4**-2.0))*factor
	flux2_c5 = (flux_2_c5*(lam2_c5**-2.0))*factor
	flux2_c6 = (flux_2_c6*(lam2_c6**-2.0))*factor
	flux2_c7 = (flux_2_c7*(lam2_c7**-2.0))*factor
	flux2_c8 = (flux_2_c8*(lam2_c8**-2.0))*factor
	flux2_c9 = (flux_2_c9*(lam2_c9**-2.0))*factor
	flux2_c10 = (flux_2_c10*(lam2_c10**-2.0))*factor
	flux2_c11 = (flux_2_c11*(lam2_c11**-2.0))*factor
	flux2_c12 = (flux_2_c12*(lam2_c12**-2.0))*factor
	flux2_c13 = (flux_2_c13*(lam2_c13**-2.0))*factor
	flux2_c14 = (flux_2_c14*(lam2_c14**-2.0))*factor
	flux2_c15 = (flux_2_c15*(lam2_c15**-2.0))*factor
	flux2_c16 = (flux_2_c16*(lam2_c16**-2.0))*factor
	flux2_c17 = (flux_2_c17*(lam2_c17**-2.0))*factor
	flux2_c18 = (flux_2_c18*(lam2_c18**-2.0))*factor
	flux2_c19 = (flux_2_c19*(lam2_c19**-2.0))*factor
	flux2_c20 = (flux_2_c20*(lam2_c20**-2.0))*factor
	flux2_c21 = (flux_2_c21*(lam2_c21**-2.0))*factor
	flux2_c22 = (flux_2_c22*(lam2_c22**-2.0))*factor
	flux2_c23 = (flux_2_c23*(lam2_c23**-2.0))*factor
	flux2_c24 = (flux_2_c24*(lam2_c24**-2.0))*factor
	flux2_c25 = (flux_2_c25*(lam2_c25**-2.0))*factor
	flux2_c26 = (flux_2_c26*(lam2_c26**-2.0))*factor
	
	
	data3_c1 = ascii.read("490..temp_sed")
	data3_c2 = ascii.read("1769..temp_sed")
	data3_c3 = ascii.read("2049..temp_sed")
	data3_c4 = ascii.read("3182..temp_sed")
	data3_c5 = ascii.read("3206..temp_sed")
	data3_c6 = ascii.read("5473..temp_sed")
	data3_c7 = ascii.read("5530..temp_sed")
	data3_c8 = ascii.read("6159..temp_sed")
	data3_c9 = ascii.read("7884..temp_sed")
	data3_c10 = ascii.read("7951..temp_sed")
	data3_c11 = ascii.read("9871..temp_sed")
	data3_c12 = ascii.read("11314..temp_sed")
	data3_c13 = ascii.read("11337..temp_sed")
	data3_c14 = ascii.read("11363..temp_sed")
	data3_c15 = ascii.read("11494..temp_sed")
	data3_c16 = ascii.read("12020..temp_sed")
	data3_c17 = ascii.read("12995..temp_sed")
	data3_c18 = ascii.read("13083..temp_sed")
	data3_c19 = ascii.read("13174..temp_sed")
	data3_c20 = ascii.read("16419..temp_sed")
	data3_c21 = ascii.read("19090..temp_sed")
	data3_c22 = ascii.read("19153..temp_sed")
	data3_c23 = ascii.read("20668..temp_sed")
	data3_c24 = ascii.read("22995..temp_sed")
	data3_c25 = ascii.read("23021..temp_sed")
	data3_c26 = ascii.read("23673..temp_sed")
	data3_c27 = ascii.read("25515..temp_sed")
	data3_c28 = ascii.read("26039..temp_sed")
	data3_c29 = ascii.read("26338..temp_sed")
	data3_c30 = ascii.read("26957..temp_sed")
	data3_c31 = ascii.read("27289..temp_sed")
	data3_c32 = ascii.read("28344..temp_sed")
	data3_c33 = ascii.read("28565..temp_sed")
	data3_c34 = ascii.read("31922..temp_sed")
	data3_c35 = ascii.read("33199..temp_sed")
	
	
	lam3_c1 = data3_c1["lambda"]
	lam3_c2 = data3_c2["lambda"]
	lam3_c3 = data3_c3["lambda"]
	lam3_c4 = data3_c4["lambda"]
	lam3_c5 = data3_c5["lambda"]
	lam3_c6 = data3_c6["lambda"]
	lam3_c7 = data3_c7["lambda"]
	lam3_c8 = data3_c8["lambda"]
	lam3_c9 = data3_c9["lambda"]
	lam3_c10 = data3_c10["lambda"]
	lam3_c11 = data3_c11["lambda"]
	lam3_c12 = data3_c12["lambda"]
	lam3_c13 = data3_c13["lambda"]
	lam3_c14 = data3_c14["lambda"]
	lam3_c15 = data3_c15["lambda"]
	lam3_c16 = data3_c16["lambda"]
	lam3_c17 = data3_c17["lambda"]
	lam3_c18 = data3_c18["lambda"]
	lam3_c19 = data3_c19["lambda"]
	lam3_c20 = data3_c20["lambda"]
	lam3_c21 = data3_c21["lambda"]
	lam3_c22 = data3_c22["lambda"]
	lam3_c23 = data3_c23["lambda"]
	lam3_c24 = data3_c24["lambda"]
	lam3_c25 = data3_c25["lambda"]
	lam3_c26 = data3_c26["lambda"]
	lam3_c27 = data3_c27["lambda"]
	lam3_c28 = data3_c28["lambda"]
	lam3_c29 = data3_c29["lambda"]
	lam3_c30 = data3_c30["lambda"]
	lam3_c31 = data3_c31["lambda"]
	lam3_c32 = data3_c32["lambda"]
	lam3_c33 = data3_c33["lambda"]
	lam3_c34 = data3_c34["lambda"]
	lam3_c35 = data3_c35["lambda"]
	
	
	flux_3_c1 = data3_c1["tempflux"]
	flux_3_c2 = data3_c2["tempflux"]
	flux_3_c3 = data3_c3["tempflux"]
	flux_3_c4 = data3_c4["tempflux"]
	flux_3_c5 = data3_c5["tempflux"]
	flux_3_c6 = data3_c6["tempflux"]
	flux_3_c7 = data3_c7["tempflux"]
	flux_3_c8 = data3_c8["tempflux"]
	flux_3_c9 = data3_c9["tempflux"]
	flux_3_c10 = data3_c10["tempflux"]
	flux_3_c11 = data3_c11["tempflux"]
	flux_3_c12 = data3_c12["tempflux"]
	flux_3_c13 = data3_c13["tempflux"]
	flux_3_c14 = data3_c14["tempflux"]
	flux_3_c15 = data3_c15["tempflux"]
	flux_3_c16 = data3_c16["tempflux"]
	flux_3_c17 = data3_c17["tempflux"]
	flux_3_c18 = data3_c18["tempflux"]
	flux_3_c19 = data3_c19["tempflux"]
	flux_3_c20 = data3_c20["tempflux"]
	flux_3_c21 = data3_c21["tempflux"]
	flux_3_c22 = data3_c22["tempflux"]
	flux_3_c23 = data3_c23["tempflux"]
	flux_3_c24 = data3_c24["tempflux"]
	flux_3_c25 = data3_c25["tempflux"]
	flux_3_c26 = data3_c26["tempflux"]
	flux_3_c27 = data3_c27["tempflux"]
	flux_3_c28 = data3_c28["tempflux"]
	flux_3_c29 = data3_c29["tempflux"]
	flux_3_c30 = data3_c30["tempflux"]
	flux_3_c31 = data3_c31["tempflux"]
	flux_3_c32 = data3_c32["tempflux"]
	flux_3_c33 = data3_c33["tempflux"]
	flux_3_c34 = data3_c34["tempflux"]
	flux_3_c35 = data3_c35["tempflux"]
	
	
	flux3_c1 = (flux_3_c1*(lam3_c1**-2.0))*factor
	flux3_c2 = (flux_3_c2*(lam3_c2**-2.0))*factor
	flux3_c3 = (flux_3_c3*(lam3_c3**-2.0))*factor
	flux3_c4 = (flux_3_c4*(lam3_c4**-2.0))*factor
	flux3_c5 = (flux_3_c5*(lam3_c5**-2.0))*factor
	flux3_c6 = (flux_3_c6*(lam3_c6**-2.0))*factor
	flux3_c7 = (flux_3_c7*(lam3_c7**-2.0))*factor
	flux3_c8 = (flux_3_c8*(lam3_c8**-2.0))*factor
	flux3_c9 = (flux_3_c9*(lam3_c9**-2.0))*factor
	flux3_c10 = (flux_3_c10*(lam3_c10**-2.0))*factor
	flux3_c11 = (flux_3_c11*(lam3_c11**-2.0))*factor
	flux3_c12 = (flux_3_c12*(lam3_c12**-2.0))*factor
	flux3_c13 = (flux_3_c13*(lam3_c13**-2.0))*factor
	flux3_c14 = (flux_3_c14*(lam3_c14**-2.0))*factor
	flux3_c15 = (flux_3_c15*(lam3_c15**-2.0))*factor
	flux3_c16 = (flux_3_c16*(lam3_c16**-2.0))*factor
	flux3_c17 = (flux_3_c17*(lam3_c17**-2.0))*factor
	flux3_c18 = (flux_3_c18*(lam3_c18**-2.0))*factor
	flux3_c19 = (flux_3_c19*(lam3_c19**-2.0))*factor
	flux3_c20 = (flux_3_c20*(lam3_c20**-2.0))*factor
	flux3_c21 = (flux_3_c21*(lam3_c21**-2.0))*factor
	flux3_c22 = (flux_3_c22*(lam3_c22**-2.0))*factor
	flux3_c23 = (flux_3_c23*(lam3_c23**-2.0))*factor
	flux3_c24 = (flux_3_c24*(lam3_c24**-2.0))*factor
	flux3_c25 = (flux_3_c25*(lam3_c25**-2.0))*factor
	flux3_c26 = (flux_3_c26*(lam3_c26**-2.0))*factor
	flux3_c27 = (flux_3_c27*(lam3_c27**-2.0))*factor
	flux3_c28 = (flux_3_c28*(lam3_c28**-2.0))*factor
	flux3_c29 = (flux_3_c29*(lam3_c29**-2.0))*factor
	flux3_c30 = (flux_3_c30*(lam3_c30**-2.0))*factor
	flux3_c31 = (flux_3_c31*(lam3_c31**-2.0))*factor
	flux3_c32 = (flux_3_c32*(lam3_c32**-2.0))*factor
	flux3_c33 = (flux_3_c33*(lam3_c33**-2.0))*factor
	flux3_c34 = (flux_3_c34*(lam3_c34**-2.0))*factor
	flux3_c35 = (flux_3_c35*(lam3_c35**-2.0))*factor
	
	
	
	os.chdir("/Volumes/TOSHIBA EXT/3d_hst/noah_massive/goodsn_massive")
	
	data1_n1 = ascii.read("57.temp_sed")
	data1_n2 = ascii.read("128.temp_sed")
	data1_n3 = ascii.read("576.temp_sed")
	data1_n4 = ascii.read("1749.temp_sed")
	data1_n5 = ascii.read("2265.temp_sed")
	data1_n6 = ascii.read("35299.temp_sed")
	data1_n7 = ascii.read("2868.temp_sed")
	data1_n8 = ascii.read("3133.temp_sed")
	data1_n9 = ascii.read("4711.temp_sed")
	data1_n10 = ascii.read("7372.temp_sed")
	data1_n11 = ascii.read("9056.temp_sed")
	data1_n12 = ascii.read("10280.temp_sed")
	data1_n13 = ascii.read("10606.temp_sed")
	data1_n14 = ascii.read("11706.temp_sed")
	data1_n15 = ascii.read("11826.temp_sed")
	data1_n16 = ascii.read("12342.temp_sed")
	data1_n17 = ascii.read("12561.temp_sed")
	data1_n18 = ascii.read("13971.temp_sed")
	data1_n19 = ascii.read("17270.temp_sed")
	data1_n20 = ascii.read("21156.temp_sed")
	data1_n21 = ascii.read("23564.temp_sed")
	data1_n22 = ascii.read("25216.temp_sed")
	data1_n23 = ascii.read("25813.temp_sed")
	data1_n24 = ascii.read("32162.temp_sed")
	data1_n25 = ascii.read("35090.temp_sed")
	
	
	lam1_n1 = data1_n1["lambda"]
	lam1_n2 = data1_n2["lambda"]
	lam1_n3 = data1_n3["lambda"]
	lam1_n4 = data1_n4["lambda"]
	lam1_n5 = data1_n5["lambda"]
	lam1_n6 = data1_n6["lambda"]
	lam1_n7 = data1_n7["lambda"]
	lam1_n8 = data1_n8["lambda"]
	lam1_n9 = data1_n9["lambda"]
	lam1_n10 = data1_n10["lambda"]
	lam1_n11 = data1_n11["lambda"]
	lam1_n12 = data1_n12["lambda"]
	lam1_n13 = data1_n13["lambda"]
	lam1_n14 = data1_n14["lambda"]
	lam1_n15 = data1_n15["lambda"]
	lam1_n16 = data1_n16["lambda"]
	lam1_n17 = data1_n17["lambda"]
	lam1_n18 = data1_n18["lambda"]
	lam1_n19 = data1_n19["lambda"]
	lam1_n20 = data1_n20["lambda"]
	lam1_n21 = data1_n21["lambda"]
	lam1_n22 = data1_n22["lambda"]
	lam1_n23 = data1_n23["lambda"]
	lam1_n24 = data1_n24["lambda"]
	lam1_n25 = data1_n25["lambda"]
	
	
	flux_1_n1 = data1_n1["tempflux"]
	flux_1_n2 = data1_n2["tempflux"]
	flux_1_n3 = data1_n3["tempflux"]
	flux_1_n4 = data1_n4["tempflux"]
	flux_1_n5 = data1_n5["tempflux"]
	flux_1_n6 = data1_n6["tempflux"]
	flux_1_n7 = data1_n7["tempflux"]
	flux_1_n8 = data1_n8["tempflux"]
	flux_1_n9 = data1_n9["tempflux"]
	flux_1_n10 = data1_n10["tempflux"]
	flux_1_n11 = data1_n11["tempflux"]
	flux_1_n12 = data1_n12["tempflux"]
	flux_1_n13 = data1_n13["tempflux"]
	flux_1_n14 = data1_n14["tempflux"]
	flux_1_n15 = data1_n15["tempflux"]
	flux_1_n16 = data1_n16["tempflux"]
	flux_1_n17 = data1_n17["tempflux"]
	flux_1_n18 = data1_n18["tempflux"]
	flux_1_n19 = data1_n19["tempflux"]
	flux_1_n20 = data1_n20["tempflux"]
	flux_1_n21 = data1_n21["tempflux"]
	flux_1_n22 = data1_n22["tempflux"]
	flux_1_n23 = data1_n23["tempflux"]
	flux_1_n24 = data1_n24["tempflux"]
	flux_1_n25 = data1_n25["tempflux"]
	
	
	flux1_n1 = (flux_1_n1*(lam1_n1**-2.0))*factor
	flux1_n2 = (flux_1_n2*(lam1_n2**-2.0))*factor
	flux1_n3 = (flux_1_n3*(lam1_n3**-2.0))*factor
	flux1_n4 = (flux_1_n4*(lam1_n4**-2.0))*factor
	flux1_n5 = (flux_1_n5*(lam1_n5**-2.0))*factor
	flux1_n6 = (flux_1_n6*(lam1_n6**-2.0))*factor
	flux1_n7 = (flux_1_n7*(lam1_n7**-2.0))*factor
	flux1_n8 = (flux_1_n8*(lam1_n8**-2.0))*factor
	flux1_n9 = (flux_1_n9*(lam1_n9**-2.0))*factor
	flux1_n10 = (flux_1_n10*(lam1_n10**-2.0))*factor
	flux1_n11 = (flux_1_n11*(lam1_n11**-2.0))*factor
	flux1_n12 = (flux_1_n12*(lam1_n12**-2.0))*factor
	flux1_n13 = (flux_1_n13*(lam1_n13**-2.0))*factor
	flux1_n14 = (flux_1_n14*(lam1_n14**-2.0))*factor
	flux1_n15 = (flux_1_n15*(lam1_n15**-2.0))*factor
	flux1_n16 = (flux_1_n16*(lam1_n16**-2.0))*factor
	flux1_n17 = (flux_1_n17*(lam1_n17**-2.0))*factor
	flux1_n18 = (flux_1_n18*(lam1_n18**-2.0))*factor
	flux1_n19 = (flux_1_n19*(lam1_n19**-2.0))*factor
	flux1_n20 = (flux_1_n20*(lam1_n20**-2.0))*factor
	flux1_n21 = (flux_1_n21*(lam1_n21**-2.0))*factor
	flux1_n22 = (flux_1_n22*(lam1_n22**-2.0))*factor
	flux1_n23 = (flux_1_n23*(lam1_n23**-2.0))*factor
	flux1_n24 = (flux_1_n24*(lam1_n24**-2.0))*factor
	flux1_n25 = (flux_1_n25*(lam1_n25**-2.0))*factor
	
	
	data2_n1 = ascii.read("1616.temp_sed")
	data2_n2 = ascii.read("3186.temp_sed")
	data2_n3 = ascii.read("4117.temp_sed")
	data2_n4 = ascii.read("5932.temp_sed")
	data2_n5 = ascii.read("9692.temp_sed")
	data2_n6 = ascii.read("10311.temp_sed")
	data2_n7 = ascii.read("14140.temp_sed")
	data2_n8 = ascii.read("14532.temp_sed")
	data2_n9 = ascii.read("16827.temp_sed")
	data2_n10 = ascii.read("18633.temp_sed")
	data2_n11 = ascii.read("19913.temp_sed")
	data2_n12 = ascii.read("20709.temp_sed")
	data2_n13 = ascii.read("23187.temp_sed")
	data2_n14 = ascii.read("23548.temp_sed")
	data2_n15 = ascii.read("25265.temp_sed")
	data2_n16 = ascii.read("29464.temp_sed")
	data2_n17 = ascii.read("32842.temp_sed")
	data2_n18 = ascii.read("33780.temp_sed")
	data2_n19 = ascii.read("35292.temp_sed")
	data2_n20 = ascii.read("36582.temp_sed")
	data2_n21 = ascii.read("37738.temp_sed")
	
	
	lam2_n1 = data2_n1["lambda"]
	lam2_n2 = data2_n2["lambda"]
	lam2_n3 = data2_n3["lambda"]
	lam2_n4 = data2_n4["lambda"]
	lam2_n5 = data2_n5["lambda"]
	lam2_n6 = data2_n6["lambda"]
	lam2_n7 = data2_n7["lambda"]
	lam2_n8 = data2_n8["lambda"]
	lam2_n9 = data2_n9["lambda"]
	lam2_n10 = data2_n10["lambda"]
	lam2_n11 = data2_n11["lambda"]
	lam2_n12 = data2_n12["lambda"]
	lam2_n13 = data2_n13["lambda"]
	lam2_n14 = data2_n14["lambda"]
	lam2_n15 = data2_n15["lambda"]
	lam2_n16 = data2_n16["lambda"]
	lam2_n17 = data2_n17["lambda"]
	lam2_n18 = data2_n18["lambda"]
	lam2_n19 = data2_n19["lambda"]
	lam2_n20 = data2_n20["lambda"]
	lam2_n21 = data2_n21["lambda"]
	
	
	flux_2_n1 = data2_n1["tempflux"]
	flux_2_n2 = data2_n2["tempflux"]
	flux_2_n3 = data2_n3["tempflux"]
	flux_2_n4 = data2_n4["tempflux"]
	flux_2_n5 = data2_n5["tempflux"]
	flux_2_n6 = data2_n6["tempflux"]
	flux_2_n7 = data2_n7["tempflux"]
	flux_2_n8 = data2_n8["tempflux"]
	flux_2_n9 = data2_n9["tempflux"]
	flux_2_n10 = data2_n10["tempflux"]
	flux_2_n11 = data2_n11["tempflux"]
	flux_2_n12 = data2_n12["tempflux"]
	flux_2_n13 = data2_n13["tempflux"]
	flux_2_n14 = data2_n14["tempflux"]
	flux_2_n15 = data2_n15["tempflux"]
	flux_2_n16 = data2_n16["tempflux"]
	flux_2_n17 = data2_n17["tempflux"]
	flux_2_n18 = data2_n18["tempflux"]
	flux_2_n19 = data2_n19["tempflux"]
	flux_2_n20 = data2_n20["tempflux"]
	flux_2_n21 = data2_n21["tempflux"]
	
	
	flux2_n1 = (flux_2_n1*(lam2_n1**-2.0))*factor
	flux2_n2 = (flux_2_n2*(lam2_n2**-2.0))*factor
	flux2_n3 = (flux_2_n3*(lam2_n3**-2.0))*factor
	flux2_n4 = (flux_2_n4*(lam2_n4**-2.0))*factor
	flux2_n5 = (flux_2_n5*(lam2_n5**-2.0))*factor
	flux2_n6 = (flux_2_n6*(lam2_n6**-2.0))*factor
	flux2_n7 = (flux_2_n7*(lam2_n7**-2.0))*factor
	flux2_n8 = (flux_2_n8*(lam2_n8**-2.0))*factor
	flux2_n9 = (flux_2_n9*(lam2_n9**-2.0))*factor
	flux2_n10 = (flux_2_n10*(lam2_n10**-2.0))*factor
	flux2_n11 = (flux_2_n11*(lam2_n11**-2.0))*factor
	flux2_n12 = (flux_2_n12*(lam2_n12**-2.0))*factor
	flux2_n13 = (flux_2_n13*(lam2_n13**-2.0))*factor
	flux2_n14 = (flux_2_n14*(lam2_n14**-2.0))*factor
	flux2_n15 = (flux_2_n15*(lam2_n15**-2.0))*factor
	flux2_n16 = (flux_2_n16*(lam2_n16**-2.0))*factor
	flux2_n17 = (flux_2_n17*(lam2_n17**-2.0))*factor
	flux2_n18 = (flux_2_n18*(lam2_n18**-2.0))*factor
	flux2_n19 = (flux_2_n19*(lam2_n19**-2.0))*factor
	flux2_n20 = (flux_2_n20*(lam2_n20**-2.0))*factor
	flux2_n21 = (flux_2_n21*(lam2_n21**-2.0))*factor
	
	
	data3_n1 = ascii.read("338.temp_sed")
	data3_n2 = ascii.read("764.temp_sed")
	data3_n3 = ascii.read("774.temp_sed")
	data3_n4 = ascii.read("1678.temp_sed")
	data3_n5 = ascii.read("2295.temp_sed")
	data3_n6 = ascii.read("3776.temp_sed")
	data3_n7 = ascii.read("4854.temp_sed")
	data3_n8 = ascii.read("4927.temp_sed")
	data3_n9 = ascii.read("5346.temp_sed")
	data3_n10 = ascii.read("5371.temp_sed")
	data3_n11 = ascii.read("5507.temp_sed")
	data3_n12 = ascii.read("5677.temp_sed")
	data3_n13 = ascii.read("6215.temp_sed")
	data3_n14 = ascii.read("6789.temp_sed")
	data3_n15 = ascii.read("6877.temp_sed")
	data3_n16 = ascii.read("9122.temp_sed")
	data3_n17 = ascii.read("10125.temp_sed")
	data3_n18 = ascii.read("10657.temp_sed")
	data3_n19 = ascii.read("11064.temp_sed")
	data3_n20 = ascii.read("12066.temp_sed")
	data3_n21 = ascii.read("12302.temp_sed")
	data3_n22 = ascii.read("16129.temp_sed")
	data3_n23 = ascii.read("16346.temp_sed")
	data3_n24 = ascii.read("16879.temp_sed")
	data3_n25 = ascii.read("19082.temp_sed")
	data3_n26 = ascii.read("20052.temp_sed")
	data3_n27 = ascii.read("20317.temp_sed")
	data3_n28 = ascii.read("21738.temp_sed")
	data3_n29 = ascii.read("23018.temp_sed")
	data3_n30 = ascii.read("25942.temp_sed")
	data3_n31 = ascii.read("26529.temp_sed")
	data3_n32 = ascii.read("26888.temp_sed")
	data3_n33 = ascii.read("28810.temp_sed")
	data3_n34 = ascii.read("28826.temp_sed")
	data3_n35 = ascii.read("30283.temp_sed")
	data3_n36 = ascii.read("32002.temp_sed")
	data3_n37 = ascii.read("32033.temp_sed")
	data3_n38 = ascii.read("36988.temp_sed")
	
	
	lam3_n1 = data3_n1["lambda"]
	lam3_n2 = data3_n2["lambda"]
	lam3_n3 = data3_n3["lambda"]
	lam3_n4 = data3_n4["lambda"]
	lam3_n5 = data3_n5["lambda"]
	lam3_n6 = data3_n6["lambda"]
	lam3_n7 = data3_n7["lambda"]
	lam3_n8 = data3_n8["lambda"]
	lam3_n9 = data3_n9["lambda"]
	lam3_n10 = data3_n10["lambda"]
	lam3_n11 = data3_n11["lambda"]
	lam3_n12 = data3_n12["lambda"]
	lam3_n13 = data3_n13["lambda"]
	lam3_n14 = data3_n14["lambda"]
	lam3_n15 = data3_n15["lambda"]
	lam3_n16 = data3_n16["lambda"]
	lam3_n17 = data3_n17["lambda"]
	lam3_n18 = data3_n18["lambda"]
	lam3_n19 = data3_n19["lambda"]
	lam3_n20 = data3_n20["lambda"]
	lam3_n21 = data3_n21["lambda"]
	lam3_n22 = data3_n22["lambda"]
	lam3_n23 = data3_n23["lambda"]
	lam3_n24 = data3_n24["lambda"]
	lam3_n25 = data3_n25["lambda"]
	lam3_n26 = data3_n26["lambda"]
	lam3_n27 = data3_n27["lambda"]
	lam3_n28 = data3_n28["lambda"]
	lam3_n29 = data3_n29["lambda"]
	lam3_n30 = data3_n30["lambda"]
	lam3_n31 = data3_n31["lambda"]
	lam3_n32 = data3_n32["lambda"]
	lam3_n33 = data3_n33["lambda"]
	lam3_n34 = data3_n34["lambda"]
	lam3_n35 = data3_n35["lambda"]
	lam3_n36 = data3_n36["lambda"]
	lam3_n37 = data3_n37["lambda"]
	lam3_n38 = data3_n38["lambda"]
	
	
	flux_3_n1 = data3_n1["tempflux"]
	flux_3_n2 = data3_n2["tempflux"]
	flux_3_n3 = data3_n3["tempflux"]
	flux_3_n4 = data3_n4["tempflux"]
	flux_3_n5 = data3_n5["tempflux"]
	flux_3_n6 = data3_n6["tempflux"]
	flux_3_n7 = data3_n7["tempflux"]
	flux_3_n8 = data3_n8["tempflux"]
	flux_3_n9 = data3_n9["tempflux"]
	flux_3_n10 = data3_n10["tempflux"]
	flux_3_n11 = data3_n11["tempflux"]
	flux_3_n12 = data3_n12["tempflux"]
	flux_3_n13 = data3_n13["tempflux"]
	flux_3_n14 = data3_n14["tempflux"]
	flux_3_n15 = data3_n15["tempflux"]
	flux_3_n16 = data3_n16["tempflux"]
	flux_3_n17 = data3_n17["tempflux"]
	flux_3_n18 = data3_n18["tempflux"]
	flux_3_n19 = data3_n19["tempflux"]
	flux_3_n20 = data3_n20["tempflux"]
	flux_3_n21 = data3_n21["tempflux"]
	flux_3_n22 = data3_n22["tempflux"]
	flux_3_n23 = data3_n23["tempflux"]
	flux_3_n24 = data3_n24["tempflux"]
	flux_3_n25 = data3_n25["tempflux"]
	flux_3_n26 = data3_n26["tempflux"]
	flux_3_n27 = data3_n27["tempflux"]
	flux_3_n28 = data3_n28["tempflux"]
	flux_3_n29 = data3_n29["tempflux"]
	flux_3_n30 = data3_n30["tempflux"]
	flux_3_n31 = data3_n31["tempflux"]
	flux_3_n32 = data3_n32["tempflux"]
	flux_3_n33 = data3_n33["tempflux"]
	flux_3_n34 = data3_n34["tempflux"]
	flux_3_n35 = data3_n35["tempflux"]
	flux_3_n36 = data3_n36["tempflux"]
	flux_3_n37 = data3_n37["tempflux"]
	flux_3_n38 = data3_n38["tempflux"]
	
	
	flux3_n1 = (flux_3_n1*(lam3_n1**-2.0))*factor
	flux3_n2 = (flux_3_n2*(lam3_n2**-2.0))*factor
	flux3_n3 = (flux_3_n3*(lam3_n3**-2.0))*factor
	flux3_n4 = (flux_3_n4*(lam3_n4**-2.0))*factor
	flux3_n5 = (flux_3_n5*(lam3_n5**-2.0))*factor
	flux3_n6 = (flux_3_n6*(lam3_n6**-2.0))*factor
	flux3_n7 = (flux_3_n7*(lam3_n7**-2.0))*factor
	flux3_n8 = (flux_3_n8*(lam3_n8**-2.0))*factor
	flux3_n9 = (flux_3_n9*(lam3_n9**-2.0))*factor
	flux3_n10 = (flux_3_n10*(lam3_n10**-2.0))*factor
	flux3_n11 = (flux_3_n11*(lam3_n11**-2.0))*factor
	flux3_n12 = (flux_3_n12*(lam3_n12**-2.0))*factor
	flux3_n13 = (flux_3_n13*(lam3_n13**-2.0))*factor
	flux3_n14 = (flux_3_n14*(lam3_n14**-2.0))*factor
	flux3_n15 = (flux_3_n15*(lam3_n15**-2.0))*factor
	flux3_n16 = (flux_3_n16*(lam3_n16**-2.0))*factor
	flux3_n17 = (flux_3_n17*(lam3_n17**-2.0))*factor
	flux3_n18 = (flux_3_n18*(lam3_n18**-2.0))*factor
	flux3_n19 = (flux_3_n19*(lam3_n19**-2.0))*factor
	flux3_n20 = (flux_3_n20*(lam3_n20**-2.0))*factor
	flux3_n21 = (flux_3_n21*(lam3_n21**-2.0))*factor
	flux3_n22 = (flux_3_n22*(lam3_n22**-2.0))*factor
	flux3_n23 = (flux_3_n23*(lam3_n23**-2.0))*factor
	flux3_n24 = (flux_3_n24*(lam3_n24**-2.0))*factor
	flux3_n25 = (flux_3_n25*(lam3_n25**-2.0))*factor
	flux3_n26 = (flux_3_n26*(lam3_n26**-2.0))*factor
	flux3_n27 = (flux_3_n27*(lam3_n27**-2.0))*factor
	flux3_n28 = (flux_3_n28*(lam3_n28**-2.0))*factor
	flux3_n29 = (flux_3_n29*(lam3_n29**-2.0))*factor
	flux3_n30 = (flux_3_n30*(lam3_n30**-2.0))*factor
	flux3_n31 = (flux_3_n31*(lam3_n31**-2.0))*factor
	flux3_n32 = (flux_3_n32*(lam3_n32**-2.0))*factor
	flux3_n33 = (flux_3_n33*(lam3_n33**-2.0))*factor
	flux3_n34 = (flux_3_n34*(lam3_n34**-2.0))*factor
	flux3_n35 = (flux_3_n35*(lam3_n35**-2.0))*factor
	flux3_n36 = (flux_3_n36*(lam3_n36**-2.0))*factor
	flux3_n37 = (flux_3_n37*(lam3_n37**-2.0))*factor
	flux3_n38 = (flux_3_n38*(lam3_n38**-2.0))*factor
	
	
	
	
	os.chdir("/Volumes/TOSHIBA EXT/3d_hst/noah_massive/goodss_massive")
	
	data1_s1 = ascii.read("1523.temp_sed")
	data1_s2 = ascii.read("1924.temp_sed")
	data1_s3 = ascii.read("2707.temp_sed")
	data1_s4 = ascii.read("4210.temp_sed")
	data1_s5 = ascii.read("6098.temp_sed")
	data1_s6 = ascii.read("6106.temp_sed")
	data1_s7 = ascii.read("7444.temp_sed")
	data1_s8 = ascii.read("7503.temp_sed")
	data1_s9 = ascii.read("19186.temp_sed")
	data1_s10 = ascii.read("27442.temp_sed")
	data1_s11 = ascii.read("29928.temp_sed")
	data1_s12 = ascii.read("30394.temp_sed")
	data1_s13 = ascii.read("30997.temp_sed")
	data1_s14 = ascii.read("33163.temp_sed")
	data1_s15 = ascii.read("33164.temp_sed")
	data1_s16 = ascii.read("38111.temp_sed")
	data1_s17 = ascii.read("39170.temp_sed")
	data1_s18 = ascii.read("43042.temp_sed")
	data1_s19 = ascii.read("45775.temp_sed")
	data1_s20 = ascii.read("46392.temp_sed")
	data1_s21 = ascii.read("46846.temp_sed")
	data1_s22 = ascii.read("47742.temp_sed")
	data1_s23 = ascii.read("47873.temp_sed")
	data1_s24 = ascii.read("48631.temp_sed")
	
	
	lam1_s1 = data1_s1["lambda"]
	lam1_s2 = data1_s2["lambda"]
	lam1_s3 = data1_s3["lambda"]
	lam1_s4 = data1_s4["lambda"]
	lam1_s5 = data1_s5["lambda"]
	lam1_s6 = data1_s6["lambda"]
	lam1_s7 = data1_s7["lambda"]
	lam1_s8 = data1_s8["lambda"]
	lam1_s9 = data1_s9["lambda"]
	lam1_s10 = data1_s10["lambda"]
	lam1_s11 = data1_s11["lambda"]
	lam1_s12 = data1_s12["lambda"]
	lam1_s13 = data1_s13["lambda"]
	lam1_s14 = data1_s14["lambda"]
	lam1_s15 = data1_s15["lambda"]
	lam1_s16 = data1_s16["lambda"]
	lam1_s17 = data1_s17["lambda"]
	lam1_s18 = data1_s18["lambda"]
	lam1_s19 = data1_s19["lambda"]
	lam1_s20 = data1_s20["lambda"]
	lam1_s21 = data1_s21["lambda"]
	lam1_s22 = data1_s22["lambda"]
	lam1_s23 = data1_s23["lambda"]
	lam1_s24 = data1_s24["lambda"]
	

	flux_1_s1 = data1_s1["tempflux"]
	flux_1_s2 = data1_s2["tempflux"]
	flux_1_s3 = data1_s3["tempflux"]
	flux_1_s4 = data1_s4["tempflux"]
	flux_1_s5 = data1_s5["tempflux"]
	flux_1_s6 = data1_s6["tempflux"]
	flux_1_s7 = data1_s7["tempflux"]
	flux_1_s8 = data1_s8["tempflux"]
	flux_1_s9 = data1_s9["tempflux"]
	flux_1_s10 = data1_s10["tempflux"]
	flux_1_s11 = data1_s11["tempflux"]
	flux_1_s12 = data1_s12["tempflux"]
	flux_1_s13 = data1_s13["tempflux"]
	flux_1_s14 = data1_s14["tempflux"]
	flux_1_s15 = data1_s15["tempflux"]
	flux_1_s16 = data1_s16["tempflux"]
	flux_1_s17 = data1_s17["tempflux"]
	flux_1_s18 = data1_s18["tempflux"]
	flux_1_s19 = data1_s19["tempflux"]
	flux_1_s20 = data1_s20["tempflux"]
	flux_1_s21 = data1_s21["tempflux"]
	flux_1_s22 = data1_s22["tempflux"]
	flux_1_s23 = data1_s23["tempflux"]
	flux_1_s24 = data1_s24["tempflux"]
	
	
	flux1_s1 = (flux_1_s1*(lam1_s1**-2.0))*factor
	flux1_s2 = (flux_1_s2*(lam1_s2**-2.0))*factor
	flux1_s3 = (flux_1_s3*(lam1_s3**-2.0))*factor
	flux1_s4 = (flux_1_s4*(lam1_s4**-2.0))*factor
	flux1_s5 = (flux_1_s5*(lam1_s5**-2.0))*factor
	flux1_s6 = (flux_1_s6*(lam1_s6**-2.0))*factor
	flux1_s7 = (flux_1_s7*(lam1_s7**-2.0))*factor
	flux1_s8 = (flux_1_s8*(lam1_s8**-2.0))*factor
	flux1_s9 = (flux_1_s9*(lam1_s9**-2.0))*factor
	flux1_s10 = (flux_1_s10*(lam1_s10**-2.0))*factor
	flux1_s11 = (flux_1_s11*(lam1_s11**-2.0))*factor
	flux1_s12 = (flux_1_s12*(lam1_s12**-2.0))*factor
	flux1_s13 = (flux_1_s13*(lam1_s13**-2.0))*factor
	flux1_s14 = (flux_1_s14*(lam1_s14**-2.0))*factor
	flux1_s15 = (flux_1_s15*(lam1_s15**-2.0))*factor
	flux1_s16 = (flux_1_s16*(lam1_s16**-2.0))*factor
	flux1_s17 = (flux_1_s17*(lam1_s17**-2.0))*factor
	flux1_s18 = (flux_1_s18*(lam1_s18**-2.0))*factor
	flux1_s19 = (flux_1_s19*(lam1_s19**-2.0))*factor
	flux1_s20 = (flux_1_s20*(lam1_s20**-2.0))*factor
	flux1_s21 = (flux_1_s21*(lam1_s21**-2.0))*factor
	flux1_s22 = (flux_1_s22*(lam1_s22**-2.0))*factor
	flux1_s23 = (flux_1_s23*(lam1_s23**-2.0))*factor
	flux1_s24 = (flux_1_s24*(lam1_s24**-2.0))*factor
	
	
	data2_s1 = ascii.read("2383.temp_sed")
	data2_s2 = ascii.read("4505.temp_sed")
	data2_s3 = ascii.read("7457.temp_sed")
	data2_s4 = ascii.read("8422.temp_sed")
	data2_s5 = ascii.read("9704.temp_sed")
	data2_s6 = ascii.read("10436.temp_sed")
	data2_s7 = ascii.read("13369.temp_sed")
	data2_s8 = ascii.read("13628.temp_sed")
	data2_s9 = ascii.read("14152.temp_sed")
	data2_s10 = ascii.read("14335.temp_sed")
	data2_s11 = ascii.read("14747.temp_sed")
	data2_s12 = ascii.read("15214.temp_sed")
	data2_s13 = ascii.read("16769.temp_sed")
	data2_s14 = ascii.read("16814.temp_sed")
	data2_s15 = ascii.read("24308.temp_sed")
	data2_s16 = ascii.read("26139.temp_sed")
	data2_s17 = ascii.read("27881.temp_sed")
	data2_s18 = ascii.read("29407.temp_sed")
	data2_s19 = ascii.read("29652.temp_sed")
	data2_s20 = ascii.read("29900.temp_sed")
	data2_s21 = ascii.read("31397.temp_sed")
	data2_s22 = ascii.read("32048.temp_sed")
	data2_s23 = ascii.read("32783.temp_sed")
	data2_s24 = ascii.read("33912.temp_sed")
	data2_s25 = ascii.read("34491.temp_sed")
	data2_s26 = ascii.read("34519.temp_sed")
	data2_s27 = ascii.read("34567.temp_sed")
	data2_s28 = ascii.read("35444.temp_sed")
	data2_s29 = ascii.read("36095.temp_sed")
	data2_s30 = ascii.read("39012.temp_sed")
	data2_s31 = ascii.read("39208.temp_sed")
	data2_s32 = ascii.read("39364.temp_sed")
	data2_s33 = ascii.read("40889.temp_sed")
	data2_s34 = ascii.read("41148.temp_sed")
	data2_s35 = ascii.read("42113.temp_sed")
	data2_s36 = ascii.read("42501.temp_sed")
	data2_s37 = ascii.read("42705.temp_sed")
	data2_s38 = ascii.read("42957.temp_sed")
	data2_s39 = ascii.read("43114.temp_sed")
	data2_s40 = ascii.read("44042.temp_sed")
	data2_s41 = ascii.read("44157.temp_sed")
	data2_s42 = ascii.read("48312.temp_sed")
	
	
	lam2_s1 = data2_s1["lambda"]
	lam2_s2 = data2_s2["lambda"]
	lam2_s3 = data2_s3["lambda"]
	lam2_s4 = data2_s4["lambda"]
	lam2_s5 = data2_s5["lambda"]
	lam2_s6 = data2_s6["lambda"]
	lam2_s7 = data2_s7["lambda"]
	lam2_s8 = data2_s8["lambda"]
	lam2_s9 = data2_s9["lambda"]
	lam2_s10 = data2_s10["lambda"]
	lam2_s11 = data2_s11["lambda"]
	lam2_s12 = data2_s12["lambda"]
	lam2_s13 = data2_s13["lambda"]
	lam2_s14 = data2_s14["lambda"]
	lam2_s15 = data2_s15["lambda"]
	lam2_s16 = data2_s16["lambda"]
	lam2_s17 = data2_s17["lambda"]
	lam2_s18 = data2_s18["lambda"]
	lam2_s19 = data2_s19["lambda"]
	lam2_s20 = data2_s20["lambda"]
	lam2_s21 = data2_s21["lambda"]
	lam2_s22 = data2_s22["lambda"]
	lam2_s23 = data2_s23["lambda"]
	lam2_s24 = data2_s24["lambda"]
	lam2_s25 = data2_s25["lambda"]
	lam2_s26 = data2_s26["lambda"]
	lam2_s27 = data2_s27["lambda"]
	lam2_s28 = data2_s28["lambda"]
	lam2_s29 = data2_s29["lambda"]
	lam2_s30 = data2_s30["lambda"]
	lam2_s31 = data2_s31["lambda"]
	lam2_s32 = data2_s32["lambda"]
	lam2_s33 = data2_s33["lambda"]
	lam2_s34 = data2_s34["lambda"]
	lam2_s35 = data2_s35["lambda"]
	lam2_s36 = data2_s36["lambda"]
	lam2_s37 = data2_s37["lambda"]
	lam2_s38 = data2_s38["lambda"]
	lam2_s39 = data2_s39["lambda"]
	lam2_s40 = data2_s40["lambda"]
	lam2_s41 = data2_s41["lambda"]
	lam2_s42 = data2_s42["lambda"]
	
	
	flux_2_s1 = data2_s1["tempflux"]
	flux_2_s2 = data2_s2["tempflux"]
	flux_2_s3 = data2_s3["tempflux"]
	flux_2_s4 = data2_s4["tempflux"]
	flux_2_s5 = data2_s5["tempflux"]
	flux_2_s6 = data2_s6["tempflux"]
	flux_2_s7 = data2_s7["tempflux"]
	flux_2_s8 = data2_s8["tempflux"]
	flux_2_s9 = data2_s9["tempflux"]
	flux_2_s10 = data2_s10["tempflux"]
	flux_2_s11 = data2_s11["tempflux"]
	flux_2_s12 = data2_s12["tempflux"]
	flux_2_s13 = data2_s13["tempflux"]
	flux_2_s14 = data2_s14["tempflux"]
	flux_2_s15 = data2_s15["tempflux"]
	flux_2_s16 = data2_s16["tempflux"]
	flux_2_s17 = data2_s17["tempflux"]
	flux_2_s18 = data2_s18["tempflux"]
	flux_2_s19 = data2_s19["tempflux"]
	flux_2_s20 = data2_s20["tempflux"]
	flux_2_s21 = data2_s21["tempflux"]
	flux_2_s22 = data2_s22["tempflux"]
	flux_2_s23 = data2_s23["tempflux"]
	flux_2_s24 = data2_s24["tempflux"]
	flux_2_s25 = data2_s25["tempflux"]
	flux_2_s26 = data2_s26["tempflux"]
	flux_2_s27 = data2_s27["tempflux"]
	flux_2_s28 = data2_s28["tempflux"]
	flux_2_s29 = data2_s29["tempflux"]
	flux_2_s30 = data2_s30["tempflux"]
	flux_2_s31 = data2_s31["tempflux"]
	flux_2_s32 = data2_s32["tempflux"]
	flux_2_s33 = data2_s33["tempflux"]
	flux_2_s34 = data2_s34["tempflux"]
	flux_2_s35 = data2_s35["tempflux"]
	flux_2_s36 = data2_s36["tempflux"]
	flux_2_s37 = data2_s37["tempflux"]
	flux_2_s38 = data2_s38["tempflux"]
	flux_2_s39 = data2_s39["tempflux"]
	flux_2_s40 = data2_s40["tempflux"]
	flux_2_s41 = data2_s41["tempflux"]
	flux_2_s42 = data2_s42["tempflux"]
	
	
	flux2_s1 = (flux_2_s1*(lam2_s1**-2.0))*factor
	flux2_s2 = (flux_2_s2*(lam2_s2**-2.0))*factor
	flux2_s3 = (flux_2_s3*(lam2_s3**-2.0))*factor
	flux2_s4 = (flux_2_s4*(lam2_s4**-2.0))*factor
	flux2_s5 = (flux_2_s5*(lam2_s5**-2.0))*factor
	flux2_s6 = (flux_2_s6*(lam2_s6**-2.0))*factor
	flux2_s7 = (flux_2_s7*(lam2_s7**-2.0))*factor
	flux2_s8 = (flux_2_s8*(lam2_s8**-2.0))*factor
	flux2_s9 = (flux_2_s9*(lam2_s9**-2.0))*factor
	flux2_s10 = (flux_2_s10*(lam2_s10**-2.0))*factor
	flux2_s11 = (flux_2_s11*(lam2_s11**-2.0))*factor
	flux2_s12 = (flux_2_s12*(lam2_s12**-2.0))*factor
	flux2_s13 = (flux_2_s13*(lam2_s13**-2.0))*factor
	flux2_s14 = (flux_2_s14*(lam2_s14**-2.0))*factor
	flux2_s15 = (flux_2_s15*(lam2_s15**-2.0))*factor
	flux2_s16 = (flux_2_s16*(lam2_s16**-2.0))*factor
	flux2_s17 = (flux_2_s17*(lam2_s17**-2.0))*factor
	flux2_s18 = (flux_2_s18*(lam2_s18**-2.0))*factor
	flux2_s19 = (flux_2_s19*(lam2_s19**-2.0))*factor
	flux2_s20 = (flux_2_s20*(lam2_s20**-2.0))*factor
	flux2_s21 = (flux_2_s21*(lam2_s21**-2.0))*factor
	flux2_s22 = (flux_2_s22*(lam2_s22**-2.0))*factor
	flux2_s23 = (flux_2_s23*(lam2_s23**-2.0))*factor
	flux2_s24 = (flux_2_s24*(lam2_s24**-2.0))*factor
	flux2_s25 = (flux_2_s25*(lam2_s25**-2.0))*factor
	flux2_s26 = (flux_2_s26*(lam2_s26**-2.0))*factor
	flux2_s27 = (flux_2_s27*(lam2_s27**-2.0))*factor
	flux2_s28 = (flux_2_s28*(lam2_s28**-2.0))*factor
	flux2_s29 = (flux_2_s29*(lam2_s29**-2.0))*factor
	flux2_s30 = (flux_2_s30*(lam2_s30**-2.0))*factor
	flux2_s31 = (flux_2_s31*(lam2_s31**-2.0))*factor
	flux2_s32 = (flux_2_s32*(lam2_s32**-2.0))*factor
	flux2_s33 = (flux_2_s33*(lam2_s33**-2.0))*factor
	flux2_s34 = (flux_2_s34*(lam2_s34**-2.0))*factor
	flux2_s35 = (flux_2_s35*(lam2_s35**-2.0))*factor
	flux2_s36 = (flux_2_s36*(lam2_s36**-2.0))*factor
	flux2_s37 = (flux_2_s37*(lam2_s37**-2.0))*factor
	flux2_s38 = (flux_2_s38*(lam2_s38**-2.0))*factor
	flux2_s39 = (flux_2_s39*(lam2_s39**-2.0))*factor
	flux2_s40 = (flux_2_s40*(lam2_s40**-2.0))*factor
	flux2_s41 = (flux_2_s41*(lam2_s41**-2.0))*factor
	flux2_s42 = (flux_2_s42*(lam2_s42**-2.0))*factor
	
	
	data3_s1 = ascii.read("1725.temp_sed")
	data3_s2 = ascii.read("2467.temp_sed")
	data3_s3 = ascii.read("4583.temp_sed")
	data3_s4 = ascii.read("6341.temp_sed")
	data3_s5 = ascii.read("7686.temp_sed")
	data3_s6 = ascii.read("8683.temp_sed")
	data3_s7 = ascii.read("11016.temp_sed")
	data3_s8 = ascii.read("14813.temp_sed")
	data3_s9 = ascii.read("15847.temp_sed")
	data3_s10 = ascii.read("16888.temp_sed")
	data3_s11 = ascii.read("22079.temp_sed")
	data3_s12 = ascii.read("22825.temp_sed")
	data3_s13 = ascii.read("28604.temp_sed")
	data3_s14 = ascii.read("29288.temp_sed")
	data3_s15 = ascii.read("30274.temp_sed")
	data3_s16 = ascii.read("30534.temp_sed")
	data3_s17 = ascii.read("39568.temp_sed")
	data3_s18 = ascii.read("40185.temp_sed")
	data3_s19 = ascii.read("41181.temp_sed")
	data3_s20 = ascii.read("42607.temp_sed")
	data3_s21 = ascii.read("43901.temp_sed")
	data3_s22 = ascii.read("45475.temp_sed")
	data3_s23 = ascii.read("48464.temp_sed")
	data3_s24 = ascii.read("49285.temp_sed")
	data3_s25 = ascii.read("49834.temp_sed")
	
	
	lam3_s1 = data3_s1["lambda"]
	lam3_s2 = data3_s2["lambda"]
	lam3_s3 = data3_s3["lambda"]
	lam3_s4 = data3_s4["lambda"]
	lam3_s5 = data3_s5["lambda"]
	lam3_s6 = data3_s6["lambda"]
	lam3_s7 = data3_s7["lambda"]
	lam3_s8 = data3_s8["lambda"]
	lam3_s9 = data3_s9["lambda"]
	lam3_s10 = data3_s10["lambda"]
	lam3_s11 = data3_s11["lambda"]
	lam3_s12 = data3_s12["lambda"]
	lam3_s13 = data3_s13["lambda"]
	lam3_s14 = data3_s14["lambda"]
	lam3_s15 = data3_s15["lambda"]
	lam3_s16 = data3_s16["lambda"]
	lam3_s17 = data3_s17["lambda"]
	lam3_s18 = data3_s18["lambda"]
	lam3_s19 = data3_s19["lambda"]
	lam3_s20 = data3_s20["lambda"]
	lam3_s21 = data3_s21["lambda"]
	lam3_s22 = data3_s22["lambda"]
	lam3_s23 = data3_s23["lambda"]
	lam3_s24 = data3_s24["lambda"]
	lam3_s25 = data3_s25["lambda"]
	
	
	flux_3_s1 = data3_s1["tempflux"]
	flux_3_s2 = data3_s2["tempflux"]
	flux_3_s3 = data3_s3["tempflux"]
	flux_3_s4 = data3_s4["tempflux"]
	flux_3_s5 = data3_s5["tempflux"]
	flux_3_s6 = data3_s6["tempflux"]
	flux_3_s7 = data3_s7["tempflux"]
	flux_3_s8 = data3_s8["tempflux"]
	flux_3_s9 = data3_s9["tempflux"]
	flux_3_s10 = data3_s10["tempflux"]
	flux_3_s11 = data3_s11["tempflux"]
	flux_3_s12 = data3_s12["tempflux"]
	flux_3_s13 = data3_s13["tempflux"]
	flux_3_s14 = data3_s14["tempflux"]
	flux_3_s15 = data3_s15["tempflux"]
	flux_3_s16 = data3_s16["tempflux"]
	flux_3_s17 = data3_s17["tempflux"]
	flux_3_s18 = data3_s18["tempflux"]
	flux_3_s19 = data3_s19["tempflux"]
	flux_3_s20 = data3_s20["tempflux"]
	flux_3_s21 = data3_s21["tempflux"]
	flux_3_s22 = data3_s22["tempflux"]
	flux_3_s23 = data3_s23["tempflux"]
	flux_3_s24 = data3_s24["tempflux"]
	flux_3_s25 = data3_s25["tempflux"]
	

	flux3_s1 = (flux_3_s1*(lam3_s1**-2.0))*factor
	flux3_s2 = (flux_3_s2*(lam3_s2**-2.0))*factor
	flux3_s3 = (flux_3_s3*(lam3_s3**-2.0))*factor
	flux3_s4 = (flux_3_s4*(lam3_s4**-2.0))*factor
	flux3_s5 = (flux_3_s5*(lam3_s5**-2.0))*factor
	flux3_s6 = (flux_3_s6*(lam3_s6**-2.0))*factor
	flux3_s7 = (flux_3_s7*(lam3_s7**-2.0))*factor
	flux3_s8 = (flux_3_s8*(lam3_s8**-2.0))*factor
	flux3_s9 = (flux_3_s9*(lam3_s9**-2.0))*factor
	flux3_s10 = (flux_3_s10*(lam3_s10**-2.0))*factor
	flux3_s11 = (flux_3_s11*(lam3_s11**-2.0))*factor
	flux3_s12 = (flux_3_s12*(lam3_s12**-2.0))*factor
	flux3_s13 = (flux_3_s13*(lam3_s13**-2.0))*factor
	flux3_s14 = (flux_3_s14*(lam3_s14**-2.0))*factor
	flux3_s15 = (flux_3_s15*(lam3_s15**-2.0))*factor
	flux3_s16 = (flux_3_s16*(lam3_s16**-2.0))*factor
	flux3_s17 = (flux_3_s17*(lam3_s17**-2.0))*factor
	flux3_s18 = (flux_3_s18*(lam3_s18**-2.0))*factor
	flux3_s19 = (flux_3_s19*(lam3_s19**-2.0))*factor
	flux3_s20 = (flux_3_s20*(lam3_s20**-2.0))*factor
	flux3_s21 = (flux_3_s21*(lam3_s21**-2.0))*factor
	flux3_s22 = (flux_3_s22*(lam3_s22**-2.0))*factor
	flux3_s23 = (flux_3_s23*(lam3_s23**-2.0))*factor
	flux3_s24 = (flux_3_s24*(lam3_s24**-2.0))*factor
	flux3_s25 = (flux_3_s25*(lam3_s25**-2.0))*factor
	
	
	
	
	os.chdir("/Volumes/TOSHIBA EXT/3d_hst/noah_massive/uds_massive")
	
	data1_u1 = ascii.read("1123.temp_sed")
	data1_u2 = ascii.read("2393.temp_sed")
	data1_u3 = ascii.read("2394.temp_sed")
	data1_u4 = ascii.read("5126.temp_sed")
	data1_u5 = ascii.read("5924.temp_sed")
	data1_u6 = ascii.read("6299.temp_sed")
	data1_u7 = ascii.read("6852.temp_sed")
	data1_u8 = ascii.read("7071.temp_sed")
	data1_u9 = ascii.read("7783.temp_sed")
	data1_u10 = ascii.read("9261.temp_sed")
	data1_u11 = ascii.read("10758.temp_sed")
	data1_u12 = ascii.read("11533.temp_sed")
	data1_u13 = ascii.read("13482.temp_sed")
	data1_u14 = ascii.read("14152.temp_sed")
	data1_u15 = ascii.read("14854.temp_sed")
	data1_u16 = ascii.read("15063.temp_sed")
	data1_u17 = ascii.read("16239.temp_sed")
	data1_u18 = ascii.read("17879.temp_sed")
	data1_u19 = ascii.read("19765.temp_sed")
	data1_u20 = ascii.read("19954.temp_sed")
	data1_u21 = ascii.read("21513.temp_sed")
	data1_u22 = ascii.read("23590.temp_sed")
	data1_u23 = ascii.read("24953.temp_sed")
	data1_u24 = ascii.read("26552.temp_sed")
	data1_u25 = ascii.read("26875.temp_sed")
	data1_u26 = ascii.read("28773.temp_sed")
	data1_u27 = ascii.read("30057.temp_sed")
	data1_u28 = ascii.read("30192.temp_sed")
	data1_u29 = ascii.read("30255.temp_sed")
	data1_u30 = ascii.read("32077.temp_sed")
	data1_u31 = ascii.read("32256.temp_sed")
	data1_u32 = ascii.read("32691.temp_sed")
	data1_u33 = ascii.read("32777.temp_sed")
	data1_u34 = ascii.read("32921.temp_sed")
	data1_u35 = ascii.read("32986.temp_sed")
	data1_u36 = ascii.read("34353.temp_sed")
	data1_u37 = ascii.read("34641.temp_sed")
	data1_u38 = ascii.read("35071.temp_sed")
	data1_u39 = ascii.read("35356.temp_sed")
	data1_u40 = ascii.read("36013.temp_sed")
	data1_u41 = ascii.read("40631.temp_sed")
	data1_u42 = ascii.read("41412.temp_sed")
	data1_u43 = ascii.read("41671.temp_sed")
	data1_u44 = ascii.read("41835.temp_sed")
	
	
	lam1_u1 = data1_u1["lambda"]
	lam1_u2 = data1_u2["lambda"]
	lam1_u3 = data1_u3["lambda"]
	lam1_u4 = data1_u4["lambda"]
	lam1_u5 = data1_u5["lambda"]
	lam1_u6 = data1_u6["lambda"]
	lam1_u7 = data1_u7["lambda"]
	lam1_u8 = data1_u8["lambda"]
	lam1_u9 = data1_u9["lambda"]
	lam1_u10 = data1_u10["lambda"]
	lam1_u11 = data1_u11["lambda"]
	lam1_u12 = data1_u12["lambda"]
	lam1_u13 = data1_u13["lambda"]
	lam1_u14 = data1_u14["lambda"]
	lam1_u15 = data1_u15["lambda"]
	lam1_u16 = data1_u16["lambda"]
	lam1_u17 = data1_u17["lambda"]
	lam1_u18 = data1_u18["lambda"]
	lam1_u19 = data1_u19["lambda"]
	lam1_u20 = data1_u20["lambda"]
	lam1_u21 = data1_u21["lambda"]
	lam1_u22 = data1_u22["lambda"]
	lam1_u23 = data1_u23["lambda"]
	lam1_u24 = data1_u24["lambda"]
	lam1_u25 = data1_u25["lambda"]
	lam1_u26 = data1_u26["lambda"]
	lam1_u27 = data1_u27["lambda"]
	lam1_u28 = data1_u28["lambda"]
	lam1_u29 = data1_u29["lambda"]
	lam1_u30 = data1_u30["lambda"]
	lam1_u31 = data1_u31["lambda"]
	lam1_u32 = data1_u32["lambda"]
	lam1_u33 = data1_u33["lambda"]
	lam1_u34 = data1_u34["lambda"]
	lam1_u35 = data1_u35["lambda"]
	lam1_u36 = data1_u36["lambda"]
	lam1_u37 = data1_u37["lambda"]
	lam1_u38 = data1_u38["lambda"]
	lam1_u39 = data1_u39["lambda"]
	lam1_u40 = data1_u40["lambda"]
	lam1_u41 = data1_u41["lambda"]
	lam1_u42 = data1_u42["lambda"]
	lam1_u43 = data1_u43["lambda"]
	lam1_u44 = data1_u44["lambda"]
	

	flux_1_u1 = data1_u1["tempflux"]
	flux_1_u2 = data1_u2["tempflux"]
	flux_1_u3 = data1_u3["tempflux"]
	flux_1_u4 = data1_u4["tempflux"]
	flux_1_u5 = data1_u5["tempflux"]
	flux_1_u6 = data1_u6["tempflux"]
	flux_1_u7 = data1_u7["tempflux"]
	flux_1_u8 = data1_u8["tempflux"]
	flux_1_u9 = data1_u9["tempflux"]
	flux_1_u10 = data1_u10["tempflux"]
	flux_1_u11 = data1_u11["tempflux"]
	flux_1_u12 = data1_u12["tempflux"]
	flux_1_u13 = data1_u13["tempflux"]
	flux_1_u14 = data1_u14["tempflux"]
	flux_1_u15 = data1_u15["tempflux"]
	flux_1_u16 = data1_u16["tempflux"]
	flux_1_u17 = data1_u17["tempflux"]
	flux_1_u18 = data1_u18["tempflux"]
	flux_1_u19 = data1_u19["tempflux"]
	flux_1_u20 = data1_u20["tempflux"]
	flux_1_u21 = data1_u21["tempflux"]
	flux_1_u22 = data1_u22["tempflux"]
	flux_1_u23 = data1_u23["tempflux"]
	flux_1_u24 = data1_u24["tempflux"]
	flux_1_u25 = data1_u25["tempflux"]
	flux_1_u26 = data1_u26["tempflux"]
	flux_1_u27 = data1_u27["tempflux"]
	flux_1_u28 = data1_u28["tempflux"]
	flux_1_u29 = data1_u29["tempflux"]
	flux_1_u30 = data1_u30["tempflux"]
	flux_1_u31 = data1_u31["tempflux"]
	flux_1_u32 = data1_u32["tempflux"]
	flux_1_u33 = data1_u33["tempflux"]
	flux_1_u34 = data1_u34["tempflux"]
	flux_1_u35 = data1_u35["tempflux"]
	flux_1_u36 = data1_u36["tempflux"]
	flux_1_u37 = data1_u37["tempflux"]
	flux_1_u38 = data1_u38["tempflux"]
	flux_1_u39 = data1_u39["tempflux"]
	flux_1_u40 = data1_u40["tempflux"]
	flux_1_u41 = data1_u41["tempflux"]
	flux_1_u42 = data1_u42["tempflux"]
	flux_1_u43 = data1_u43["tempflux"]
	flux_1_u44 = data1_u44["tempflux"]
	
	
	flux1_u1 = (flux_1_u1*(lam1_u1**-2.0))*factor
	flux1_u2 = (flux_1_u2*(lam1_u2**-2.0))*factor
	flux1_u3 = (flux_1_u3*(lam1_u3**-2.0))*factor
	flux1_u4 = (flux_1_u4*(lam1_u4**-2.0))*factor
	flux1_u5 = (flux_1_u5*(lam1_u5**-2.0))*factor
	flux1_u6 = (flux_1_u6*(lam1_u6**-2.0))*factor
	flux1_u7 = (flux_1_u7*(lam1_u7**-2.0))*factor
	flux1_u8 = (flux_1_u8*(lam1_u8**-2.0))*factor
	flux1_u9 = (flux_1_u9*(lam1_u9**-2.0))*factor
	flux1_u10 = (flux_1_u10*(lam1_u10**-2.0))*factor
	flux1_u11 = (flux_1_u11*(lam1_u11**-2.0))*factor
	flux1_u12 = (flux_1_u12*(lam1_u12**-2.0))*factor
	flux1_u13 = (flux_1_u13*(lam1_u13**-2.0))*factor
	flux1_u14 = (flux_1_u14*(lam1_u14**-2.0))*factor
	flux1_u15 = (flux_1_u15*(lam1_u15**-2.0))*factor
	flux1_u16 = (flux_1_u16*(lam1_u16**-2.0))*factor
	flux1_u17 = (flux_1_u17*(lam1_u17**-2.0))*factor
	flux1_u18 = (flux_1_u18*(lam1_u18**-2.0))*factor
	flux1_u19 = (flux_1_u19*(lam1_u19**-2.0))*factor
	flux1_u20 = (flux_1_u20*(lam1_u20**-2.0))*factor
	flux1_u21 = (flux_1_u21*(lam1_u21**-2.0))*factor
	flux1_u22 = (flux_1_u22*(lam1_u22**-2.0))*factor
	flux1_u23 = (flux_1_u23*(lam1_u23**-2.0))*factor
	flux1_u24 = (flux_1_u24*(lam1_u24**-2.0))*factor
	flux1_u25 = (flux_1_u25*(lam1_u25**-2.0))*factor
	flux1_u26 = (flux_1_u26*(lam1_u26**-2.0))*factor
	flux1_u27 = (flux_1_u27*(lam1_u27**-2.0))*factor
	flux1_u28 = (flux_1_u28*(lam1_u28**-2.0))*factor
	flux1_u29 = (flux_1_u29*(lam1_u29**-2.0))*factor
	flux1_u30 = (flux_1_u30*(lam1_u30**-2.0))*factor
	flux1_u31 = (flux_1_u31*(lam1_u31**-2.0))*factor
	flux1_u32 = (flux_1_u32*(lam1_u32**-2.0))*factor
	flux1_u33 = (flux_1_u33*(lam1_u33**-2.0))*factor
	flux1_u34 = (flux_1_u34*(lam1_u34**-2.0))*factor
	flux1_u35 = (flux_1_u35*(lam1_u35**-2.0))*factor
	flux1_u36 = (flux_1_u36*(lam1_u36**-2.0))*factor
	flux1_u37 = (flux_1_u37*(lam1_u37**-2.0))*factor
	flux1_u38 = (flux_1_u38*(lam1_u38**-2.0))*factor
	flux1_u39 = (flux_1_u39*(lam1_u39**-2.0))*factor
	flux1_u40 = (flux_1_u40*(lam1_u40**-2.0))*factor
	flux1_u41 = (flux_1_u41*(lam1_u41**-2.0))*factor
	flux1_u42 = (flux_1_u42*(lam1_u42**-2.0))*factor
	flux1_u43 = (flux_1_u43*(lam1_u43**-2.0))*factor
	flux1_u44 = (flux_1_u44*(lam1_u44**-2.0))*factor
	
	
	data2_u1 = ascii.read("922.temp_sed")
	data2_u2 = ascii.read("1513.temp_sed")
	data2_u3 = ascii.read("1831.temp_sed")
	data2_u4 = ascii.read("1854.temp_sed")
	data2_u5 = ascii.read("2294.temp_sed")
	data2_u6 = ascii.read("3445.temp_sed")
	data2_u7 = ascii.read("4721.temp_sed")
	data2_u8 = ascii.read("6590.temp_sed")
	data2_u9 = ascii.read("6764.temp_sed")
	data2_u10 = ascii.read("7258.temp_sed")
	data2_u11 = ascii.read("9073.temp_sed")
	data2_u12 = ascii.read("10237.temp_sed")
	data2_u13 = ascii.read("10604.temp_sed")
	data2_u14 = ascii.read("12441.temp_sed")
	data2_u15 = ascii.read("12778.temp_sed")
	data2_u16 = ascii.read("14723.temp_sed")
	data2_u17 = ascii.read("15270.temp_sed")
	data2_u18 = ascii.read("18803.temp_sed")
	data2_u19 = ascii.read("19572.temp_sed")
	data2_u20 = ascii.read("19703.temp_sed")
	data2_u21 = ascii.read("19708.temp_sed")
	data2_u22 = ascii.read("19850.temp_sed")
	data2_u23 = ascii.read("20529.temp_sed")
	data2_u24 = ascii.read("20917.temp_sed")
	data2_u25 = ascii.read("20941.temp_sed")
	data2_u26 = ascii.read("21031.temp_sed")
	data2_u27 = ascii.read("21267.temp_sed")
	data2_u28 = ascii.read("21665.temp_sed")
	data2_u29 = ascii.read("22480.temp_sed")
	data2_u30 = ascii.read("25206.temp_sed")
	data2_u31 = ascii.read("25394.temp_sed")
	data2_u32 = ascii.read("25630.temp_sed")
	data2_u33 = ascii.read("27672.temp_sed")
	data2_u34 = ascii.read("28791.temp_sed")
	data2_u35 = ascii.read("30133.temp_sed")
	data2_u36 = ascii.read("30737.temp_sed")
	data2_u37 = ascii.read("31684.temp_sed")
	data2_u38 = ascii.read("32468.temp_sed")
	data2_u39 = ascii.read("32707.temp_sed")
	data2_u40 = ascii.read("33422.temp_sed")
	data2_u41 = ascii.read("33527.temp_sed")
	data2_u42 = ascii.read("35616.temp_sed")
	data2_u43 = ascii.read("35829.temp_sed")
	data2_u44 = ascii.read("36010.temp_sed")
	data2_u45 = ascii.read("36685.temp_sed")
	data2_u46 = ascii.read("37182.temp_sed")
	data2_u47 = ascii.read("37775.temp_sed")
	data2_u48 = ascii.read("38246.temp_sed")
	data2_u49 = ascii.read("38288.temp_sed")
	data2_u50 = ascii.read("38631.temp_sed")
	data2_u51 = ascii.read("39349.temp_sed")
	data2_u52 = ascii.read("39487.temp_sed")
	data2_u53 = ascii.read("40420.temp_sed")
	data2_u54 = ascii.read("40472.temp_sed")
	data2_u55 = ascii.read("40849.temp_sed")
	data2_u56 = ascii.read("41302.temp_sed")
	data2_u57 = ascii.read("41456.temp_sed")
	data2_u58 = ascii.read("42319.temp_sed")
	data2_u59 = ascii.read("43367.temp_sed")
	
	
	lam2_u1 = data2_u1["lambda"]
	lam2_u2 = data2_u2["lambda"]
	lam2_u3 = data2_u3["lambda"]
	lam2_u4 = data2_u4["lambda"]
	lam2_u5 = data2_u5["lambda"]
	lam2_u6 = data2_u6["lambda"]
	lam2_u7 = data2_u7["lambda"]
	lam2_u8 = data2_u8["lambda"]
	lam2_u9 = data2_u9["lambda"]
	lam2_u10 = data2_u10["lambda"]
	lam2_u11 = data2_u11["lambda"]
	lam2_u12 = data2_u12["lambda"]
	lam2_u13 = data2_u13["lambda"]
	lam2_u14 = data2_u14["lambda"]
	lam2_u15 = data2_u15["lambda"]
	lam2_u16 = data2_u16["lambda"]
	lam2_u17 = data2_u17["lambda"]
	lam2_u18 = data2_u18["lambda"]
	lam2_u19 = data2_u19["lambda"]
	lam2_u20 = data2_u20["lambda"]
	lam2_u21 = data2_u21["lambda"]
	lam2_u22 = data2_u22["lambda"]
	lam2_u23 = data2_u23["lambda"]
	lam2_u24 = data2_u24["lambda"]
	lam2_u25 = data2_u25["lambda"]
	lam2_u26 = data2_u26["lambda"]
	lam2_u27 = data2_u27["lambda"]
	lam2_u28 = data2_u28["lambda"]
	lam2_u29 = data2_u29["lambda"]
	lam2_u30 = data2_u30["lambda"]
	lam2_u31 = data2_u31["lambda"]
	lam2_u32 = data2_u32["lambda"]
	lam2_u33 = data2_u33["lambda"]
	lam2_u34 = data2_u34["lambda"]
	lam2_u35 = data2_u35["lambda"]
	lam2_u36 = data2_u36["lambda"]
	lam2_u37 = data2_u37["lambda"]
	lam2_u38 = data2_u38["lambda"]
	lam2_u39 = data2_u39["lambda"]
	lam2_u40 = data2_u40["lambda"]
	lam2_u41 = data2_u41["lambda"]
	lam2_u42 = data2_u42["lambda"]
	lam2_u43 = data2_u43["lambda"]
	lam2_u44 = data2_u44["lambda"]
	lam2_u45 = data2_u45["lambda"]
	lam2_u46 = data2_u46["lambda"]
	lam2_u47 = data2_u47["lambda"]
	lam2_u48 = data2_u48["lambda"]
	lam2_u49 = data2_u49["lambda"]
	lam2_u50 = data2_u50["lambda"]
	lam2_u51 = data2_u51["lambda"]
	lam2_u52 = data2_u52["lambda"]
	lam2_u53 = data2_u53["lambda"]
	lam2_u54 = data2_u54["lambda"]
	lam2_u55 = data2_u55["lambda"]
	lam2_u56 = data2_u56["lambda"]
	lam2_u57 = data2_u57["lambda"]
	lam2_u58 = data2_u58["lambda"]
	lam2_u59 = data2_u59["lambda"]
	
	
	flux_2_u1 = data2_u1["tempflux"]
	flux_2_u2 = data2_u2["tempflux"]
	flux_2_u3 = data2_u3["tempflux"]
	flux_2_u4 = data2_u4["tempflux"]
	flux_2_u5 = data2_u5["tempflux"]
	flux_2_u6 = data2_u6["tempflux"]
	flux_2_u7 = data2_u7["tempflux"]
	flux_2_u8 = data2_u8["tempflux"]
	flux_2_u9 = data2_u9["tempflux"]
	flux_2_u10 = data2_u10["tempflux"]
	flux_2_u11 = data2_u11["tempflux"]
	flux_2_u12 = data2_u12["tempflux"]
	flux_2_u13 = data2_u13["tempflux"]
	flux_2_u14 = data2_u14["tempflux"]
	flux_2_u15 = data2_u15["tempflux"]
	flux_2_u16 = data2_u16["tempflux"]
	flux_2_u17 = data2_u17["tempflux"]
	flux_2_u18 = data2_u18["tempflux"]
	flux_2_u19 = data2_u19["tempflux"]
	flux_2_u20 = data2_u20["tempflux"]
	flux_2_u21 = data2_u21["tempflux"]
	flux_2_u22 = data2_u22["tempflux"]
	flux_2_u23 = data2_u23["tempflux"]
	flux_2_u24 = data2_u24["tempflux"]
	flux_2_u25 = data2_u25["tempflux"]
	flux_2_u26 = data2_u26["tempflux"]
	flux_2_u27 = data2_u27["tempflux"]
	flux_2_u28 = data2_u28["tempflux"]
	flux_2_u29 = data2_u29["tempflux"]
	flux_2_u30 = data2_u30["tempflux"]
	flux_2_u31 = data2_u31["tempflux"]
	flux_2_u32 = data2_u32["tempflux"]
	flux_2_u33 = data2_u33["tempflux"]
	flux_2_u34 = data2_u34["tempflux"]
	flux_2_u35 = data2_u35["tempflux"]
	flux_2_u36 = data2_u36["tempflux"]
	flux_2_u37 = data2_u37["tempflux"]
	flux_2_u38 = data2_u38["tempflux"]
	flux_2_u39 = data2_u39["tempflux"]
	flux_2_u40 = data2_u40["tempflux"]
	flux_2_u41 = data2_u41["tempflux"]
	flux_2_u42 = data2_u42["tempflux"]
	flux_2_u43 = data2_u43["tempflux"]
	flux_2_u44 = data2_u44["tempflux"]
	flux_2_u45 = data2_u45["tempflux"]
	flux_2_u46 = data2_u46["tempflux"]
	flux_2_u47 = data2_u47["tempflux"]
	flux_2_u48 = data2_u48["tempflux"]
	flux_2_u49 = data2_u49["tempflux"]
	flux_2_u50 = data2_u50["tempflux"]
	flux_2_u51 = data2_u51["tempflux"]
	flux_2_u52 = data2_u52["tempflux"]
	flux_2_u53 = data2_u53["tempflux"]
	flux_2_u54 = data2_u54["tempflux"]
	flux_2_u55 = data2_u55["tempflux"]
	flux_2_u56 = data2_u56["tempflux"]
	flux_2_u57 = data2_u57["tempflux"]
	flux_2_u58 = data2_u58["tempflux"]
	flux_2_u59 = data2_u59["tempflux"]

	
	flux2_u1 = (flux_2_u1*(lam2_u1**-2.0))*factor
	flux2_u2 = (flux_2_u2*(lam2_u2**-2.0))*factor
	flux2_u3 = (flux_2_u3*(lam2_u3**-2.0))*factor
	flux2_u4 = (flux_2_u4*(lam2_u4**-2.0))*factor
	flux2_u5 = (flux_2_u5*(lam2_u5**-2.0))*factor
	flux2_u6 = (flux_2_u6*(lam2_u6**-2.0))*factor
	flux2_u7 = (flux_2_u7*(lam2_u7**-2.0))*factor
	flux2_u8 = (flux_2_u8*(lam2_u8**-2.0))*factor
	flux2_u9 = (flux_2_u9*(lam2_u9**-2.0))*factor
	flux2_u10 = (flux_2_u10*(lam2_u10**-2.0))*factor
	flux2_u11 = (flux_2_u11*(lam2_u11**-2.0))*factor
	flux2_u12 = (flux_2_u12*(lam2_u12**-2.0))*factor
	flux2_u13 = (flux_2_u13*(lam2_u13**-2.0))*factor
	flux2_u14 = (flux_2_u14*(lam2_u14**-2.0))*factor
	flux2_u15 = (flux_2_u15*(lam2_u15**-2.0))*factor
	flux2_u16 = (flux_2_u16*(lam2_u16**-2.0))*factor
	flux2_u17 = (flux_2_u17*(lam2_u17**-2.0))*factor
	flux2_u18 = (flux_2_u18*(lam2_u18**-2.0))*factor
	flux2_u19 = (flux_2_u19*(lam2_u19**-2.0))*factor
	flux2_u20 = (flux_2_u20*(lam2_u20**-2.0))*factor
	flux2_u21 = (flux_2_u21*(lam2_u21**-2.0))*factor
	flux2_u22 = (flux_2_u22*(lam2_u22**-2.0))*factor
	flux2_u23 = (flux_2_u23*(lam2_u23**-2.0))*factor
	flux2_u24 = (flux_2_u24*(lam2_u24**-2.0))*factor
	flux2_u25 = (flux_2_u25*(lam2_u25**-2.0))*factor
	flux2_u26 = (flux_2_u26*(lam2_u26**-2.0))*factor
	flux2_u27 = (flux_2_u27*(lam2_u27**-2.0))*factor
	flux2_u28 = (flux_2_u28*(lam2_u28**-2.0))*factor
	flux2_u29 = (flux_2_u29*(lam2_u29**-2.0))*factor
	flux2_u30 = (flux_2_u30*(lam2_u30**-2.0))*factor
	flux2_u31 = (flux_2_u31*(lam2_u31**-2.0))*factor
	flux2_u32 = (flux_2_u32*(lam2_u32**-2.0))*factor
	flux2_u33 = (flux_2_u33*(lam2_u33**-2.0))*factor
	flux2_u34 = (flux_2_u34*(lam2_u34**-2.0))*factor
	flux2_u35 = (flux_2_u35*(lam2_u35**-2.0))*factor
	flux2_u36 = (flux_2_u36*(lam2_u36**-2.0))*factor
	flux2_u37 = (flux_2_u37*(lam2_u37**-2.0))*factor
	flux2_u38 = (flux_2_u38*(lam2_u38**-2.0))*factor
	flux2_u39 = (flux_2_u39*(lam2_u39**-2.0))*factor
	flux2_u40 = (flux_2_u40*(lam2_u40**-2.0))*factor
	flux2_u41 = (flux_2_u41*(lam2_u41**-2.0))*factor
	flux2_u42 = (flux_2_u42*(lam2_u42**-2.0))*factor
	flux2_u43 = (flux_2_u43*(lam2_u43**-2.0))*factor
	flux2_u44 = (flux_2_u44*(lam2_u44**-2.0))*factor
	flux2_u45 = (flux_2_u45*(lam2_u45**-2.0))*factor
	flux2_u46 = (flux_2_u46*(lam2_u46**-2.0))*factor
	flux2_u47 = (flux_2_u47*(lam2_u47**-2.0))*factor
	flux2_u48 = (flux_2_u48*(lam2_u48**-2.0))*factor
	flux2_u49 = (flux_2_u49*(lam2_u49**-2.0))*factor
	flux2_u50 = (flux_2_u50*(lam2_u50**-2.0))*factor
	flux2_u51 = (flux_2_u51*(lam2_u51**-2.0))*factor
	flux2_u52 = (flux_2_u52*(lam2_u52**-2.0))*factor
	flux2_u53 = (flux_2_u53*(lam2_u53**-2.0))*factor
	flux2_u54 = (flux_2_u54*(lam2_u54**-2.0))*factor
	flux2_u55 = (flux_2_u55*(lam2_u55**-2.0))*factor
	flux2_u56 = (flux_2_u56*(lam2_u56**-2.0))*factor
	flux2_u57 = (flux_2_u57*(lam2_u57**-2.0))*factor
	flux2_u58 = (flux_2_u58*(lam2_u58**-2.0))*factor
	flux2_u59 = (flux_2_u59*(lam2_u59**-2.0))*factor

	
	data3_u1 = ascii.read("190.temp_sed")
	data3_u2 = ascii.read("394.temp_sed")
	data3_u3 = ascii.read("1620.temp_sed")
	data3_u4 = ascii.read("2166.temp_sed")
	data3_u5 = ascii.read("2211.temp_sed")
	data3_u6 = ascii.read("4059.temp_sed")
	data3_u7 = ascii.read("4128.temp_sed")
	data3_u8 = ascii.read("4701.temp_sed")
	data3_u9 = ascii.read("4706.temp_sed")
	data3_u10 = ascii.read("5155.temp_sed")
	data3_u11 = ascii.read("7516.temp_sed")
	data3_u12 = ascii.read("9207.temp_sed")
	data3_u13 = ascii.read("9367.temp_sed")
	data3_u14 = ascii.read("11558.temp_sed")
	data3_u15 = ascii.read("12010.temp_sed")
	data3_u16 = ascii.read("13108.temp_sed")
	data3_u17 = ascii.read("14409.temp_sed")
	data3_u18 = ascii.read("14996.temp_sed")
	data3_u19 = ascii.read("15598.temp_sed")
	data3_u20 = ascii.read("16022.temp_sed")
	data3_u21 = ascii.read("16709.temp_sed")
	data3_u22 = ascii.read("17838.temp_sed")
	data3_u23 = ascii.read("19068.temp_sed")
	data3_u24 = ascii.read("19126.temp_sed")
	data3_u25 = ascii.read("20694.temp_sed")
	data3_u26 = ascii.read("20704.temp_sed")
	data3_u27 = ascii.read("20770.temp_sed")
	data3_u28 = ascii.read("21998.temp_sed")
	data3_u29 = ascii.read("22227.temp_sed")
	data3_u30 = ascii.read("22416.temp_sed")
	data3_u31 = ascii.read("22685.temp_sed")
	data3_u32 = ascii.read("23692.temp_sed")
	data3_u33 = ascii.read("26581.temp_sed")
	data3_u34 = ascii.read("28087.temp_sed")
	data3_u35 = ascii.read("29179.temp_sed")
	data3_u36 = ascii.read("29461.temp_sed")
	data3_u37 = ascii.read("30196.temp_sed")
	data3_u38 = ascii.read("30916.temp_sed")
	data3_u39 = ascii.read("31615.temp_sed")
	data3_u40 = ascii.read("32147.temp_sed")
	data3_u41 = ascii.read("32351.temp_sed")
	data3_u42 = ascii.read("32947.temp_sed")
	data3_u43 = ascii.read("34150.temp_sed")
	data3_u44 = ascii.read("34817.temp_sed")
	data3_u45 = ascii.read("36247.temp_sed")
	data3_u46 = ascii.read("38640.temp_sed")
	data3_u47 = ascii.read("39126.temp_sed")
	data3_u48 = ascii.read("39624.temp_sed")
	data3_u49 = ascii.read("42529.temp_sed")
	data3_u50 = ascii.read("42571.temp_sed")
	data3_u51 = ascii.read("42812.temp_sed")
	data3_u52 = ascii.read("43667.temp_sed")
	

	lam3_u1 = data3_u1["lambda"]
	lam3_u2 = data3_u2["lambda"]
	lam3_u3 = data3_u3["lambda"]
	lam3_u4 = data3_u4["lambda"]
	lam3_u5 = data3_u5["lambda"]
	lam3_u6 = data3_u6["lambda"]
	lam3_u7 = data3_u7["lambda"]
	lam3_u8 = data3_u8["lambda"]
	lam3_u9 = data3_u9["lambda"]
	lam3_u10 = data3_u10["lambda"]
	lam3_u11 = data3_u11["lambda"]
	lam3_u12 = data3_u12["lambda"]
	lam3_u13 = data3_u13["lambda"]
	lam3_u14 = data3_u14["lambda"]
	lam3_u15 = data3_u15["lambda"]
	lam3_u16 = data3_u16["lambda"]
	lam3_u17 = data3_u17["lambda"]
	lam3_u18 = data3_u18["lambda"]
	lam3_u19 = data3_u19["lambda"]
	lam3_u20 = data3_u20["lambda"]
	lam3_u21 = data3_u21["lambda"]
	lam3_u22 = data3_u22["lambda"]
	lam3_u23 = data3_u23["lambda"]
	lam3_u24 = data3_u24["lambda"]
	lam3_u25 = data3_u25["lambda"]
	lam3_u26 = data3_u26["lambda"]
	lam3_u27 = data3_u27["lambda"]
	lam3_u28 = data3_u28["lambda"]
	lam3_u29 = data3_u29["lambda"]
	lam3_u30 = data3_u30["lambda"]
	lam3_u31 = data3_u31["lambda"]
	lam3_u32 = data3_u32["lambda"]
	lam3_u33 = data3_u33["lambda"]
	lam3_u34 = data3_u34["lambda"]
	lam3_u35 = data3_u35["lambda"]
	lam3_u36 = data3_u36["lambda"]
	lam3_u37 = data3_u37["lambda"]
	lam3_u38 = data3_u38["lambda"]
	lam3_u39 = data3_u39["lambda"]
	lam3_u40 = data3_u40["lambda"]
	lam3_u41 = data3_u41["lambda"]
	lam3_u42 = data3_u42["lambda"]
	lam3_u43 = data3_u43["lambda"]
	lam3_u44 = data3_u44["lambda"]
	lam3_u45 = data3_u45["lambda"]
	lam3_u46 = data3_u46["lambda"]
	lam3_u47 = data3_u47["lambda"]
	lam3_u48 = data3_u48["lambda"]
	lam3_u49 = data3_u49["lambda"]
	lam3_u50 = data3_u50["lambda"]
	lam3_u51 = data3_u51["lambda"]
	lam3_u52 = data3_u52["lambda"]
	
	
	flux_3_u1 = data3_u1["tempflux"]
	flux_3_u2 = data3_u2["tempflux"]
	flux_3_u3 = data3_u3["tempflux"]
	flux_3_u4 = data3_u4["tempflux"]
	flux_3_u5 = data3_u5["tempflux"]
	flux_3_u6 = data3_u6["tempflux"]
	flux_3_u7 = data3_u7["tempflux"]
	flux_3_u8 = data3_u8["tempflux"]
	flux_3_u9 = data3_u9["tempflux"]
	flux_3_u10 = data3_u10["tempflux"]
	flux_3_u11 = data3_u11["tempflux"]
	flux_3_u12 = data3_u12["tempflux"]
	flux_3_u13 = data3_u13["tempflux"]
	flux_3_u14 = data3_u14["tempflux"]
	flux_3_u15 = data3_u15["tempflux"]
	flux_3_u16 = data3_u16["tempflux"]
	flux_3_u17 = data3_u17["tempflux"]
	flux_3_u18 = data3_u18["tempflux"]
	flux_3_u19 = data3_u19["tempflux"]
	flux_3_u20 = data3_u20["tempflux"]
	flux_3_u21 = data3_u21["tempflux"]
	flux_3_u22 = data3_u22["tempflux"]
	flux_3_u23 = data3_u23["tempflux"]
	flux_3_u24 = data3_u24["tempflux"]
	flux_3_u25 = data3_u25["tempflux"]
	flux_3_u26 = data3_u26["tempflux"]
	flux_3_u27 = data3_u27["tempflux"]
	flux_3_u28 = data3_u28["tempflux"]
	flux_3_u29 = data3_u29["tempflux"]
	flux_3_u30 = data3_u30["tempflux"]
	flux_3_u31 = data3_u31["tempflux"]
	flux_3_u32 = data3_u32["tempflux"]
	flux_3_u33 = data3_u33["tempflux"]
	flux_3_u34 = data3_u34["tempflux"]
	flux_3_u35 = data3_u35["tempflux"]
	flux_3_u36 = data3_u36["tempflux"]
	flux_3_u37 = data3_u37["tempflux"]
	flux_3_u38 = data3_u38["tempflux"]
	flux_3_u39 = data3_u39["tempflux"]
	flux_3_u40 = data3_u40["tempflux"]
	flux_3_u41 = data3_u41["tempflux"]
	flux_3_u42 = data3_u42["tempflux"]
	flux_3_u43 = data3_u43["tempflux"]
	flux_3_u44 = data3_u44["tempflux"]
	flux_3_u45 = data3_u45["tempflux"]
	flux_3_u46 = data3_u46["tempflux"]
	flux_3_u47 = data3_u47["tempflux"]
	flux_3_u48 = data3_u48["tempflux"]
	flux_3_u49 = data3_u49["tempflux"]
	flux_3_u50 = data3_u50["tempflux"]
	flux_3_u51 = data3_u51["tempflux"]
	flux_3_u52 = data3_u52["tempflux"]
	

	flux3_u1 = (flux_3_u1*(lam3_u1**-2.0))*factor
	flux3_u2 = (flux_3_u2*(lam3_u2**-2.0))*factor
	flux3_u3 = (flux_3_u3*(lam3_u3**-2.0))*factor
	flux3_u4 = (flux_3_u4*(lam3_u4**-2.0))*factor
	flux3_u5 = (flux_3_u5*(lam3_u5**-2.0))*factor
	flux3_u6 = (flux_3_u6*(lam3_u6**-2.0))*factor
	flux3_u7 = (flux_3_u7*(lam3_u7**-2.0))*factor
	flux3_u8 = (flux_3_u8*(lam3_u8**-2.0))*factor
	flux3_u9 = (flux_3_u9*(lam3_u9**-2.0))*factor
	flux3_u10 = (flux_3_u10*(lam3_u10**-2.0))*factor
	flux3_u11 = (flux_3_u11*(lam3_u11**-2.0))*factor
	flux3_u12 = (flux_3_u12*(lam3_u12**-2.0))*factor
	flux3_u13 = (flux_3_u13*(lam3_u13**-2.0))*factor
	flux3_u14 = (flux_3_u14*(lam3_u14**-2.0))*factor
	flux3_u15 = (flux_3_u15*(lam3_u15**-2.0))*factor
	flux3_u16 = (flux_3_u16*(lam3_u16**-2.0))*factor
	flux3_u17 = (flux_3_u17*(lam3_u17**-2.0))*factor
	flux3_u18 = (flux_3_u18*(lam3_u18**-2.0))*factor
	flux3_u19 = (flux_3_u19*(lam3_u19**-2.0))*factor
	flux3_u20 = (flux_3_u20*(lam3_u20**-2.0))*factor
	flux3_u21 = (flux_3_u21*(lam3_u21**-2.0))*factor
	flux3_u22 = (flux_3_u22*(lam3_u22**-2.0))*factor
	flux3_u23 = (flux_3_u23*(lam3_u23**-2.0))*factor
	flux3_u24 = (flux_3_u24*(lam3_u24**-2.0))*factor
	flux3_u25 = (flux_3_u25*(lam3_u25**-2.0))*factor
	flux3_u26 = (flux_3_u26*(lam3_u26**-2.0))*factor
	flux3_u27 = (flux_3_u27*(lam3_u27**-2.0))*factor
	flux3_u28 = (flux_3_u28*(lam3_u28**-2.0))*factor
	flux3_u29 = (flux_3_u29*(lam3_u29**-2.0))*factor
	flux3_u30 = (flux_3_u30*(lam3_u30**-2.0))*factor
	flux3_u31 = (flux_3_u31*(lam3_u31**-2.0))*factor
	flux3_u32 = (flux_3_u32*(lam3_u32**-2.0))*factor
	flux3_u33 = (flux_3_u33*(lam3_u33**-2.0))*factor
	flux3_u34 = (flux_3_u34*(lam3_u34**-2.0))*factor
	flux3_u35 = (flux_3_u35*(lam3_u35**-2.0))*factor
	flux3_u36 = (flux_3_u36*(lam3_u36**-2.0))*factor
	flux3_u37 = (flux_3_u37*(lam3_u37**-2.0))*factor
	flux3_u38 = (flux_3_u38*(lam3_u38**-2.0))*factor
	flux3_u39 = (flux_3_u39*(lam3_u39**-2.0))*factor
	flux3_u40 = (flux_3_u40*(lam3_u40**-2.0))*factor
	flux3_u41 = (flux_3_u41*(lam3_u41**-2.0))*factor
	flux3_u42 = (flux_3_u42*(lam3_u42**-2.0))*factor
	flux3_u43 = (flux_3_u43*(lam3_u43**-2.0))*factor
	flux3_u44 = (flux_3_u44*(lam3_u44**-2.0))*factor
	flux3_u45 = (flux_3_u45*(lam3_u45**-2.0))*factor
	flux3_u46 = (flux_3_u46*(lam3_u46**-2.0))*factor
	flux3_u47 = (flux_3_u47*(lam3_u47**-2.0))*factor
	flux3_u48 = (flux_3_u48*(lam3_u48**-2.0))*factor
	flux3_u49 = (flux_3_u49*(lam3_u49**-2.0))*factor
	flux3_u50 = (flux_3_u50*(lam3_u50**-2.0))*factor
	flux3_u51 = (flux_3_u51*(lam3_u51**-2.0))*factor
	flux3_u52 = (flux_3_u52*(lam3_u52**-2.0))*factor
	
	
	
	
	fig,pylab.axes = pylab.subplots(3, 1, sharex=True)

	a1 = pylab.axes[0]
	a2 = pylab.axes[1]
	a3 = pylab.axes[2]
	
	a1.plot(lam1_a1, flux1_a1, color="g", alpha=0.7)
	a1.plot(lam1_a2, flux1_a2, color="g", alpha=0.7)
	a1.plot(lam1_a3, flux1_a3, color="g", alpha=0.7)
	a1.plot(lam1_a4, flux1_a4, color="g", alpha=0.7)
	a1.plot(lam1_a5, flux1_a5, color="g", alpha=0.7)
	a1.plot(lam1_a6, flux1_a6, color="g", alpha=0.7)
	a1.plot(lam1_a7, flux1_a7, color="g", alpha=0.7)
	a1.plot(lam1_a8, flux1_a8, color="g", alpha=0.7)
	a1.plot(lam1_a9, flux1_a9, color="g", alpha=0.7)
	a1.plot(lam1_a10, flux1_a10, color="g", alpha=0.7)
	a1.plot(lam1_a11, flux1_a11, color="g", alpha=0.7)
	a1.plot(lam1_a12, flux1_a12, color="g", alpha=0.7)
	a1.plot(lam1_a13, flux1_a13, color="g", alpha=0.7)
	a1.plot(lam1_a14, flux1_a14, color="g", alpha=0.7)
	a1.plot(lam1_a15, flux1_a15, color="g", alpha=0.7)
	a1.plot(lam1_a16, flux1_a16, color="g", alpha=0.7)
	a1.plot(lam1_a17, flux1_a17, color="g", alpha=0.7)
	a1.plot(lam1_a18, flux1_a18, color="g", alpha=0.7)
	a1.plot(lam1_a19, flux1_a19, color="g", alpha=0.7)
	a1.plot(lam1_a20, flux1_a20, color="g", alpha=0.7)
	a1.plot(lam1_a21, flux1_a21, color="g", alpha=0.7)
	a1.plot(lam1_a22, flux1_a22, color="g", alpha=0.7)
	a1.plot(lam1_a23, flux1_a23, color="g", alpha=0.7)
	a1.plot(lam1_a24, flux1_a24, color="g", alpha=0.7)
	a1.plot(lam1_a25, flux1_a25, color="g", alpha=0.7)
	a1.plot(lam1_a26, flux1_a26, color="g", alpha=0.7)
	a1.plot(lam1_a27, flux1_a27, color="g", alpha=0.7)
	a1.plot(lam1_a28, flux1_a28, color="g", alpha=0.7)
	a1.plot(lam1_a29, flux1_a29, color="g", alpha=0.7)
	a1.plot(lam1_a30, flux1_a30, color="g", alpha=0.7)
	a1.plot(lam1_a31, flux1_a31, color="g", alpha=0.7)
	a1.plot(lam1_a32, flux1_a32, color="g", alpha=0.7)
	a1.plot(lam1_a33, flux1_a33, color="g", alpha=0.7)
	a1.plot(lam1_a34, flux1_a34, color="g", alpha=0.7)
	a1.plot(lam1_a35, flux1_a35, color="g", alpha=0.7)
	a1.plot(lam1_a36, flux1_a36, color="g", alpha=0.7)
	a1.plot(lam1_a37, flux1_a37, color="g", alpha=0.7)
	a1.plot(lam1_a38, flux1_a38, color="g", alpha=0.7)
	a1.plot(lam1_a39, flux1_a39, color="g", alpha=0.7)
	a1.plot(lam1_a40, flux1_a40, color="g", alpha=0.7)
	a1.plot(lam1_a41, flux1_a41, color="g", alpha=0.7)
	a1.plot(lam1_a42, flux1_a42, color="g", alpha=0.7)
	a1.plot(lam1_a43, flux1_a43, color="g", alpha=0.7)
	a1.plot(lam1_a44, flux1_a44, color="g", alpha=0.7)
	a1.plot(lam1_a45, flux1_a45, color="g", alpha=0.7)
	a1.plot(lam1_a46, flux1_a46, color="g", alpha=0.7)
	a1.plot(lam1_a47, flux1_a47, color="g", alpha=0.7)
	a1.plot(lam1_a48, flux1_a48, color="g", alpha=0.7)
	a1.plot(lam1_a49, flux1_a49, color="g", alpha=0.7)
	a1.plot(lam1_a50, flux1_a50, color="g", alpha=0.7)
	a1.plot(lam1_a51, flux1_a51, color="g", alpha=0.7)
	a1.plot(lam1_a52, flux1_a52, color="g", alpha=0.7)
	a1.plot(lam1_a53, flux1_a53, color="g", alpha=0.7)
	
	
	a2.plot(lam2_a1, flux2_a1, color="b", alpha=0.7)
	a2.plot(lam2_a2, flux2_a2, color="b", alpha=0.7)
	a2.plot(lam2_a3, flux2_a3, color="b", alpha=0.7)
	a2.plot(lam2_a4, flux2_a4, color="b", alpha=0.7)
	a2.plot(lam2_a5, flux2_a5, color="b", alpha=0.7)
	a2.plot(lam2_a6, flux2_a6, color="b", alpha=0.7)
	a2.plot(lam2_a7, flux2_a7, color="b", alpha=0.7)
	a2.plot(lam2_a8, flux2_a8, color="b", alpha=0.7)
	a2.plot(lam2_a9, flux2_a9, color="b", alpha=0.7)
	a2.plot(lam2_a10, flux2_a10, color="b", alpha=0.7)
	a2.plot(lam2_a11, flux2_a11, color="b", alpha=0.7)
	a2.plot(lam2_a12, flux2_a12, color="b", alpha=0.7)
	a2.plot(lam2_a13, flux2_a13, color="b", alpha=0.7)
	a2.plot(lam2_a14, flux2_a14, color="b", alpha=0.7)
	a2.plot(lam2_a15, flux2_a15, color="b", alpha=0.7)
	a2.plot(lam2_a16, flux2_a16, color="b", alpha=0.7)
	a2.plot(lam2_a17, flux2_a17, color="b", alpha=0.7)
	a2.plot(lam2_a18, flux2_a18, color="b", alpha=0.7)
	a2.plot(lam2_a19, flux2_a19, color="b", alpha=0.7)
	a2.plot(lam2_a20, flux2_a20, color="b", alpha=0.7)
	a2.plot(lam2_a21, flux2_a21, color="b", alpha=0.7)
	a2.plot(lam2_a22, flux2_a22, color="b", alpha=0.7)
	a2.plot(lam2_a23, flux2_a23, color="b", alpha=0.7)
	a2.plot(lam2_a24, flux2_a24, color="b", alpha=0.7)
	a2.plot(lam2_a25, flux2_a25, color="b", alpha=0.7)
	a2.plot(lam2_a26, flux2_a26, color="b", alpha=0.7)
	a2.plot(lam2_a27, flux2_a27, color="b", alpha=0.7)
	a2.plot(lam2_a28, flux2_a28, color="b", alpha=0.7)
	a2.plot(lam2_a29, flux2_a29, color="b", alpha=0.7)
	a2.plot(lam2_a30, flux2_a30, color="b", alpha=0.7)
	a2.plot(lam2_a31, flux2_a31, color="b", alpha=0.7)
	a2.plot(lam2_a32, flux2_a32, color="b", alpha=0.7)
	a2.plot(lam2_a33, flux2_a33, color="b", alpha=0.7)
	a2.plot(lam2_a34, flux2_a34, color="b", alpha=0.7)
	a2.plot(lam2_a35, flux2_a35, color="b", alpha=0.7)
	a2.plot(lam2_a36, flux2_a36, color="b", alpha=0.7)
	a2.plot(lam2_a37, flux2_a37, color="b", alpha=0.7)
	a2.plot(lam2_a38, flux2_a38, color="b", alpha=0.7)
	a2.plot(lam2_a39, flux2_a39, color="b", alpha=0.7)
	a2.plot(lam2_a40, flux2_a40, color="b", alpha=0.7)
	a2.plot(lam2_a41, flux2_a41, color="b", alpha=0.7)
	a2.plot(lam2_a42, flux2_a42, color="b", alpha=0.7)
	a2.plot(lam2_a43, flux2_a43, color="b", alpha=0.7)
	a2.plot(lam2_a44, flux2_a44, color="b", alpha=0.7)
	a2.plot(lam2_a45, flux2_a45, color="b", alpha=0.7)
	a2.plot(lam2_a46, flux2_a46, color="b", alpha=0.7)
	a2.plot(lam2_a47, flux2_a47, color="b", alpha=0.7)
	a2.plot(lam2_a48, flux2_a48, color="b", alpha=0.7)
	a2.plot(lam2_a49, flux2_a49, color="b", alpha=0.7)
	a2.plot(lam2_a50, flux2_a50, color="b", alpha=0.7)
	a2.plot(lam2_a51, flux2_a51, color="b", alpha=0.7)
	a2.plot(lam2_a52, flux2_a52, color="b", alpha=0.7)
	a2.plot(lam2_a53, flux2_a53, color="b", alpha=0.7)
	a2.plot(lam2_a54, flux2_a54, color="b", alpha=0.7)
	a2.plot(lam2_a55, flux2_a55, color="b", alpha=0.7)
	a2.plot(lam2_a56, flux2_a56, color="b", alpha=0.7)
	a2.plot(lam2_a57, flux2_a57, color="b", alpha=0.7)
	a2.plot(lam2_a58, flux2_a58, color="b", alpha=0.7)
	a2.plot(lam2_a59, flux2_a59, color="b", alpha=0.7)
	a2.plot(lam2_a60, flux2_a60, color="b", alpha=0.7)
	a2.plot(lam2_a61, flux2_a51, color="b", alpha=0.7)
	a2.plot(lam2_a62, flux2_a52, color="b", alpha=0.7)
	a2.plot(lam2_a63, flux2_a53, color="b", alpha=0.7)
	a2.plot(lam2_a64, flux2_a54, color="b", alpha=0.7)
	a2.plot(lam2_a65, flux2_a55, color="b", alpha=0.7)
	a2.plot(lam2_a66, flux2_a56, color="b", alpha=0.7)
	a2.plot(lam2_a67, flux2_a57, color="b", alpha=0.7)
	a2.plot(lam2_a68, flux2_a58, color="b", alpha=0.7)
	a2.plot(lam2_a69, flux2_a59, color="b", alpha=0.7)
	a2.plot(lam2_a70, flux2_a60, color="b", alpha=0.7)
	
	a3.plot(lam3_a1, flux3_a1, color="purple", alpha=0.7)
	a3.plot(lam3_a2, flux3_a2, color="purple", alpha=0.7)
	a3.plot(lam3_a3, flux3_a3, color="purple", alpha=0.7)
	a3.plot(lam3_a4, flux3_a4, color="purple", alpha=0.7)
	a3.plot(lam3_a5, flux3_a5, color="purple", alpha=0.7)
	a3.plot(lam3_a6, flux3_a6, color="purple", alpha=0.7)
	a3.plot(lam3_a7, flux3_a7, color="purple", alpha=0.7)
	a3.plot(lam3_a8, flux3_a8, color="purple", alpha=0.7)
	a3.plot(lam3_a9, flux3_a9, color="purple", alpha=0.7)
	a3.plot(lam3_a10, flux3_a10, color="purple", alpha=0.7)
	a3.plot(lam3_a11, flux3_a11, color="purple", alpha=0.7)
	a3.plot(lam3_a12, flux3_a12, color="purple", alpha=0.7)
	a3.plot(lam3_a13, flux3_a13, color="purple", alpha=0.7)
	a3.plot(lam3_a14, flux3_a14, color="purple", alpha=0.7)
	a3.plot(lam3_a15, flux3_a15, color="purple", alpha=0.7)
	a3.plot(lam3_a16, flux3_a16, color="purple", alpha=0.7)
	a3.plot(lam3_a17, flux3_a17, color="purple", alpha=0.7)
	a3.plot(lam3_a18, flux3_a18, color="purple", alpha=0.7)
	a3.plot(lam3_a19, flux3_a19, color="purple", alpha=0.7)
	a3.plot(lam3_a20, flux3_a20, color="purple", alpha=0.7)
	a3.plot(lam3_a21, flux3_a21, color="purple", alpha=0.7)
	a3.plot(lam3_a22, flux3_a22, color="purple", alpha=0.7)
	a3.plot(lam3_a23, flux3_a23, color="purple", alpha=0.7)
	a3.plot(lam3_a24, flux3_a24, color="purple", alpha=0.7)
	a3.plot(lam3_a25, flux3_a25, color="purple", alpha=0.7)
	a3.plot(lam3_a26, flux3_a26, color="purple", alpha=0.7)
	a3.plot(lam3_a27, flux3_a27, color="purple", alpha=0.7)
	a3.plot(lam3_a28, flux3_a28, color="purple", alpha=0.7)
	a3.plot(lam3_a29, flux3_a29, color="purple", alpha=0.7)
	a3.plot(lam3_a30, flux3_a30, color="purple", alpha=0.7)
	a3.plot(lam3_a31, flux3_a31, color="purple", alpha=0.7)
	a3.plot(lam3_a32, flux3_a32, color="purple", alpha=0.7)
	a3.plot(lam3_a33, flux3_a33, color="purple", alpha=0.7)
	a3.plot(lam3_a34, flux3_a34, color="purple", alpha=0.7)
	a3.plot(lam3_a35, flux3_a35, color="purple", alpha=0.7)
	a3.plot(lam3_a36, flux3_a36, color="purple", alpha=0.7)
	a3.plot(lam3_a37, flux3_a37, color="purple", alpha=0.7)
	a3.plot(lam3_a38, flux3_a38, color="purple", alpha=0.7)
	a3.plot(lam3_a39, flux3_a39, color="purple", alpha=0.7)
	a3.plot(lam3_a40, flux3_a40, color="purple", alpha=0.7)
	
	
	
	
	a1.plot(lam1_c1, flux1_c1, color="g", alpha=0.7)
	a1.plot(lam1_c2, flux1_c2, color="g", alpha=0.7)
	a1.plot(lam1_c3, flux1_c3, color="g", alpha=0.7)
	a1.plot(lam1_c4, flux1_c4, color="g", alpha=0.7)
	a1.plot(lam1_c5, flux1_c5, color="g", alpha=0.7)
	a1.plot(lam1_c6, flux1_c6, color="g", alpha=0.7)
	a1.plot(lam1_c7, flux1_c7, color="g", alpha=0.7)
	a1.plot(lam1_c8, flux1_c8, color="g", alpha=0.7)
	a1.plot(lam1_c9, flux1_c9, color="g", alpha=0.7)
	a1.plot(lam1_c10, flux1_c10, color="g", alpha=0.7)
	a1.plot(lam1_c11, flux1_c11, color="g", alpha=0.7)
	a1.plot(lam1_c12, flux1_c12, color="g", alpha=0.7)
	a1.plot(lam1_c13, flux1_c13, color="g", alpha=0.7)
	a1.plot(lam1_c14, flux1_c14, color="g", alpha=0.7)
	a1.plot(lam1_c15, flux1_c15, color="g", alpha=0.7)
	a1.plot(lam1_c16, flux1_c16, color="g", alpha=0.7)
	a1.plot(lam1_c17, flux1_c17, color="g", alpha=0.7)
	a1.plot(lam1_c18, flux1_c18, color="g", alpha=0.7)
	
	
	a2.plot(lam2_c1, flux2_c1, color="b", alpha=0.7)
	a2.plot(lam2_c2, flux2_c2, color="b", alpha=0.7)
	a2.plot(lam2_c3, flux2_c3, color="b", alpha=0.7)
	a2.plot(lam2_c4, flux2_c4, color="b", alpha=0.7)
	a2.plot(lam2_c5, flux2_c5, color="b", alpha=0.7)
	a2.plot(lam2_c6, flux2_c6, color="b", alpha=0.7)
	a2.plot(lam2_c7, flux2_c7, color="b", alpha=0.7)
	a2.plot(lam2_c8, flux2_c8, color="b", alpha=0.7)
	a2.plot(lam2_c9, flux2_c9, color="b", alpha=0.7)
	a2.plot(lam2_c10, flux2_c10, color="b", alpha=0.7)
	a2.plot(lam2_c11, flux2_c11, color="b", alpha=0.7)
	a2.plot(lam2_c12, flux2_c12, color="b", alpha=0.7)
	a2.plot(lam2_c13, flux2_c13, color="b", alpha=0.7)
	a2.plot(lam2_c14, flux2_c14, color="b", alpha=0.7)
	a2.plot(lam2_c15, flux2_c15, color="b", alpha=0.7)
	a2.plot(lam2_c16, flux2_c16, color="b", alpha=0.7)
	a2.plot(lam2_c17, flux2_c17, color="b", alpha=0.7)
	a2.plot(lam2_c18, flux2_c18, color="b", alpha=0.7)
	a2.plot(lam2_c19, flux2_c19, color="b", alpha=0.7)
	a2.plot(lam2_c20, flux2_c20, color="b", alpha=0.7)
	a2.plot(lam2_c21, flux2_c21, color="b", alpha=0.7)
	a2.plot(lam2_c22, flux2_c22, color="b", alpha=0.7)
	a2.plot(lam2_c23, flux2_c23, color="b", alpha=0.7)
	a2.plot(lam2_c24, flux2_c24, color="b", alpha=0.7)
	a2.plot(lam2_c25, flux2_c25, color="b", alpha=0.7)
	a2.plot(lam2_c26, flux2_c26, color="b", alpha=0.7)
	
	
	a3.plot(lam3_c1, flux3_c1, color="purple", alpha=0.7)
	a3.plot(lam3_c2, flux3_c2, color="purple", alpha=0.7)
	a3.plot(lam3_c3, flux3_c3, color="purple", alpha=0.7)
	a3.plot(lam3_c4, flux3_c4, color="purple", alpha=0.7)
	a3.plot(lam3_c5, flux3_c5, color="purple", alpha=0.7)
	a3.plot(lam3_c6, flux3_c6, color="purple", alpha=0.7)
	a3.plot(lam3_c7, flux3_c7, color="purple", alpha=0.7)
	a3.plot(lam3_c8, flux3_c8, color="purple", alpha=0.7)
	a3.plot(lam3_c9, flux3_c9, color="purple", alpha=0.7)
	a3.plot(lam3_c10, flux3_c10, color="purple", alpha=0.7)
	a3.plot(lam3_c11, flux3_c11, color="purple", alpha=0.7)
	a3.plot(lam3_c12, flux3_c12, color="purple", alpha=0.7)
	a3.plot(lam3_c13, flux3_c13, color="purple", alpha=0.7)
	a3.plot(lam3_c14, flux3_c14, color="purple", alpha=0.7)
	a3.plot(lam3_c15, flux3_c15, color="purple", alpha=0.7)
	a3.plot(lam3_c16, flux3_c16, color="purple", alpha=0.7)
	a3.plot(lam3_c17, flux3_c17, color="purple", alpha=0.7)
	a3.plot(lam3_c18, flux3_c18, color="purple", alpha=0.7)
	a3.plot(lam3_c19, flux3_c19, color="purple", alpha=0.7)
	a3.plot(lam3_c20, flux3_c20, color="purple", alpha=0.7)
	a3.plot(lam3_c21, flux3_c21, color="purple", alpha=0.7)
	a3.plot(lam3_c22, flux3_c22, color="purple", alpha=0.7)
	a3.plot(lam3_c23, flux3_c23, color="purple", alpha=0.7)
	a3.plot(lam3_c24, flux3_c24, color="purple", alpha=0.7)
	a3.plot(lam3_c25, flux3_c25, color="purple", alpha=0.7)
	a3.plot(lam3_c26, flux3_c26, color="purple", alpha=0.7)
	a3.plot(lam3_c27, flux3_c27, color="purple", alpha=0.7)
	a3.plot(lam3_c28, flux3_c28, color="purple", alpha=0.7)
	a3.plot(lam3_c29, flux3_c29, color="purple", alpha=0.7)
	a3.plot(lam3_c30, flux3_c30, color="purple", alpha=0.7)
	a3.plot(lam3_c31, flux3_c31, color="purple", alpha=0.7)
	a3.plot(lam3_c32, flux3_c32, color="purple", alpha=0.7)
	a3.plot(lam3_c33, flux3_c33, color="purple", alpha=0.7)
	a3.plot(lam3_c34, flux3_c34, color="purple", alpha=0.7)
	a3.plot(lam3_c35, flux3_c35, color="purple", alpha=0.7)
	
	
	
	a1.plot(lam1_n1, flux1_n1, color="g", alpha=0.7)
	a1.plot(lam1_n2, flux1_n2, color="g", alpha=0.7)
	a1.plot(lam1_n3, flux1_n3, color="g", alpha=0.7)
	a1.plot(lam1_n4, flux1_n4, color="g", alpha=0.7)
	a1.plot(lam1_n5, flux1_n5, color="g", alpha=0.7)
	a1.plot(lam1_n6, flux1_n6, color="g", alpha=0.7)
	a1.plot(lam1_n7, flux1_n7, color="g", alpha=0.7)
	a1.plot(lam1_n8, flux1_n8, color="g", alpha=0.7)
	a1.plot(lam1_n9, flux1_n9, color="g", alpha=0.7)
	a1.plot(lam1_n10, flux1_n10, color="g", alpha=0.7)
	a1.plot(lam1_n11, flux1_n11, color="g", alpha=0.7)
	a1.plot(lam1_n12, flux1_n12, color="g", alpha=0.7)
	a1.plot(lam1_n13, flux1_n13, color="g", alpha=0.7)
	a1.plot(lam1_n14, flux1_n14, color="g", alpha=0.7)
	a1.plot(lam1_n15, flux1_n15, color="g", alpha=0.7)
	a1.plot(lam1_n16, flux1_n16, color="g", alpha=0.7)
	a1.plot(lam1_n17, flux1_n17, color="g", alpha=0.7)
	a1.plot(lam1_n18, flux1_n18, color="g", alpha=0.7)
	a1.plot(lam1_n19, flux1_n19, color="g", alpha=0.7)
	a1.plot(lam1_n20, flux1_n20, color="g", alpha=0.7)
	a1.plot(lam1_n21, flux1_n21, color="g", alpha=0.7)
	a1.plot(lam1_n22, flux1_n22, color="g", alpha=0.7)
	a1.plot(lam1_n23, flux1_n23, color="g", alpha=0.7)
	a1.plot(lam1_n24, flux1_n24, color="g", alpha=0.7)
	a1.plot(lam1_n25, flux1_n25, color="g", alpha=0.7)
	

	a2.plot(lam2_n1, flux2_n1, color="b", alpha=0.7)
	a2.plot(lam2_n2, flux2_n2, color="b", alpha=0.7)
	a2.plot(lam2_n3, flux2_n3, color="b", alpha=0.7)
	a2.plot(lam2_n4, flux2_n4, color="b", alpha=0.7)
	a2.plot(lam2_n5, flux2_n5, color="b", alpha=0.7)
	a2.plot(lam2_n6, flux2_n6, color="b", alpha=0.7)
	a2.plot(lam2_n7, flux2_n7, color="b", alpha=0.7)
	a2.plot(lam2_n8, flux2_n8, color="b", alpha=0.7)
	a2.plot(lam2_n9, flux2_n9, color="b", alpha=0.7)
	a2.plot(lam2_n10, flux2_n10, color="b", alpha=0.7)
	a2.plot(lam2_n11, flux2_n11, color="b", alpha=0.7)
	a2.plot(lam2_n12, flux2_n12, color="b", alpha=0.7)
	a2.plot(lam2_n13, flux2_n13, color="b", alpha=0.7)
	a2.plot(lam2_n14, flux2_n14, color="b", alpha=0.7)
	a2.plot(lam2_n15, flux2_n15, color="b", alpha=0.7)
	a2.plot(lam2_n16, flux2_n16, color="b", alpha=0.7)
	a2.plot(lam2_n17, flux2_n17, color="b", alpha=0.7)
	a2.plot(lam2_n18, flux2_n18, color="b", alpha=0.7)
	a2.plot(lam2_n19, flux2_n19, color="b", alpha=0.7)
	a2.plot(lam2_n20, flux2_n20, color="b", alpha=0.7)
	a2.plot(lam2_n21, flux2_n21, color="b", alpha=0.7)
	
	
	a3.plot(lam3_n1, flux3_n1, color="purple", alpha=0.7)
	a3.plot(lam3_n2, flux3_n2, color="purple", alpha=0.7)
	a3.plot(lam3_n3, flux3_n3, color="purple", alpha=0.7)
	a3.plot(lam3_n4, flux3_n4, color="purple", alpha=0.7)
	a3.plot(lam3_n5, flux3_n5, color="purple", alpha=0.7)
	a3.plot(lam3_n6, flux3_n6, color="purple", alpha=0.7)
	a3.plot(lam3_n7, flux3_n7, color="purple", alpha=0.7)
	a3.plot(lam3_n8, flux3_n8, color="purple", alpha=0.7)
	a3.plot(lam3_n9, flux3_n9, color="purple", alpha=0.7)
	a3.plot(lam3_n10, flux3_n10, color="purple", alpha=0.7)
	a3.plot(lam3_n11, flux3_n11, color="purple", alpha=0.7)
	a3.plot(lam3_n12, flux3_n12, color="purple", alpha=0.7)
	a3.plot(lam3_n13, flux3_n13, color="purple", alpha=0.7)
	a3.plot(lam3_n14, flux3_n14, color="purple", alpha=0.7)
	a3.plot(lam3_n15, flux3_n15, color="purple", alpha=0.7)
	a3.plot(lam3_n16, flux3_n16, color="purple", alpha=0.7)
	a3.plot(lam3_n17, flux3_n17, color="purple", alpha=0.7)
	a3.plot(lam3_n18, flux3_n18, color="purple", alpha=0.7)
	a3.plot(lam3_n19, flux3_n19, color="purple", alpha=0.7)
	a3.plot(lam3_n20, flux3_n20, color="purple", alpha=0.7)
	a3.plot(lam3_n21, flux3_n21, color="purple", alpha=0.7)
	a3.plot(lam3_n22, flux3_n22, color="purple", alpha=0.7)
	a3.plot(lam3_n23, flux3_n23, color="purple", alpha=0.7)
	a3.plot(lam3_n24, flux3_n24, color="purple", alpha=0.7)
	a3.plot(lam3_n25, flux3_n25, color="purple", alpha=0.7)
	a3.plot(lam3_n26, flux3_n26, color="purple", alpha=0.7)
	a3.plot(lam3_n27, flux3_n27, color="purple", alpha=0.7)
	a3.plot(lam3_n28, flux3_n28, color="purple", alpha=0.7)
	a3.plot(lam3_n29, flux3_n29, color="purple", alpha=0.7)
	a3.plot(lam3_n30, flux3_n30, color="purple", alpha=0.7)
	a3.plot(lam3_n31, flux3_n31, color="purple", alpha=0.7)
	a3.plot(lam3_n32, flux3_n32, color="purple", alpha=0.7)
	a3.plot(lam3_n33, flux3_n33, color="purple", alpha=0.7)
	a3.plot(lam3_n34, flux3_n34, color="purple", alpha=0.7)
	a3.plot(lam3_n35, flux3_n35, color="purple", alpha=0.7)
	a3.plot(lam3_n36, flux3_n36, color="purple", alpha=0.7)
	a3.plot(lam3_n37, flux3_n37, color="purple", alpha=0.7)
	a3.plot(lam3_n38, flux3_n38, color="purple", alpha=0.7)
	
	
	
	
	a1.plot(lam1_s1, flux1_s1, color="g", alpha=0.7)
	a1.plot(lam1_s2, flux1_s2, color="g", alpha=0.7)
	a1.plot(lam1_s3, flux1_s3, color="g", alpha=0.7)
	a1.plot(lam1_s4, flux1_s4, color="g", alpha=0.7)
	a1.plot(lam1_s5, flux1_s5, color="g", alpha=0.7)
	a1.plot(lam1_s6, flux1_s6, color="g", alpha=0.7)
	a1.plot(lam1_s7, flux1_s7, color="g", alpha=0.7)
	a1.plot(lam1_s8, flux1_s8, color="g", alpha=0.7)
	a1.plot(lam1_s9, flux1_s9, color="g", alpha=0.7)
	a1.plot(lam1_s10, flux1_s10, color="g", alpha=0.7)
	a1.plot(lam1_s11, flux1_s11, color="g", alpha=0.7)
	a1.plot(lam1_s12, flux1_s12, color="g", alpha=0.7)
	a1.plot(lam1_s13, flux1_s13, color="g", alpha=0.7)
	a1.plot(lam1_s14, flux1_s14, color="g", alpha=0.7)
	a1.plot(lam1_s15, flux1_s15, color="g", alpha=0.7)
	a1.plot(lam1_s16, flux1_s16, color="g", alpha=0.7)
	a1.plot(lam1_s17, flux1_s17, color="g", alpha=0.7)
	a1.plot(lam1_s18, flux1_s18, color="g", alpha=0.7)
	a1.plot(lam1_s19, flux1_s19, color="g", alpha=0.7)
	a1.plot(lam1_s20, flux1_s20, color="g", alpha=0.7)
	a1.plot(lam1_s21, flux1_s21, color="g", alpha=0.7)
	a1.plot(lam1_s22, flux1_s22, color="g", alpha=0.7)
	a1.plot(lam1_s23, flux1_s23, color="g", alpha=0.7)
	a1.plot(lam1_s24, flux1_s24, color="g", alpha=0.7)
	

	a2.plot(lam2_s1, flux2_s1, color="b", alpha=0.7)
	a2.plot(lam2_s2, flux2_s2, color="b", alpha=0.7)
	a2.plot(lam2_s3, flux2_s3, color="b", alpha=0.7)
	a2.plot(lam2_s4, flux2_s4, color="b", alpha=0.7)
	a2.plot(lam2_s5, flux2_s5, color="b", alpha=0.7)
	a2.plot(lam2_s6, flux2_s6, color="b", alpha=0.7)
	a2.plot(lam2_s7, flux2_s7, color="b", alpha=0.7)
	a2.plot(lam2_s8, flux2_s8, color="b", alpha=0.7)
	a2.plot(lam2_s9, flux2_s9, color="b", alpha=0.7)
	a2.plot(lam2_s10, flux2_s10, color="b", alpha=0.7)
	a2.plot(lam2_s11, flux2_s11, color="b", alpha=0.7)
	a2.plot(lam2_s12, flux2_s12, color="b", alpha=0.7)
	a2.plot(lam2_s13, flux2_s13, color="b", alpha=0.7)
	a2.plot(lam2_s14, flux2_s14, color="b", alpha=0.7)
	a2.plot(lam2_s15, flux2_s15, color="b", alpha=0.7)
	a2.plot(lam2_s16, flux2_s16, color="b", alpha=0.7)
	a2.plot(lam2_s17, flux2_s17, color="b", alpha=0.7)
	a2.plot(lam2_s18, flux2_s18, color="b", alpha=0.7)
	a2.plot(lam2_s19, flux2_s19, color="b", alpha=0.7)
	a2.plot(lam2_s20, flux2_s20, color="b", alpha=0.7)
	a2.plot(lam2_s21, flux2_s21, color="b", alpha=0.7)
	a2.plot(lam2_s22, flux2_s22, color="b", alpha=0.7)
	a2.plot(lam2_s23, flux2_s23, color="b", alpha=0.7)
	a2.plot(lam2_s24, flux2_s24, color="b", alpha=0.7)
	a2.plot(lam2_s25, flux2_s25, color="b", alpha=0.7)
	a2.plot(lam2_s26, flux2_s26, color="b", alpha=0.7)
	a2.plot(lam2_s27, flux2_s27, color="b", alpha=0.7)
	a2.plot(lam2_s28, flux2_s28, color="b", alpha=0.7)
	a2.plot(lam2_s29, flux2_s29, color="b", alpha=0.7)
	a2.plot(lam2_s30, flux2_s30, color="b", alpha=0.7)
	a2.plot(lam2_s31, flux2_s31, color="b", alpha=0.7)
	a2.plot(lam2_s32, flux2_s32, color="b", alpha=0.7)
	a2.plot(lam2_s33, flux2_s33, color="b", alpha=0.7)
	a2.plot(lam2_s34, flux2_s34, color="b", alpha=0.7)
	a2.plot(lam2_s35, flux2_s35, color="b", alpha=0.7)
	a2.plot(lam2_s36, flux2_s36, color="b", alpha=0.7)
	a2.plot(lam2_s37, flux2_s37, color="b", alpha=0.7)
	a2.plot(lam2_s38, flux2_s38, color="b", alpha=0.7)
	a2.plot(lam2_s39, flux2_s39, color="b", alpha=0.7)
	a2.plot(lam2_s40, flux2_s40, color="b", alpha=0.7)
	a2.plot(lam2_s41, flux2_s41, color="b", alpha=0.7)
	a2.plot(lam2_s42, flux2_s42, color="b", alpha=0.7)
	
	
	a3.plot(lam3_s1, flux3_s1, color="purple", alpha=0.7)
	a3.plot(lam3_s2, flux3_s2, color="purple", alpha=0.7)
	a3.plot(lam3_s3, flux3_s3, color="purple", alpha=0.7)
	a3.plot(lam3_s4, flux3_s4, color="purple", alpha=0.7)
	a3.plot(lam3_s5, flux3_s5, color="purple", alpha=0.7)
	a3.plot(lam3_s6, flux3_s6, color="purple", alpha=0.7)
	a3.plot(lam3_s7, flux3_s7, color="purple", alpha=0.7)
	a3.plot(lam3_s8, flux3_s8, color="purple", alpha=0.7)
	a3.plot(lam3_s9, flux3_s9, color="purple", alpha=0.7)
	a3.plot(lam3_s10, flux3_s10, color="purple", alpha=0.7)
	a3.plot(lam3_s11, flux3_s11, color="purple", alpha=0.7)
	a3.plot(lam3_s12, flux3_s12, color="purple", alpha=0.7)
	a3.plot(lam3_s13, flux3_s13, color="purple", alpha=0.7)
	a3.plot(lam3_s14, flux3_s14, color="purple", alpha=0.7)
	a3.plot(lam3_s15, flux3_s15, color="purple", alpha=0.7)
	a3.plot(lam3_s16, flux3_s16, color="purple", alpha=0.7)
	a3.plot(lam3_s17, flux3_s17, color="purple", alpha=0.7)
	a3.plot(lam3_s18, flux3_s18, color="purple", alpha=0.7)
	a3.plot(lam3_s19, flux3_s19, color="purple", alpha=0.7)
	a3.plot(lam3_s20, flux3_s20, color="purple", alpha=0.7)
	a3.plot(lam3_s21, flux3_s21, color="purple", alpha=0.7)
	a3.plot(lam3_s22, flux3_s22, color="purple", alpha=0.7)
	a3.plot(lam3_s23, flux3_s23, color="purple", alpha=0.7)
	a3.plot(lam3_s24, flux3_s24, color="purple", alpha=0.7)
	a3.plot(lam3_s25, flux3_s25, color="purple", alpha=0.7)
	
	
	
	
	a1.plot(lam1_u1, flux1_u1, color="g", label="1 < z < 1.5", alpha=0.7)
	a1.plot(lam1_u2, flux1_u2, color="g", alpha=0.7)
	a1.plot(lam1_u3, flux1_u3, color="g", alpha=0.7)
	a1.plot(lam1_u4, flux1_u4, color="g", alpha=0.7)
	a1.plot(lam1_u5, flux1_u5, color="g", alpha=0.7)
	a1.plot(lam1_u6, flux1_u6, color="g", alpha=0.7)
	a1.plot(lam1_u7, flux1_u7, color="g", alpha=0.7)
	a1.plot(lam1_u8, flux1_u8, color="g", alpha=0.7)
	a1.plot(lam1_u9, flux1_u9, color="g", alpha=0.7)
	a1.plot(lam1_u10, flux1_u10, color="g", alpha=0.7)
	a1.plot(lam1_u11, flux1_u11, color="g", alpha=0.7)
	a1.plot(lam1_u12, flux1_u12, color="g", alpha=0.7)
	a1.plot(lam1_u13, flux1_u13, color="g", alpha=0.7)
	a1.plot(lam1_u14, flux1_u14, color="g", alpha=0.7)
	a1.plot(lam1_u15, flux1_u15, color="g", alpha=0.7)
	a1.plot(lam1_u16, flux1_u16, color="g", alpha=0.7)
	a1.plot(lam1_u17, flux1_u17, color="g", alpha=0.7)
	a1.plot(lam1_u18, flux1_u18, color="g", alpha=0.7)
	a1.plot(lam1_u19, flux1_u19, color="g", alpha=0.7)
	a1.plot(lam1_u20, flux1_u20, color="g", alpha=0.7)
	a1.plot(lam1_u21, flux1_u21, color="g", alpha=0.7)
	a1.plot(lam1_u22, flux1_u22, color="g", alpha=0.7)
	a1.plot(lam1_u23, flux1_u23, color="g", alpha=0.7)
	a1.plot(lam1_u24, flux1_u24, color="g", alpha=0.7)
	a1.plot(lam1_u25, flux1_u25, color="g", alpha=0.7)
	a1.plot(lam1_u26, flux1_u26, color="g", alpha=0.7)
	a1.plot(lam1_u27, flux1_u27, color="g", alpha=0.7)
	a1.plot(lam1_u28, flux1_u28, color="g", alpha=0.7)
	a1.plot(lam1_u29, flux1_u29, color="g", alpha=0.7)
	a1.plot(lam1_u30, flux1_u30, color="g", alpha=0.7)
	a1.plot(lam1_u31, flux1_u31, color="g", alpha=0.7)
	a1.plot(lam1_u32, flux1_u32, color="g", alpha=0.7)
	a1.plot(lam1_u33, flux1_u33, color="g", alpha=0.7)
	a1.plot(lam1_u34, flux1_u34, color="g", alpha=0.7)
	a1.plot(lam1_u35, flux1_u35, color="g", alpha=0.7)
	a1.plot(lam1_u36, flux1_u36, color="g", alpha=0.7)
	a1.plot(lam1_u37, flux1_u37, color="g", alpha=0.7)
	a1.plot(lam1_u38, flux1_u38, color="g", alpha=0.7)
	a1.plot(lam1_u39, flux1_u39, color="g", alpha=0.7)
	a1.plot(lam1_u40, flux1_u40, color="g", alpha=0.7)
	a1.plot(lam1_u41, flux1_u41, color="g", alpha=0.7)
	a1.plot(lam1_u42, flux1_u42, color="g", alpha=0.7)
	a1.plot(lam1_u43, flux1_u43, color="g", alpha=0.7)
	a1.plot(lam1_u44, flux1_u44, color="g", alpha=0.7)
	
	
	a2.plot(lam2_u1, flux2_u1, color="b", label="1.5 < z < 2", alpha=0.7)
	a2.plot(lam2_u2, flux2_u2, color="b", alpha=0.7)
	a2.plot(lam2_u3, flux2_u3, color="b", alpha=0.7)
	a2.plot(lam2_u4, flux2_u4, color="b", alpha=0.7)
	a2.plot(lam2_u5, flux2_u5, color="b", alpha=0.7)
	a2.plot(lam2_u6, flux2_u6, color="b", alpha=0.7)
	a2.plot(lam2_u7, flux2_u7, color="b", alpha=0.7)
	a2.plot(lam2_u8, flux2_u8, color="b", alpha=0.7)
	a2.plot(lam2_u9, flux2_u9, color="b", alpha=0.7)
	a2.plot(lam2_u10, flux2_u10, color="b", alpha=0.7)
	a2.plot(lam2_u11, flux2_u11, color="b", alpha=0.7)
	a2.plot(lam2_u12, flux2_u12, color="b", alpha=0.7)
	a2.plot(lam2_u13, flux2_u13, color="b", alpha=0.7)
	a2.plot(lam2_u14, flux2_u14, color="b", alpha=0.7)
	a2.plot(lam2_u15, flux2_u15, color="b", alpha=0.7)
	a2.plot(lam2_u16, flux2_u16, color="b", alpha=0.7)
	a2.plot(lam2_u17, flux2_u17, color="b", alpha=0.7)
	a2.plot(lam2_u18, flux2_u18, color="b", alpha=0.7)
	a2.plot(lam2_u19, flux2_u19, color="b", alpha=0.7)
	a2.plot(lam2_u20, flux2_u20, color="b", alpha=0.7)
	a2.plot(lam2_u21, flux2_u21, color="b", alpha=0.7)
	a2.plot(lam2_u22, flux2_u22, color="b", alpha=0.7)
	a2.plot(lam2_u23, flux2_u23, color="b", alpha=0.7)
	a2.plot(lam2_u24, flux2_u24, color="b", alpha=0.7)
	a2.plot(lam2_u25, flux2_u25, color="b", alpha=0.7)
	a2.plot(lam2_u26, flux2_u26, color="b", alpha=0.7)
	a2.plot(lam2_u27, flux2_u27, color="b", alpha=0.7)
	a2.plot(lam2_u28, flux2_u28, color="b", alpha=0.7)
	a2.plot(lam2_u29, flux2_u29, color="b", alpha=0.7)
	a2.plot(lam2_u30, flux2_u30, color="b", alpha=0.7)
	a2.plot(lam2_u31, flux2_u31, color="b", alpha=0.7)
	a2.plot(lam2_u32, flux2_u32, color="b", alpha=0.7)
	a2.plot(lam2_u33, flux2_u33, color="b", alpha=0.7)
	a2.plot(lam2_u34, flux2_u34, color="b", alpha=0.7)
	a2.plot(lam2_u35, flux2_u35, color="b", alpha=0.7)
	a2.plot(lam2_u36, flux2_u36, color="b", alpha=0.7)
	a2.plot(lam2_u37, flux2_u37, color="b", alpha=0.7)
	a2.plot(lam2_u38, flux2_u38, color="b", alpha=0.7)
	a2.plot(lam2_u39, flux2_u39, color="b", alpha=0.7)
	a2.plot(lam2_u40, flux2_u40, color="b", alpha=0.7)
	a2.plot(lam2_u41, flux2_u41, color="b", alpha=0.7)
	a2.plot(lam2_u42, flux2_u42, color="b", alpha=0.7)
	a2.plot(lam2_u43, flux2_u43, color="b", alpha=0.7)
	a2.plot(lam2_u44, flux2_u44, color="b", alpha=0.7)
	a2.plot(lam2_u45, flux2_u45, color="b", alpha=0.7)
	a2.plot(lam2_u46, flux2_u46, color="b", alpha=0.7)
	a2.plot(lam2_u47, flux2_u47, color="b", alpha=0.7)
	a2.plot(lam2_u48, flux2_u48, color="b", alpha=0.7)
	a2.plot(lam2_u49, flux2_u49, color="b", alpha=0.7)
	a2.plot(lam2_u50, flux2_u50, color="b", alpha=0.7)
	a2.plot(lam2_u51, flux2_u51, color="b", alpha=0.7)
	a2.plot(lam2_u52, flux2_u52, color="b", alpha=0.7)
	a2.plot(lam2_u53, flux2_u53, color="b", alpha=0.7)
	a2.plot(lam2_u54, flux2_u54, color="b", alpha=0.7)
	a2.plot(lam2_u55, flux2_u55, color="b", alpha=0.7)
	a2.plot(lam2_u56, flux2_u56, color="b", alpha=0.7)
	a2.plot(lam2_u57, flux2_u57, color="b", alpha=0.7)
	a2.plot(lam2_u58, flux2_u58, color="b", alpha=0.7)
	a2.plot(lam2_u59, flux2_u59, color="b", alpha=0.7)
	
	
	a3.plot(lam3_u1, flux3_u1, color="purple", label="2 < z < 2.5", alpha=0.7)
	a3.plot(lam3_u2, flux3_u2, color="purple", alpha=0.7)
	a3.plot(lam3_u3, flux3_u3, color="purple", alpha=0.7)
	a3.plot(lam3_u4, flux3_u4, color="purple", alpha=0.7)
	a3.plot(lam3_u5, flux3_u5, color="purple", alpha=0.7)
	a3.plot(lam3_u6, flux3_u6, color="purple", alpha=0.7)
	a3.plot(lam3_u7, flux3_u7, color="purple", alpha=0.7)
	a3.plot(lam3_u8, flux3_u8, color="purple", alpha=0.7)
	a3.plot(lam3_u9, flux3_u9, color="purple", alpha=0.7)
	a3.plot(lam3_u10, flux3_u10, color="purple", alpha=0.7)
	a3.plot(lam3_u11, flux3_u11, color="purple", alpha=0.7)
	a3.plot(lam3_u12, flux3_u12, color="purple", alpha=0.7)
	a3.plot(lam3_u13, flux3_u13, color="purple", alpha=0.7)
	a3.plot(lam3_u14, flux3_u14, color="purple", alpha=0.7)
	a3.plot(lam3_u15, flux3_u15, color="purple", alpha=0.7)
	a3.plot(lam3_u16, flux3_u16, color="purple", alpha=0.7)
	a3.plot(lam3_u17, flux3_u17, color="purple", alpha=0.7)
	a3.plot(lam3_u18, flux3_u18, color="purple", alpha=0.7)
	a3.plot(lam3_u19, flux3_u19, color="purple", alpha=0.7)
	a3.plot(lam3_u20, flux3_u20, color="purple", alpha=0.7)
	a3.plot(lam3_u21, flux3_u21, color="purple", alpha=0.7)
	a3.plot(lam3_u22, flux3_u22, color="purple", alpha=0.7)
	a3.plot(lam3_u23, flux3_u23, color="purple", alpha=0.7)
	a3.plot(lam3_u24, flux3_u24, color="purple", alpha=0.7)
	a3.plot(lam3_u25, flux3_u25, color="purple", alpha=0.7)
	a3.plot(lam3_u26, flux3_u26, color="purple", alpha=0.7)
	a3.plot(lam3_u27, flux3_u27, color="purple", alpha=0.7)
	a3.plot(lam3_u28, flux3_u28, color="purple", alpha=0.7)
	a3.plot(lam3_u29, flux3_u29, color="purple", alpha=0.7)
	a3.plot(lam3_u30, flux3_u30, color="purple", alpha=0.7)
	a3.plot(lam3_u31, flux3_u31, color="purple", alpha=0.7)
	a3.plot(lam3_u32, flux3_u32, color="purple", alpha=0.7)
	a3.plot(lam3_u33, flux3_u33, color="purple", alpha=0.7)
	a3.plot(lam3_u34, flux3_u34, color="purple", alpha=0.7)
	a3.plot(lam3_u35, flux3_u35, color="purple", alpha=0.7)
	a3.plot(lam3_u36, flux3_u36, color="purple", alpha=0.7)
	a3.plot(lam3_u37, flux3_u37, color="purple", alpha=0.7)
	a3.plot(lam3_u38, flux3_u38, color="purple", alpha=0.7)
	a3.plot(lam3_u39, flux3_u39, color="purple", alpha=0.7)
	a3.plot(lam3_u40, flux3_u40, color="purple", alpha=0.7)
	a3.plot(lam3_u41, flux3_u41, color="purple", alpha=0.7)
	a3.plot(lam3_u42, flux3_u42, color="purple", alpha=0.7)
	a3.plot(lam3_u43, flux3_u43, color="purple", alpha=0.7)
	a3.plot(lam3_u44, flux3_u44, color="purple", alpha=0.7)
	a3.plot(lam3_u45, flux3_u45, color="purple", alpha=0.7)
	a3.plot(lam3_u46, flux3_u46, color="purple", alpha=0.7)
	a3.plot(lam3_u47, flux3_u47, color="purple", alpha=0.7)
	a3.plot(lam3_u48, flux3_u48, color="purple", alpha=0.7)
	a3.plot(lam3_u49, flux3_u49, color="purple", alpha=0.7)
	a3.plot(lam3_u50, flux3_u50, color="purple", alpha=0.7)
	a3.plot(lam3_u51, flux3_u51, color="purple", alpha=0.7)
	a3.plot(lam3_u52, flux3_u52, color="purple", alpha=0.7)
	
	
	a1.legend(loc=2)
	a2.legend(loc=2)
	a3.legend(loc=2)
	a1.set_xlim([2000,140000])
	a1.set_ylim([0,2.4])
	a1.set_xscale("log")
	a2.set_xlim([2000,140000])
	a2.set_ylim([0,2.4])
	a2.set_xscale("log")
	a3.set_xlim([2000,140000])
	a3.set_ylim([0,2.4])
	a3.set_xscale("log")
	pylab.suptitle("Wavelength vs Flux", fontsize=19)
	fig.text(0.5, 0.01, "Wavelength", fontsize=18)
	fig.text(0.01, 0.5, "Flux", rotation = "vertical", fontsize=18)
	fig.subplots_adjust(hspace=0, wspace=0)
	
	os.chdir("/Volumes/TOSHIBA EXT/3d_hst")
	
	pylab.ion()
	pylab.show()
		










print "end"