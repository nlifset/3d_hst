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
import pandas as pd
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
	
	chunk1_a = data_fast_flagged1[(data_z_flagged1["z_peak"] < 1)]
	chunk2_a = data_fast_flagged1[(data_z_flagged1["z_peak"] >= 1) & (data_z_flagged1["z_peak"] < 1.5)]
	chunk3_a = data_fast_flagged1[(data_z_flagged1["z_peak"] >= 1.5) & (data_z_flagged1["z_peak"] < 2)]
	chunk4_a = data_fast_flagged1[(data_z_flagged1["z_peak"] >= 2) & (data_z_flagged1["z_peak"] <2.5)]
	chunk5_a = data_fast_flagged1[(data_z_flagged1["z_peak"] >= 2.5)]
	
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
	
	chunk_1_a = data_z_flagged1[(data_z_flagged1["z_peak"] < 1)]
	chunk_2_a = data_z_flagged1[(data_z_flagged1["z_peak"] >= 1) & (data_z_flagged1["z_peak"] < 1.5)]
	chunk_3_a = data_z_flagged1[(data_z_flagged1["z_peak"] >= 1.5) & (data_z_flagged1["z_peak"] < 2)]
	chunk_4_a = data_z_flagged1[(data_z_flagged1["z_peak"] >= 2) & (data_z_flagged1["z_peak"] <2.5)]
	chunk_5_a = data_z_flagged1[(data_z_flagged1["z_peak"] >= 2.5)]
	
	massed_1_a = chunk_1_a[(chunk1_a["lmass"] >= 11)]
	massed_2_a = chunk_2_a[(chunk2_a["lmass"] >= 11)]
	massed_3_a = chunk_3_a[(chunk3_a["lmass"] >= 11)]
	massed_4_a = chunk_4_a[(chunk4_a["lmass"] >= 11)]
	massed_5_a = chunk_5_a[(chunk5_a["lmass"] >= 11)]
	unmassed_a_ = data_z_flagged1[(data_fast_flagged1["lmass"] < 11)]
	
	z1_a = massed_1_a["z_peak"]
	z2_a = massed_2_a["z_peak"]
	z3_a = massed_3_a["z_peak"]
	z4_a = massed_4_a["z_peak"]
	z5_a = massed_5_a["z_peak"]
	unz_a = unmassed_a_["z_peak"]
	
	
	
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
	
	chunk1_c = data_fast_flagged2[(data_z_flagged2["z_peak"] < 1)]
	chunk2_c = data_fast_flagged2[(data_z_flagged2["z_peak"] >= 1) & (data_z_flagged2["z_peak"] < 1.5)]
	chunk3_c = data_fast_flagged2[(data_z_flagged2["z_peak"] >= 1.5) & (data_z_flagged2["z_peak"] < 2)]
	chunk4_c = data_fast_flagged2[(data_z_flagged2["z_peak"] >= 2) & (data_z_flagged2["z_peak"] <2.5)]
	chunk5_c = data_fast_flagged2[(data_z_flagged2["z_peak"] >= 2.5)]
	
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
	
	chunk_1_c = data_z_flagged2[(data_z_flagged2["z_peak"] < 1)]
	chunk_2_c = data_z_flagged2[(data_z_flagged2["z_peak"] >= 1) & (data_z_flagged2["z_peak"] < 1.5)]
	chunk_3_c = data_z_flagged2[(data_z_flagged2["z_peak"] >= 1.5) & (data_z_flagged2["z_peak"] < 2)]
	chunk_4_c = data_z_flagged2[(data_z_flagged2["z_peak"] >= 2) & (data_z_flagged2["z_peak"] <2.5)]
	chunk_5_c = data_z_flagged2[(data_z_flagged2["z_peak"] >= 2.5)]
	
	massed_1_c = chunk_1_c[(chunk1_c["lmass"] >= 11)]
	massed_2_c = chunk_2_c[(chunk2_c["lmass"] >= 11)]
	massed_3_c = chunk_3_c[(chunk3_c["lmass"] >= 11)]
	massed_4_c = chunk_4_c[(chunk4_c["lmass"] >= 11)]
	massed_5_c = chunk_5_c[(chunk5_c["lmass"] >= 11)]
	unmassed_c_ = data_z_flagged2[(data_fast_flagged2["lmass"] < 11)]
	
	z1_c = massed_1_c["z_peak"]
	z2_c = massed_2_c["z_peak"]
	z3_c = massed_3_c["z_peak"]
	z4_c = massed_4_c["z_peak"]
	z5_c = massed_5_c["z_peak"]
	unz_c = unmassed_c_["z_peak"]
	
	
	
	
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
	
	chunk1_n = data_fast_flagged3[(data_z_flagged3["z_peak"] < 1)]
	chunk2_n = data_fast_flagged3[(data_z_flagged3["z_peak"] >= 1) & (data_z_flagged3["z_peak"] < 1.5)]
	chunk3_n = data_fast_flagged3[(data_z_flagged3["z_peak"] >= 1.5) & (data_z_flagged3["z_peak"] < 2)]
	chunk4_n = data_fast_flagged3[(data_z_flagged3["z_peak"] >= 2) & (data_z_flagged3["z_peak"] <2.5)]
	chunk5_n = data_fast_flagged3[(data_z_flagged3["z_peak"] >= 2.5)]
	
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
	
	chunk_1_n = data_z_flagged3[(data_z_flagged3["z_peak"] < 1)]
	chunk_2_n = data_z_flagged3[(data_z_flagged3["z_peak"] >= 1) & (data_z_flagged3["z_peak"] < 1.5)]
	chunk_3_n = data_z_flagged3[(data_z_flagged3["z_peak"] >= 1.5) & (data_z_flagged3["z_peak"] < 2)]
	chunk_4_n = data_z_flagged3[(data_z_flagged3["z_peak"] >= 2) & (data_z_flagged3["z_peak"] <2.5)]
	chunk_5_n = data_z_flagged3[(data_z_flagged3["z_peak"] >= 2.5)]
	
	massed_1_n = chunk_1_n[(chunk1_n["lmass"] >= 11)]
	massed_2_n = chunk_2_n[(chunk2_n["lmass"] >= 11)]
	massed_3_n = chunk_3_n[(chunk3_n["lmass"] >= 11)]
	massed_4_n = chunk_4_n[(chunk4_n["lmass"] >= 11)]
	massed_5_n = chunk_5_n[(chunk5_n["lmass"] >= 11)]
	unmassed_n_ = data_z_flagged3[(data_fast_flagged3["lmass"] < 11)]
	
	z1_n = massed_1_n["z_peak"]
	z2_n = massed_2_n["z_peak"]
	z3_n = massed_3_n["z_peak"]
	z4_n = massed_4_n["z_peak"]
	z5_n = massed_5_n["z_peak"]
	unz_n = unmassed_n_["z_peak"]
	
	
	
	
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
	
	chunk1_s = data_fast_flagged4[(data_z_flagged4["z_peak"] < 1)]
	chunk2_s = data_fast_flagged4[(data_z_flagged4["z_peak"] >= 1) & (data_z_flagged4["z_peak"] < 1.5)]
	chunk3_s = data_fast_flagged4[(data_z_flagged4["z_peak"] >= 1.5) & (data_z_flagged4["z_peak"] < 2)]
	chunk4_s = data_fast_flagged4[(data_z_flagged4["z_peak"] >= 2) & (data_z_flagged4["z_peak"] <2.5)]
	chunk5_s = data_fast_flagged4[(data_z_flagged4["z_peak"] >= 2.5)]
	
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
	
	chunk_1_s = data_z_flagged4[(data_z_flagged4["z_peak"] < 1)]
	chunk_2_s = data_z_flagged4[(data_z_flagged4["z_peak"] >= 1) & (data_z_flagged4["z_peak"] < 1.5)]
	chunk_3_s = data_z_flagged4[(data_z_flagged4["z_peak"] >= 1.5) & (data_z_flagged4["z_peak"] < 2)]
	chunk_4_s = data_z_flagged4[(data_z_flagged4["z_peak"] >= 2) & (data_z_flagged4["z_peak"] <2.5)]
	chunk_5_s = data_z_flagged4[(data_z_flagged4["z_peak"] >= 2.5)]
	
	massed_1_s = chunk_1_s[(chunk1_s["lmass"] >= 11)]
	massed_2_s = chunk_2_s[(chunk2_s["lmass"] >= 11)]
	massed_3_s = chunk_3_s[(chunk3_s["lmass"] >= 11)]
	massed_4_s = chunk_4_s[(chunk4_s["lmass"] >= 11)]
	massed_5_s = chunk_5_s[(chunk5_s["lmass"] >= 11)]
	unmassed_s_ = data_z_flagged4[(data_fast_flagged4["lmass"] < 11)]
	
	z1_s = massed_1_s["z_peak"]
	z2_s = massed_2_s["z_peak"]
	z3_s = massed_3_s["z_peak"]
	z4_s = massed_4_s["z_peak"]
	z5_s = massed_5_s["z_peak"]
	unz_s = unmassed_s_["z_peak"]
	
	
	
	
	
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
	
	chunk1_u = data_fast_flagged5[(data_z_flagged5["z_peak"] < 1)]
	chunk2_u = data_fast_flagged5[(data_z_flagged5["z_peak"] >= 1) & (data_z_flagged5["z_peak"] < 1.5)]
	chunk3_u = data_fast_flagged5[(data_z_flagged5["z_peak"] >= 1.5) & (data_z_flagged5["z_peak"] < 2)]
	chunk4_u = data_fast_flagged5[(data_z_flagged5["z_peak"] >= 2) & (data_z_flagged5["z_peak"] <2.5)]
	chunk5_u = data_fast_flagged5[(data_z_flagged5["z_peak"] >= 2.5)]
	
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
	
	chunk_1_u = data_z_flagged5[(data_z_flagged5["z_peak"] < 1)]
	chunk_2_u = data_z_flagged5[(data_z_flagged5["z_peak"] >= 1) & (data_z_flagged5["z_peak"] < 1.5)]
	chunk_3_u = data_z_flagged5[(data_z_flagged5["z_peak"] >= 1.5) & (data_z_flagged5["z_peak"] < 2)]
	chunk_4_u = data_z_flagged5[(data_z_flagged5["z_peak"] >= 2) & (data_z_flagged5["z_peak"] <2.5)]
	chunk_5_u = data_z_flagged5[(data_z_flagged5["z_peak"] >= 2.5)]
	
	massed_1_u = chunk_1_u[(chunk1_u["lmass"] >= 11)]
	massed_2_u = chunk_2_u[(chunk2_u["lmass"] >= 11)]
	massed_3_u = chunk_3_u[(chunk3_u["lmass"] >= 11)]
	massed_4_u = chunk_4_u[(chunk4_u["lmass"] >= 11)]
	massed_5_u = chunk_5_u[(chunk5_u["lmass"] >= 11)]
	unmassed_u_ = data_z_flagged5[(data_fast_flagged5["lmass"] < 11)]
	
	z1_u = massed_1_u["z_peak"]
	z2_u = massed_2_u["z_peak"]
	z3_u = massed_3_u["z_peak"]
	z4_u = massed_4_u["z_peak"]
	z5_u = massed_5_u["z_peak"]
	unz_u = unmassed_u_["z_peak"]
	
	
	
	
	pylab.scatter(z1_a, lmass1_a, color="0.7", alpha=0.3, markeredgecolor="none", edgecolor="none")
	pylab.scatter(z2_a, lmass2_a, color="b", label="chunk1 points, 1 < z < 1.5", alpha=0.7, markeredgecolor="none", edgecolor="none")
	pylab.scatter(z3_a, lmass3_a, color="g", label="chunk2 points, 1.5 < z < 2", alpha=0.7, markeredgecolor="none", edgecolor="none")
	pylab.scatter(z4_a, lmass4_a, color="purple", label="chunk3 points, 2 < z < 2.5", alpha=0.7, markeredgecolor="none", edgecolor="none")
	pylab.scatter(z5_a, lmass5_a, color="0.7", alpha=0.3, markeredgecolor="none", edgecolor="none")
	pylab.scatter(unz_a, lunmass_a, color="0.7", alpha=0.3, markeredgecolor="none", edgecolor="none", label="unused points")
	
	pylab.scatter(z1_c, lmass1_c, color="0.7", alpha=0.3, markeredgecolor="none", edgecolor="none")
	pylab.scatter(z2_c, lmass2_c, color="b", alpha=0.7, markeredgecolor="none", edgecolor="none")
	pylab.scatter(z3_c, lmass3_c, color="g", alpha=0.7, markeredgecolor="none", edgecolor="none")
	pylab.scatter(z4_c, lmass4_c, color="purple", alpha=0.7, markeredgecolor="none", edgecolor="none")
	pylab.scatter(z5_c, lmass5_c, color="0.7", alpha=0.3, markeredgecolor="none", edgecolor="none")
	pylab.scatter(unz_c, lunmass_c, color="0.7", alpha=0.3, markeredgecolor="none", edgecolor="none")
	
	pylab.scatter(z1_n, lmass1_n, color="0.7", alpha=0.3, markeredgecolor="none", edgecolor="none")
	pylab.scatter(z2_n, lmass2_n, color="b", alpha=0.7, markeredgecolor="none", edgecolor="none")
	pylab.scatter(z3_n, lmass3_n, color="g", alpha=0.7, markeredgecolor="none", edgecolor="none")
	pylab.scatter(z4_n, lmass4_n, color="purple", alpha=0.7, markeredgecolor="none", edgecolor="none")
	pylab.scatter(z5_n, lmass5_n, color="0.7", alpha=0.3, markeredgecolor="none", edgecolor="none")
	pylab.scatter(unz_n, lunmass_n, color="0.7", alpha=0.3, markeredgecolor="none", edgecolor="none")
	
	pylab.scatter(z1_s, lmass1_s, color="0.7", alpha=0.3, markeredgecolor="none", edgecolor="none")
	pylab.scatter(z2_s, lmass2_s, color="b", alpha=0.7, markeredgecolor="none", edgecolor="none")
	pylab.scatter(z3_s, lmass3_s, color="g", alpha=0.7, markeredgecolor="none", edgecolor="none")
	pylab.scatter(z4_s, lmass4_s, color="purple", alpha=0.7, markeredgecolor="none", edgecolor="none")
	pylab.scatter(z5_s, lmass5_s, color="0.7", alpha=0.3, markeredgecolor="none", edgecolor="none")
	pylab.scatter(unz_s, lunmass_s, color="0.7", alpha=0.3, markeredgecolor="none", edgecolor="none")
	
	pylab.scatter(z1_u, lmass1_u, color="0.7", alpha=0.3, markeredgecolor="none", edgecolor="none")
	pylab.scatter(z2_u, lmass2_u, color="b", alpha=0.7, markeredgecolor="none", edgecolor="none")
	pylab.scatter(z3_u, lmass3_u, color="g", alpha=0.7, markeredgecolor="none", edgecolor="none")
	pylab.scatter(z4_u, lmass4_u, color="purple", alpha=0.7, markeredgecolor="none", edgecolor="none")
	pylab.scatter(z5_u, lmass5_u, color="0.7", alpha=0.3, markeredgecolor="none", edgecolor="none")
	pylab.scatter(unz_u, lunmass_u, color="0.7", alpha=0.3, markeredgecolor="none", edgecolor="none")
	
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
	
	chunk2_a = data_fast_flagged1[(data_z_flagged1["z_peak"] >= 1) & (data_z_flagged1["z_peak"] < 1.5)]
	chunk3_a = data_fast_flagged1[(data_z_flagged1["z_peak"] >= 1.5) & (data_z_flagged1["z_peak"] < 2)]
	chunk4_a = data_fast_flagged1[(data_z_flagged1["z_peak"] >= 2) & (data_z_flagged1["z_peak"] <2.5)]
	
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
	
	chunk2_c = data_fast_flagged2[(data_z_flagged2["z_peak"] >= 1) & (data_z_flagged2["z_peak"] < 1.5)]
	chunk3_c = data_fast_flagged2[(data_z_flagged2["z_peak"] >= 1.5) & (data_z_flagged2["z_peak"] < 2)]
	chunk4_c = data_fast_flagged2[(data_z_flagged2["z_peak"] >= 2) & (data_z_flagged2["z_peak"] <2.5)]
	
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
	
	chunk2_n = data_fast_flagged3[(data_z_flagged3["z_peak"] >= 1) & (data_z_flagged3["z_peak"] < 1.5)]
	chunk3_n = data_fast_flagged3[(data_z_flagged3["z_peak"] >= 1.5) & (data_z_flagged3["z_peak"] < 2)]
	chunk4_n = data_fast_flagged3[(data_z_flagged3["z_peak"] >= 2) & (data_z_flagged3["z_peak"] <2.5)]
	
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
	
	chunk2_s = data_fast_flagged4[(data_z_flagged4["z_peak"] >= 1) & (data_z_flagged4["z_peak"] < 1.5)]
	chunk3_s = data_fast_flagged4[(data_z_flagged4["z_peak"] >= 1.5) & (data_z_flagged4["z_peak"] < 2)]
	chunk4_s = data_fast_flagged4[(data_z_flagged4["z_peak"] >= 2) & (data_z_flagged4["z_peak"] <2.5)]
	
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
	
	chunk2_u = data_fast_flagged5[(data_z_flagged5["z_peak"] >= 1) & (data_z_flagged5["z_peak"] < 1.5)]
	chunk3_u = data_fast_flagged5[(data_z_flagged5["z_peak"] >= 1.5) & (data_z_flagged5["z_peak"] < 2)]
	chunk4_u = data_fast_flagged5[(data_z_flagged5["z_peak"] >= 2) & (data_z_flagged5["z_peak"] <2.5)]
	
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
	
	a1.scatter(lage2_a, lmass2_a, color="b", label="chunk1 points, 1 < z < 1.5", alpha=0.7, markeredgecolor="none")
	a2.scatter(lage3_a, lmass3_a, color="g", label="chunk2 points, 1.5 < z < 2", alpha=0.7, markeredgecolor="none")
	a3.scatter(lage4_a, lmass4_a, color="purple", label="chunk3 points, 2 < z < 2.5", alpha=0.7, markeredgecolor="none")
	
	a1.scatter(lage2_c, lmass2_c, color="b", alpha=0.7, markeredgecolor="none")
	a2.scatter(lage3_c, lmass3_c, color="g", alpha=0.7, markeredgecolor="none")
	a3.scatter(lage4_c, lmass4_c, color="purple", alpha=0.7, markeredgecolor="none")
	
	a1.scatter(lage2_n, lmass2_n, color="b", alpha=0.7, markeredgecolor="none")
	a2.scatter(lage3_n, lmass3_n, color="g", alpha=0.7, markeredgecolor="none")
	a3.scatter(lage4_n, lmass4_n, color="purple", alpha=0.7, markeredgecolor="none")
	
	a1.scatter(lage2_s, lmass2_s, color="b", alpha=0.7, markeredgecolor="none")
	a2.scatter(lage3_s, lmass3_s, color="g", alpha=0.7, markeredgecolor="none")
	a3.scatter(lage4_s, lmass4_s, color="purple", alpha=0.7, markeredgecolor="none")
	
	a1.scatter(lage2_u, lmass2_u, color="b", alpha=0.7, markeredgecolor="none")
	a2.scatter(lage3_u, lmass3_u, color="g", alpha=0.7, markeredgecolor="none")
	a3.scatter(lage4_u, lmass4_u, color="purple", alpha=0.7, markeredgecolor="none")
	

	
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
	
	chunk2_a = data_fast_flagged1[(data_z_flagged1["z_peak"] >= 1) & (data_z_flagged1["z_peak"] < 1.5)]
	chunk3_a = data_fast_flagged1[(data_z_flagged1["z_peak"] >= 1.5) & (data_z_flagged1["z_peak"] < 2)]
	chunk4_a = data_fast_flagged1[(data_z_flagged1["z_peak"] >= 2) & (data_z_flagged1["z_peak"] <2.5)]
	
	massed2_a = chunk2_a[(chunk2_a["lmass"] >= 11)]
	massed3_a = chunk3_a[(chunk3_a["lmass"] >= 11)]
	massed4_a = chunk4_a[(chunk4_a["lmass"] >= 11)]
	
	lmass2_a = massed2_a["lmass"]
	lmass3_a = massed3_a["lmass"]
	lmass4_a = massed4_a["lmass"]
	
	chunk_2_a = dataz_flagged1[(data_z_flagged1["z_peak"] >= 1) & (data_z_flagged1["z_peak"] < 1.5)]
	chunk_3_a = dataz_flagged1[(data_z_flagged1["z_peak"] >= 1.5) & (data_z_flagged1["z_peak"] < 2)]
	chunk_4_a = dataz_flagged1[(data_z_flagged1["z_peak"] >= 2) & (data_z_flagged1["z_peak"] <2.5)]
	
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
	
	chunk2_c = data_fast_flagged2[(data_z_flagged2["z_peak"] >= 1) & (data_z_flagged2["z_peak"] < 1.5)]
	chunk3_c = data_fast_flagged2[(data_z_flagged2["z_peak"] >= 1.5) & (data_z_flagged2["z_peak"] < 2)]
	chunk4_c = data_fast_flagged2[(data_z_flagged2["z_peak"] >= 2) & (data_z_flagged2["z_peak"] <2.5)]
	
	massed2_c = chunk2_c[(chunk2_c["lmass"] >= 11)]
	massed3_c = chunk3_c[(chunk3_c["lmass"] >= 11)]
	massed4_c = chunk4_c[(chunk4_c["lmass"] >= 11)]
	
	lmass2_c = massed2_c["lmass"]
	lmass3_c = massed3_c["lmass"]
	lmass4_c = massed4_c["lmass"]

	chunk_2_c = dataz_flagged2[(data_z_flagged2["z_peak"] >= 1) & (data_z_flagged2["z_peak"] < 1.5)]
	chunk_3_c = dataz_flagged2[(data_z_flagged2["z_peak"] >= 1.5) & (data_z_flagged2["z_peak"] < 2)]
	chunk_4_c = dataz_flagged2[(data_z_flagged2["z_peak"] >= 2) & (data_z_flagged2["z_peak"] <2.5)]
	
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
	
	chunk2_n = data_fast_flagged3[(data_z_flagged3["z_peak"] >= 1) & (data_z_flagged3["z_peak"] < 1.5)]
	chunk3_n = data_fast_flagged3[(data_z_flagged3["z_peak"] >= 1.5) & (data_z_flagged3["z_peak"] < 2)]
	chunk4_n = data_fast_flagged3[(data_z_flagged3["z_peak"] >= 2) & (data_z_flagged3["z_peak"] <2.5)]
	
	massed2_n = chunk2_n[(chunk2_n["lmass"] >= 11)]
	massed3_n = chunk3_n[(chunk3_n["lmass"] >= 11)]
	massed4_n = chunk4_n[(chunk4_n["lmass"] >= 11)]
	
	lmass2_n = massed2_n["lmass"]
	lmass3_n = massed3_n["lmass"]
	lmass4_n = massed4_n["lmass"]

	chunk_2_n = dataz_flagged3[(data_z_flagged3["z_peak"] >= 1) & (data_z_flagged3["z_peak"] < 1.5)]
	chunk_3_n = dataz_flagged3[(data_z_flagged3["z_peak"] >= 1.5) & (data_z_flagged3["z_peak"] < 2)]
	chunk_4_n = dataz_flagged3[(data_z_flagged3["z_peak"] >= 2) & (data_z_flagged3["z_peak"] <2.5)]
	
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
	
	chunk2_s = data_fast_flagged4[(data_z_flagged4["z_peak"] >= 1) & (data_z_flagged4["z_peak"] < 1.5)]
	chunk3_s = data_fast_flagged4[(data_z_flagged4["z_peak"] >= 1.5) & (data_z_flagged4["z_peak"] < 2)]
	chunk4_s = data_fast_flagged4[(data_z_flagged4["z_peak"] >= 2) & (data_z_flagged4["z_peak"] <2.5)]
	
	massed2_s = chunk2_s[(chunk2_s["lmass"] >= 11)]
	massed3_s = chunk3_s[(chunk3_s["lmass"] >= 11)]
	massed4_s = chunk4_s[(chunk4_s["lmass"] >= 11)]
	
	lmass2_s = massed2_s["lmass"]
	lmass3_s = massed3_s["lmass"]
	lmass4_s = massed4_s["lmass"]

	chunk_2_s = dataz_flagged4[(data_z_flagged4["z_peak"] >= 1) & (data_z_flagged4["z_peak"] < 1.5)]
	chunk_3_s = dataz_flagged4[(data_z_flagged4["z_peak"] >= 1.5) & (data_z_flagged4["z_peak"] < 2)]
	chunk_4_s = dataz_flagged4[(data_z_flagged4["z_peak"] >= 2) & (data_z_flagged4["z_peak"] <2.5)]
	
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
	
	chunk2_u = data_fast_flagged5[(data_z_flagged5["z_peak"] >= 1) & (data_z_flagged5["z_peak"] < 1.5)]
	chunk3_u = data_fast_flagged5[(data_z_flagged5["z_peak"] >= 1.5) & (data_z_flagged5["z_peak"] < 2)]
	chunk4_u = data_fast_flagged5[(data_z_flagged5["z_peak"] >= 2) & (data_z_flagged5["z_peak"] <2.5)]
	
	massed2_u = chunk2_u[(chunk2_u["lmass"] >= 11)]
	massed3_u = chunk3_u[(chunk3_u["lmass"] >= 11)]
	massed4_u = chunk4_u[(chunk4_u["lmass"] >= 11)]
	
	lmass2_u = massed2_u["lmass"]
	lmass3_u = massed3_u["lmass"]
	lmass4_u = massed4_u["lmass"]

	chunk_2_u = dataz_flagged5[(data_z_flagged5["z_peak"] >= 1) & (data_z_flagged5["z_peak"] < 1.5)]
	chunk_3_u = dataz_flagged5[(data_z_flagged5["z_peak"] >= 1.5) & (data_z_flagged5["z_peak"] < 2)]
	chunk_4_u = dataz_flagged5[(data_z_flagged5["z_peak"] >= 2) & (data_z_flagged5["z_peak"] <2.5)]
	
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
	
	a1.scatter(n2_a, lmass2_a, color="b", label="chunk1 points, 1 < z < 1.5", alpha=0.7, markeredgecolor="none")
	a2.scatter(n3_a, lmass3_a, color="g", label="chunk2 points, 1.5 < z < 2", alpha=0.7, markeredgecolor="none")
	a3.scatter(n4_a, lmass4_a, color="purple", label="chunk3 points, 2 < z < 2.5", alpha=0.7, markeredgecolor="none")
	
	a1.scatter(n2_c, lmass2_c, color="b", alpha=0.7, markeredgecolor="none")
	a2.scatter(n3_c, lmass3_c, color="g", alpha=0.7, markeredgecolor="none")
	a3.scatter(n4_c, lmass4_c, color="purple", alpha=0.7, markeredgecolor="none")
	
	a1.scatter(n2_n, lmass2_n, color="b", alpha=0.7, markeredgecolor="none")
	a2.scatter(n3_n, lmass3_n, color="g", alpha=0.7, markeredgecolor="none")
	a3.scatter(n4_n, lmass4_n, color="purple", alpha=0.7, markeredgecolor="none")
	
	a1.scatter(n2_s, lmass2_s, color="b", alpha=0.7, markeredgecolor="none")
	a2.scatter(n3_s, lmass3_s, color="g", alpha=0.7, markeredgecolor="none")
	a3.scatter(n4_s, lmass4_s, color="purple", alpha=0.7, markeredgecolor="none")
	
	a1.scatter(n2_u, lmass2_u, color="b", alpha=0.7, markeredgecolor="none")
	a2.scatter(n3_u, lmass3_u, color="g", alpha=0.7, markeredgecolor="none")
	a3.scatter(n4_u, lmass4_u, color="purple", alpha=0.7, markeredgecolor="none")
	

	
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
	
	chunk2_a = data_fast_flagged1[(data_z_flagged1["z_peak"] >= 1) & (data_z_flagged1["z_peak"] < 1.5)]
	chunk3_a = data_fast_flagged1[(data_z_flagged1["z_peak"] >= 1.5) & (data_z_flagged1["z_peak"] < 2)]
	chunk4_a = data_fast_flagged1[(data_z_flagged1["z_peak"] >= 2) & (data_z_flagged1["z_peak"] <2.5)]
	
	massed2_a = chunk2_a[(chunk2_a["lmass"] >= 11)]
	massed3_a = chunk3_a[(chunk3_a["lmass"] >= 11)]
	massed4_a = chunk4_a[(chunk4_a["lmass"] >= 11)]
	
	lmass2_a = massed2_a["lmass"]
	lmass3_a = massed3_a["lmass"]
	lmass4_a = massed4_a["lmass"]
	
	chunk_2_a = dataz_flagged1[(data_z_flagged1["z_peak"] >= 1) & (data_z_flagged1["z_peak"] < 1.5)]
	chunk_3_a = dataz_flagged1[(data_z_flagged1["z_peak"] >= 1.5) & (data_z_flagged1["z_peak"] < 2)]
	chunk_4_a = dataz_flagged1[(data_z_flagged1["z_peak"] >= 2) & (data_z_flagged1["z_peak"] <2.5)]
	
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
	
	chunk2_c = data_fast_flagged2[(data_z_flagged2["z_peak"] >= 1) & (data_z_flagged2["z_peak"] < 1.5)]
	chunk3_c = data_fast_flagged2[(data_z_flagged2["z_peak"] >= 1.5) & (data_z_flagged2["z_peak"] < 2)]
	chunk4_c = data_fast_flagged2[(data_z_flagged2["z_peak"] >= 2) & (data_z_flagged2["z_peak"] <2.5)]
	
	massed2_c = chunk2_c[(chunk2_c["lmass"] >= 11)]
	massed3_c = chunk3_c[(chunk3_c["lmass"] >= 11)]
	massed4_c = chunk4_c[(chunk4_c["lmass"] >= 11)]
	
	lmass2_c = massed2_c["lmass"]
	lmass3_c = massed3_c["lmass"]
	lmass4_c = massed4_c["lmass"]

	chunk_2_c = dataz_flagged2[(data_z_flagged2["z_peak"] >= 1) & (data_z_flagged2["z_peak"] < 1.5)]
	chunk_3_c = dataz_flagged2[(data_z_flagged2["z_peak"] >= 1.5) & (data_z_flagged2["z_peak"] < 2)]
	chunk_4_c = dataz_flagged2[(data_z_flagged2["z_peak"] >= 2) & (data_z_flagged2["z_peak"] <2.5)]
	
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
	
	chunk2_n = data_fast_flagged3[(data_z_flagged3["z_peak"] >= 1) & (data_z_flagged3["z_peak"] < 1.5)]
	chunk3_n = data_fast_flagged3[(data_z_flagged3["z_peak"] >= 1.5) & (data_z_flagged3["z_peak"] < 2)]
	chunk4_n = data_fast_flagged3[(data_z_flagged3["z_peak"] >= 2) & (data_z_flagged3["z_peak"] <2.5)]
	
	massed2_n = chunk2_n[(chunk2_n["lmass"] >= 11)]
	massed3_n = chunk3_n[(chunk3_n["lmass"] >= 11)]
	massed4_n = chunk4_n[(chunk4_n["lmass"] >= 11)]
	
	lmass2_n = massed2_n["lmass"]
	lmass3_n = massed3_n["lmass"]
	lmass4_n = massed4_n["lmass"]

	chunk_2_n = dataz_flagged3[(data_z_flagged3["z_peak"] >= 1) & (data_z_flagged3["z_peak"] < 1.5)]
	chunk_3_n = dataz_flagged3[(data_z_flagged3["z_peak"] >= 1.5) & (data_z_flagged3["z_peak"] < 2)]
	chunk_4_n = dataz_flagged3[(data_z_flagged3["z_peak"] >= 2) & (data_z_flagged3["z_peak"] <2.5)]
	
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
	
	chunk2_s = data_fast_flagged4[(data_z_flagged4["z_peak"] >= 1) & (data_z_flagged4["z_peak"] < 1.5)]
	chunk3_s = data_fast_flagged4[(data_z_flagged4["z_peak"] >= 1.5) & (data_z_flagged4["z_peak"] < 2)]
	chunk4_s = data_fast_flagged4[(data_z_flagged4["z_peak"] >= 2) & (data_z_flagged4["z_peak"] <2.5)]
	
	massed2_s = chunk2_s[(chunk2_s["lmass"] >= 11)]
	massed3_s = chunk3_s[(chunk3_s["lmass"] >= 11)]
	massed4_s = chunk4_s[(chunk4_s["lmass"] >= 11)]
	
	lmass2_s = massed2_s["lmass"]
	lmass3_s = massed3_s["lmass"]
	lmass4_s = massed4_s["lmass"]

	chunk_2_s = dataz_flagged4[(data_z_flagged4["z_peak"] >= 1) & (data_z_flagged4["z_peak"] < 1.5)]
	chunk_3_s = dataz_flagged4[(data_z_flagged4["z_peak"] >= 1.5) & (data_z_flagged4["z_peak"] < 2)]
	chunk_4_s = dataz_flagged4[(data_z_flagged4["z_peak"] >= 2) & (data_z_flagged4["z_peak"] <2.5)]
	
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
	
	chunk2_u = data_fast_flagged5[(data_z_flagged5["z_peak"] >= 1) & (data_z_flagged5["z_peak"] < 1.5)]
	chunk3_u = data_fast_flagged5[(data_z_flagged5["z_peak"] >= 1.5) & (data_z_flagged5["z_peak"] < 2)]
	chunk4_u = data_fast_flagged5[(data_z_flagged5["z_peak"] >= 2) & (data_z_flagged5["z_peak"] <2.5)]
	
	massed2_u = chunk2_u[(chunk2_u["lmass"] >= 11)]
	massed3_u = chunk3_u[(chunk3_u["lmass"] >= 11)]
	massed4_u = chunk4_u[(chunk4_u["lmass"] >= 11)]
	
	lmass2_u = massed2_u["lmass"]
	lmass3_u = massed3_u["lmass"]
	lmass4_u = massed4_u["lmass"]

	chunk_2_u = dataz_flagged5[(data_z_flagged5["z_peak"] >= 1) & (data_z_flagged5["z_peak"] < 1.5)]
	chunk_3_u = dataz_flagged5[(data_z_flagged5["z_peak"] >= 1.5) & (data_z_flagged5["z_peak"] < 2)]
	chunk_4_u = dataz_flagged5[(data_z_flagged5["z_peak"] >= 2) & (data_z_flagged5["z_peak"] <2.5)]
	
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
	
	
	a1.scatter(n2_a, lmass2_a, color="b", label="chunk1 points, 1 < z < 1.5", alpha=0.7, markeredgecolor="none")
	a2.scatter(n3_a, lmass3_a, color="g", label="chunk2 points, 1.5 < z < 2", alpha=0.7, markeredgecolor="none")
	a3.scatter(n4_a, lmass4_a, color="purple", label="chunk3 points, 2 < z < 2.5", alpha=0.7, markeredgecolor="none")
	
	a1.scatter(n2_c, lmass2_c, color="b", alpha=0.7, markeredgecolor="none")
	a2.scatter(n3_c, lmass3_c, color="g", alpha=0.7, markeredgecolor="none")
	a3.scatter(n4_c, lmass4_c, color="purple", alpha=0.7, markeredgecolor="none")
	
	a1.scatter(n2_n, lmass2_n, color="b", alpha=0.7, markeredgecolor="none")
	a2.scatter(n3_n, lmass3_n, color="g", alpha=0.7, markeredgecolor="none")
	a3.scatter(n4_n, lmass4_n, color="purple", alpha=0.7, markeredgecolor="none")
	
	a1.scatter(n2_s, lmass2_s, color="b", alpha=0.7, markeredgecolor="none")
	a2.scatter(n3_s, lmass3_s, color="g", alpha=0.7, markeredgecolor="none")
	a3.scatter(n4_s, lmass4_s, color="purple", alpha=0.7, markeredgecolor="none")
	
	a1.scatter(n2_u, lmass2_u, color="b", alpha=0.7, markeredgecolor="none")
	a2.scatter(n3_u, lmass3_u, color="g", alpha=0.7, markeredgecolor="none")
	a3.scatter(n4_u, lmass4_u, color="purple", alpha=0.7, markeredgecolor="none")
	

	
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
	
	chunk2_a = data_fast_flagged1[(data_z_flagged1["z_peak"] >= 1) & (data_z_flagged1["z_peak"] < 1.5)]
	chunk3_a = data_fast_flagged1[(data_z_flagged1["z_peak"] >= 1.5) & (data_z_flagged1["z_peak"] < 2)]
	chunk4_a = data_fast_flagged1[(data_z_flagged1["z_peak"] >= 2) & (data_z_flagged1["z_peak"] <2.5)]
	
	massed2_a = chunk2_a[(chunk2_a["lmass"] >= 11)]
	massed3_a = chunk3_a[(chunk3_a["lmass"] >= 11)]
	massed4_a = chunk4_a[(chunk4_a["lmass"] >= 11)]
	
	lmass2_a = massed2_a["lmass"]
	lmass3_a = massed3_a["lmass"]
	lmass4_a = massed4_a["lmass"]
	
	chunk_2_a = data_z_flagged1[(data_z_flagged1["z_peak"] >= 1) & (data_z_flagged1["z_peak"] < 1.5)]
	chunk_3_a = data_z_flagged1[(data_z_flagged1["z_peak"] >= 1.5) & (data_z_flagged1["z_peak"] < 2)]
	chunk_4_a = data_z_flagged1[(data_z_flagged1["z_peak"] >= 2) & (data_z_flagged1["z_peak"] <2.5)]
	
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
	
	chunk2_c = data_fast_flagged2[(data_z_flagged2["z_peak"] >= 1) & (data_z_flagged2["z_peak"] < 1.5)]
	chunk3_c = data_fast_flagged2[(data_z_flagged2["z_peak"] >= 1.5) & (data_z_flagged2["z_peak"] < 2)]
	chunk4_c = data_fast_flagged2[(data_z_flagged2["z_peak"] >= 2) & (data_z_flagged2["z_peak"] <2.5)]
	
	massed2_c = chunk2_c[(chunk2_c["lmass"] >= 11)]
	massed3_c = chunk3_c[(chunk3_c["lmass"] >= 11)]
	massed4_c = chunk4_c[(chunk4_c["lmass"] >= 11)]
	
	lmass2_c = massed2_c["lmass"]
	lmass3_c = massed3_c["lmass"]
	lmass4_c = massed4_c["lmass"]
	
	chunk_2_c = data_z_flagged2[(data_z_flagged2["z_peak"] >= 1) & (data_z_flagged2["z_peak"] < 1.5)]
	chunk_3_c = data_z_flagged2[(data_z_flagged2["z_peak"] >= 1.5) & (data_z_flagged2["z_peak"] < 2)]
	chunk_4_c = data_z_flagged2[(data_z_flagged2["z_peak"] >= 2) & (data_z_flagged2["z_peak"] <2.5)]
	
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
	
	chunk2_n = data_fast_flagged3[(data_z_flagged3["z_peak"] >= 1) & (data_z_flagged3["z_peak"] < 1.5)]
	chunk3_n = data_fast_flagged3[(data_z_flagged3["z_peak"] >= 1.5) & (data_z_flagged3["z_peak"] < 2)]
	chunk4_n = data_fast_flagged3[(data_z_flagged3["z_peak"] >= 2) & (data_z_flagged3["z_peak"] <2.5)]
	
	massed2_n = chunk2_n[(chunk2_n["lmass"] >= 11)]
	massed3_n = chunk3_n[(chunk3_n["lmass"] >= 11)]
	massed4_n = chunk4_n[(chunk4_n["lmass"] >= 11)]
	
	lmass2_n = massed2_n["lmass"]
	lmass3_n = massed3_n["lmass"]
	lmass4_n = massed4_n["lmass"]
	
	chunk_2_n = data_z_flagged3[(data_z_flagged3["z_peak"] >= 1) & (data_z_flagged3["z_peak"] < 1.5)]
	chunk_3_n = data_z_flagged3[(data_z_flagged3["z_peak"] >= 1.5) & (data_z_flagged3["z_peak"] < 2)]
	chunk_4_n = data_z_flagged3[(data_z_flagged3["z_peak"] >= 2) & (data_z_flagged3["z_peak"] <2.5)]
	
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
	
	chunk2_s = data_fast_flagged4[(data_z_flagged4["z_peak"] >= 1) & (data_fast_flagged4["z_peak"] < 1.5)]
	chunk3_s = data_fast_flagged4[(data_fast_flagged4["z_peak"] >= 1.5) & (data_z_flagged4["z_peak"] < 2)]
	chunk4_s = data_fast_flagged4[(data_z_flagged4["z_peak"] >= 2) & (data_z_flagged4["z_peak"] <2.5)]
	
	massed2_s = chunk2_s[(chunk2_s["lmass"] >= 11)]
	massed3_s = chunk3_s[(chunk3_s["lmass"] >= 11)]
	massed4_s = chunk4_s[(chunk4_s["lmass"] >= 11)]
	
	lmass2_s = massed2_s["lmass"]
	lmass3_s = massed3_s["lmass"]
	lmass4_s = massed4_s["lmass"]
	
	chunk_2_s = data_z_flagged4[(data_z_flagged4["z_peak"] >= 1) & (data_fast_flagged4["z_peak"] < 1.5)]
	chunk_3_s = data_z_flagged4[(data_fast_flagged4["z_peak"] >= 1.5) & (data_z_flagged4["z_peak"] < 2)]
	chunk_4_s = data_z_flagged4[(data_z_flagged4["z_peak"] >= 2) & (data_z_flagged4["z_peak"] <2.5)]
	
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
	
	chunk2_u = data_fast_flagged5[(data_z_flagged5["z_peak"] >= 1) & (data_z_flagged5["z_peak"] < 1.5)]
	chunk3_u = data_fast_flagged5[(data_z_flagged5["z_peak"] >= 1.5) & (data_z_flagged5["z_peak"] < 2)]
	chunk4_u = data_fast_flagged5[(data_z_flagged5["z_peak"] >= 2) & (data_z_flagged5["z_peak"] <2.5)]
	
	massed2_u = chunk2_u[(chunk2_u["lmass"] >= 11)]
	massed3_u = chunk3_u[(chunk3_u["lmass"] >= 11)]
	massed4_u = chunk4_u[(chunk4_u["lmass"] >= 11)]
	
	lmass2_u = massed2_u["lmass"]
	lmass3_u = massed3_u["lmass"]
	lmass4_u = massed4_u["lmass"]
	
	chunk_2_u = data_z_flagged5[(data_z_flagged5["z_peak"] >= 1) & (data_z_flagged5["z_peak"] < 1.5)]
	chunk_3_u = data_z_flagged5[(data_z_flagged5["z_peak"] >= 1.5) & (data_z_flagged5["z_peak"] < 2)]
	chunk_4_u = data_z_flagged5[(data_z_flagged5["z_peak"] >= 2) & (data_z_flagged5["z_peak"] <2.5)]
	
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
	
	a1.scatter(lage2_a, lmass2_a, color="b", label="chunk1 points, 1 < z < 1.5", alpha=0.7, markeredgecolor="none")
	a2.scatter(lage3_a, lmass3_a, color="g", label="chunk2 points, 1.5 < z < 2", alpha=0.7, markeredgecolor="none")
	a3.scatter(lage4_a, lmass4_a, color="purple", label="chunk3 points, 2 < z < 2.5", alpha=0.7, markeredgecolor="none")
	
	a1.scatter(lage2_c, lmass2_c, color="b", alpha=0.7, markeredgecolor="none")
	a2.scatter(lage3_c, lmass3_c, color="g", alpha=0.7, markeredgecolor="none")
	a3.scatter(lage4_c, lmass4_c, color="purple", alpha=0.7, markeredgecolor="none")
	
	a1.scatter(lage2_n, lmass2_n, color="b", alpha=0.7, markeredgecolor="none")
	a2.scatter(lage3_n, lmass3_n, color="g", alpha=0.7, markeredgecolor="none")
	a3.scatter(lage4_n, lmass4_n, color="purple", alpha=0.7, markeredgecolor="none")
	
	a1.scatter(lage2_s, lmass2_s, color="b", alpha=0.7, markeredgecolor="none")
	a2.scatter(lage3_s, lmass3_s, color="g", alpha=0.7, markeredgecolor="none")
	a3.scatter(lage4_s, lmass4_s, color="purple", alpha=0.7, markeredgecolor="none")
	
	a1.scatter(lage2_u, lmass2_u, color="b", alpha=0.7, markeredgecolor="none")
	a2.scatter(lage3_u, lmass3_u, color="g", alpha=0.7, markeredgecolor="none")
	a3.scatter(lage4_u, lmass4_u, color="purple", alpha=0.7, markeredgecolor="none")
	

	
	
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
	
	chunk2_a = data_fast_flagged1[(data_z_flagged1["z_peak"] >= 1) & (data_z_flagged1["z_peak"] < 1.5)]
	chunk3_a = data_fast_flagged1[(data_z_flagged1["z_peak"] >= 1.5) & (data_z_flagged1["z_peak"] < 2)]
	chunk4_a = data_fast_flagged1[(data_z_flagged1["z_peak"] >= 2) & (data_z_flagged1["z_peak"] <2.5)]
	
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
	
	chunk2_c = data_fast_flagged2[(data_z_flagged2["z_peak"] >= 1) & (data_z_flagged2["z_peak"] < 1.5)]
	chunk3_c = data_fast_flagged2[(data_z_flagged2["z_peak"] >= 1.5) & (data_z_flagged2["z_peak"] < 2)]
	chunk4_c = data_fast_flagged2[(data_z_flagged2["z_peak"] >= 2) & (data_z_flagged2["z_peak"] <2.5)]
	
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
	
	chunk2_n = data_fast_flagged3[(data_z_flagged3["z_peak"] >= 1) & (data_z_flagged3["z_peak"] < 1.5)]
	chunk3_n = data_fast_flagged3[(data_z_flagged3["z_peak"] >= 1.5) & (data_z_flagged3["z_peak"] < 2)]
	chunk4_n = data_fast_flagged3[(data_z_flagged3["z_peak"] >= 2) & (data_z_flagged3["z_peak"] <2.5)]
	
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
	
	chunk2_s = data_fast_flagged4[(data_z_flagged4["z_peak"] >= 1) & (data_z_flagged4["z_peak"] < 1.5)]
	chunk3_s = data_fast_flagged4[(data_z_flagged4["z_peak"] >= 1.5) & (data_z_flagged4["z_peak"] < 2)]
	chunk4_s = data_fast_flagged4[(data_z_flagged4["z_peak"] >= 2) & (data_z_flagged4["z_peak"] <2.5)]
	
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
	
	chunk2_u = data_fast_flagged5[(data_z_flagged5["z_peak"] >= 1) & (data_z_flagged5["z_peak"] < 1.5)]
	chunk3_u = data_fast_flagged5[(data_z_flagged5["z_peak"] >= 1.5) & (data_z_flagged5["z_peak"] < 2)]
	chunk4_u = data_fast_flagged5[(data_z_flagged5["z_peak"] >= 2) & (data_z_flagged5["z_peak"] <2.5)]
	
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
	
	chunk2_a = data_z_flagged1[(dataz_flagged1["z_peak"] >= 1) & (dataz_flagged1["z_peak"] < 1.5) & (data_fast_flagged1["lmass"] > 11)]
	chunk3_a = data_z_flagged1[(dataz_flagged1["z_peak"] >= 1.5) & (dataz_flagged1["z_peak"] < 2) & (data_fast_flagged1["lmass"] > 11)]
	chunk4_a = data_z_flagged1[(dataz_flagged1["z_peak"] >= 2) & (dataz_flagged1["z_peak"] < 2.5) & (data_fast_flagged1["lmass"] > 11)]
	
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
	
	chunk2_c = data_z_flagged2[(dataz_flagged2["z_peak"] >= 1) & (dataz_flagged2["z_peak"] < 1.5) & (data_fast_flagged2["lmass"] > 11)]
	chunk3_c = data_z_flagged2[(dataz_flagged2["z_peak"] >= 1.5) & (dataz_flagged2["z_peak"] < 2) & (data_fast_flagged2["lmass"] > 11)]
	chunk4_c = data_z_flagged2[(dataz_flagged2["z_peak"] >= 2) & (dataz_flagged2["z_peak"] < 2.5) & (data_fast_flagged2["lmass"] > 11)]
	
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
	
	chunk2_n = data_z_flagged3[(dataz_flagged3["z_peak"] >= 1) & (dataz_flagged3["z_peak"] < 1.5) & (data_fast_flagged3["lmass"] > 11)]
	chunk3_n = data_z_flagged3[(dataz_flagged3["z_peak"] >= 1.5) & (dataz_flagged3["z_peak"] < 2) & (data_fast_flagged3["lmass"] > 11)]
	chunk4_n = data_z_flagged3[(dataz_flagged3["z_peak"] >= 2) & (dataz_flagged3["z_peak"] < 2.5) & (data_fast_flagged3["lmass"] > 11)]
	
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
	
	chunk2_s = data_z_flagged4[(dataz_flagged4["z_peak"] >= 1) & (dataz_flagged4["z_peak"] < 1.5) & (data_fast_flagged4["lmass"] > 11)]
	chunk3_s = data_z_flagged4[(dataz_flagged4["z_peak"] >= 1.5) & (dataz_flagged4["z_peak"] < 2) & (data_fast_flagged4["lmass"] > 11)]
	chunk4_s = data_z_flagged4[(dataz_flagged4["z_peak"] >= 2) & (dataz_flagged4["z_peak"] < 2.5) & (data_fast_flagged4["lmass"] > 11)]
	
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
	
	chunk2_u = data_z_flagged5[(dataz_flagged5["z_peak"] >= 1) & (dataz_flagged5["z_peak"] < 1.5) & (data_fast_flagged5["lmass"] > 11)]
	chunk3_u = data_z_flagged5[(dataz_flagged5["z_peak"] >= 1.5) & (dataz_flagged5["z_peak"] < 2) & (data_fast_flagged5["lmass"] > 11)]
	chunk4_u = data_z_flagged5[(dataz_flagged5["z_peak"] >= 2) & (dataz_flagged5["z_peak"] < 2.5) & (data_fast_flagged5["lmass"] > 11)]
	
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
		
	aegis_1 = [2250, 2481, 6115, 6234, 6424, 6691, 7771, 8838, 11022, 13524, 14012, 14056, 14499, 14609, 16238, 17297, 17570, 17775, 18045, 18257, 19970, 20106, 20250, 20794, 21028, 21357, 22126, 22678, 23011, 23089, 24456, 26850, 26884, 27177, 28328, 29863, 30393, 30421, 30735, 30920, 31326, 32114, 32425, 33028, 33158, 34141, 34254, 34722, 35604, 37919, 38065, 38130, 38167]
	aegis_2 = [195, 766, 1420, 1821, 1925, 2016, 2289, 2427, 3106, 3311, 4826, 5016, 6262, 6310, 6880, 7601, 8646, 10858, 11465, 11730, 12049, 12219, 14481, 14495, 15069, 15088, 16491, 17691, 18280, 18922, 19743, 21156, 22423, 22713, 22741, 22887, 23027, 23045, 23234, 24059, 24333, 24448, 25346, 25969, 26079, 26649, 27315, 28310, 28843, 28864, 29144, 29329, 29399, 31169, 31680, 31839, 33196, 33551, 33770, 34185, 35002, 35677, 36347, 37556, 37592, 38174, 38290, 39028, 39239, 40470]
	aegis_3 = [531, 1606, 2578, 2918, 2957, 8635, 9128, 9870, 10755, 10893, 11416, 11773, 12227, 12479, 15709, 15871, 16065, 17754, 23040, 23645, 25300, 26508, 26952, 27802, 29106, 29178, 29861, 29987, 30967, 31353, 32014, 32686, 33863, 33925, 34685, 34918, 36104, 36574, 37853, 38187]
	cosmos_1 = [796, 2348, 5238, 9111, 10703, 11783, 11871, 12699, 12767, 13206, 13890, 15066, 17263, 18575, 25627, 27769, 31555, 32549]
	cosmos_2 = [312, 363, 728, 2616, 2816, 3200, 4536, 7216, 7411, 9667, 10128, 10592, 10989, 11973, 17089, 17406, 18688, 20983, 21723, 24462, 25534, 25581, 28492, 28523, 29222, 31090]
	cosmos_3 = [490, 1769, 2049, 3182, 3206, 5473, 5530, 6159, 7884, 7951, 9871, 11314, 11337, 11363, 11494, 12020, 12995, 13083, 13174, 16419, 19090, 19153, 20668, 22995, 23021, 23673, 25515, 26039, 26338, 26957, 27289, 28344, 28565, 31922, 33199]
	goodsn_1 = [57, 128, 576, 1749, 2265, 2868, 3133, 4711, 7372, 9056, 10280, 10606, 11706, 11826, 12342, 12561, 13971, 17270, 21156, 23564, 25216, 25813, 32162, 35090, 35299]
	goodsn_2 = [1616, 3186, 4117, 5932, 9692, 10311, 14140, 14532, 16827, 18633, 19913, 20709, 23187, 23548, 25265, 29464, 32842, 33780, 35292, 36582, 37738]
	goodsn_3 = [338, 764, 774, 1678, 2295, 3776, 4854, 4927, 5346, 5371, 5507, 5677, 6215, 6789, 6877, 9122, 10125, 10657, 11064, 12066, 12302, 16129, 16346, 16879, 19082, 20052, 20317, 21738, 23018, 25942, 26529, 26888, 28810, 28826, 30283, 32002, 32033, 36988]
	goodss_1 = [1523, 1924, 2707, 4210, 6098, 6106, 7444, 7503, 19186, 27442, 29928, 30394, 30997, 33163, 33164, 38111, 39170, 43042, 45775, 46392, 46846, 47742, 47873, 48631]
	goodss_2 = [2383, 4505, 7457, 8422, 9704, 10436, 13369, 13628, 14152, 14335, 14747, 15214, 16769, 16814, 24308, 26139, 27881, 29407, 29652, 29900, 31397, 32048, 32783, 33912, 34491, 34519, 34567, 35444, 36095, 39012, 39208, 39364, 40889, 41148, 42113, 42501, 42705, 42957, 43114, 44042, 44157, 48312]
	goodss_3 = [1725, 2467, 4583, 6341, 7686, 8683, 11016, 14813, 15847, 16888, 22079, 22825, 28604, 29288, 30274, 30534, 39568, 40185, 41181, 42607, 43901, 45475, 48464, 49285, 49834]
	uds_1 = [1123, 2393, 2394, 5126, 5924, 6299, 6852, 7071, 7783, 9261, 10758, 11533, 13482, 14152, 14854, 15063, 16239, 17879, 19765, 19954, 21513, 23590, 24953, 26552, 26875, 28773, 30057, 30192, 30255, 32077, 32256, 32691, 32777, 32921, 32986, 34353, 34641, 35071, 35356, 36013, 40631, 41412, 41671, 41835]
	uds_2 = [922, 1513, 1831, 1854, 2294, 3445, 4721, 6590, 6764, 7258, 9073, 10237, 10604, 12441, 12778, 14723, 15270, 18803, 19572, 19703, 19708, 19850, 20529, 20917, 20941, 21031, 21267, 21665, 22480, 25206, 25394, 25630, 27672, 28791, 30133, 30737, 31684, 32468, 32707, 33422, 33527, 35616, 35829, 36010, 36685, 37182, 37775, 38246, 38288, 38631, 39349, 39487, 40420, 40472, 40849, 41302, 41456, 42319, 43367]
	uds_3 = [190, 394, 1620, 2166, 2211, 4059, 4128, 4701, 4706, 5155, 7516, 9207, 9367, 11558, 12010, 13108, 14409, 14996, 15598, 16022, 16709, 17838, 19068, 19126, 20694, 20704, 20770, 21998, 22227, 22416, 22685, 23692, 26581, 28087, 29179, 29461, 30196, 30916, 31615, 32147, 32351, 32947, 34150, 34817, 36247, 38640, 39126, 39624, 42529, 42571, 42812, 43667]
	
	os.chdir("/Volumes/TOSHIBA EXT/3d_hst/noah_massive/aegis_massive")
	for i in aegis_1:
		"data1_a" + str(i) = ascii.read("%s..obs_sed" % (i))
		"lam1_a" + str(i) = ("data1_a" + str(i))["lambda"]
		"flux_1_a" + str(i) = ("data1_a" + (i))["flux_cat"]
		"flux1_a" + str(i) = (("flux_1_a" + str(i))*(("lam1_a" + str(i))**-2.0))*factor
	for i in aegis_2:
		"data2_a" + str(i) = ascii.read("%s..obs_sed" % (i))
		"lam2_a" + str(i) = ("data2_a" + str(i))["lambda"]
		"flux_2_a" + str(i) = ("data2_a" + str(i))["flux_cat"]
		"flux2_a" + str(i) = (("flux_2_a" + str(i))*(("lam2_a" + str(i))**-2.0))*factor
	
	for i in aegis_3:
		
		"data3_a" + str(i) = ascii.read("%s..obs_sed" % (i))
		"lam3_a" + str(i) = ("data3_a" + str(i))["lambda"]
		"flux_3_a" + str(i) = ("data3_a" + (i))["flux_cat"]
		"flux3_a" + str(i) = (("flux_3_a" + str(i))*(("lam3_a" + str(i))**-2.0))*factor
	os.chdir("/Volumes/TOSHIBA EXT/3d_hst/noah_massive/cosmos_massive")	
	for i in cosmos_1:
		"data1_c" + str(i) = ascii.read("%s..obs_sed" % (i))
		"lam1_c" + str(i) = ("data1_c" + str(i))["lambda"]
		"flux_1_c" + str(i) = ("data1_c" + (i))["flux_cat"]
		"flux1_c" + str(i) = (("flux_1_c" + str(i))*(("lam1_c" + str(i))**-2.0))*factor
	for i in cosmos_2:
		"data2_c" + str(i) = ascii.read("%s..obs_sed" % (i))
		"lam2_c" + str(i) = ("data2_c" + str(i))["lambda"]
		"flux_2_c" + str(i) = ("data2_c" + str(i))["flux_cat"]
		"flux2_c" + str(i) = (("flux_2_c" + str(i))*(("lam2_c" + str(i))**-2.0))*factor
	for i in cosmos_3:
		"data3_c" + str(i) = ascii.read("%s..obs_sed" % (i))
		"lam3_c" + str(i) = ("data3_c" + str(i))["lambda"]
		"flux_3_c" + str(i) = ("data3_c" + (i))["flux_cat"]
		"flux3_c" + str(i) = (("flux_3_c" + str(i))*(("lam3_c" + str(i))**-2.0))*factor
	os.chdir("/Volumes/TOSHIBA EXT/3d_hst/noah_massive/goodsn_massive")	
	for i in goodsn_1:
		"data1_n" + str(i) = ascii.read("%s.obs_sed" % (i))
		"lam1_n" + str(i) = ("data1_n" + str(i))["lambda"]
		"flux_1_n" + str(i) = ("data1_n" + (i))["flux_cat"]
		"flux1_n" + str(i) = (("flux_1_n" + str(i))*(("lam1_n" + str(i))**-2.0))*factor
	for i in goodsn_2:
		"data2_n" + str(i) = ascii.read("%s.obs_sed" % (i))
		"lam2_n" + str(i) = ("data2_n" + str(i))["lambda"]
		"flux_2_n" + str(i) = ("data2_n" + str(i))["flux_cat"]
		"flux2_n" + str(i) = (("flux_2_n" + str(i))*(("lam2_n" + str(i))**-2.0))*factor
	for i in goodsn_3:
		"data3_n" + str(i) = ascii.read("%s.obs_sed" % (i))
		"lam3_n" + str(i) = ("data3_n" + str(i))["lambda"]
		"flux_3_n" + str(i) = ("data3_n" + (i))["flux_cat"]
		"flux3_n" + str(i) = (("flux_3_n" + str(i))*(("lam3_n" + str(i))**-2.0))*factor
	os.chdir("/Volumes/TOSHIBA EXT/3d_hst/noah_massive/goodss_massive")	
	for i in goodss_1:
		"data1_s" + str(i) = ascii.read("%s.obs_sed" % (i))
		"lam1_s" + str(i) = ("data1_s" + str(i))["lambda"]
		"flux_1_s" + str(i) = ("data1_s" + (i))["flux_cat"]
		"flux1_s" + str(i) = (("flux_1_s" + str(i))*(("lam1_s" + str(i))**-2.0))*factor
	for i in goodss_2:
		"data2_s" + str(i) = ascii.read("%s.obs_sed" % (i))
		"lam2_s" + str(i) = ("data2_s" + str(i))["lambda"]
		"flux_2_s" + str(i) = ("data2_s" + str(i))["flux_cat"]
		"flux2_s" + str(i) = (("flux_2_s" + str(i))*(("lam2_s" + str(i))**-2.0))*factor
	for i in goodss_3:
		"data3_s" + str(i) = ascii.read("%s.obs_sed" % (i))
		"lam3_s" + str(i) = ("data3_s" + str(i))["lambda"]
		"flux_3_s" + str(i) = ("data3_s" + (i))["flux_cat"]
		"flux3_s" + str(i) = (("flux_3_s" + str(i))*(("lam3_s" + str(i))**-2.0))*factor
	os.chdir("/Volumes/TOSHIBA EXT/3d_hst/noah_massive/uds_massive")	
	for i in uds_1:
		"data1_u" + str(i) = ascii.read("%s.obs_sed" % (i))
		"lam1_u" + str(i) = ("data1_u" + str(i))["lambda"]
		"flux_1_u" + str(i) = ("data1_u" + (i))["flux_cat"]
		"flux1_u" + str(i) = (("flux_1_u" + str(i))*(("lam1_u" + str(i))**-2.0))*factor
	for i in uds_2:
		"data2_u" + str(i) = ascii.read("%s.obs_sed" % (i))
		"lam2_u" + str(i) = ("data2_u" + str(i))["lambda"]
		"flux_2_u" + str(i) = ("data2_u" + str(i))["flux_cat"]
		"flux2_u" + str(i) = (("flux_2_u" + str(i))*(("lam2_u" + str(i))**-2.0))*factor
	for i in uds_3:
		"data3_u" + str(i) = ascii.read("%s.obs_sed" % (i))
		"lam3_u" + str(i) = ("data3_u" + str(i))["lambda"]
		"flux_3_u" + str(i) = ("data3_u" + (i))["flux_cat"]
		"flux3_u" + str(i) = (("flux_3_u" + str(i))*(("lam3_u" + str(i))**-2.0))*factor
		
	
	


	
	
	
	
	
	
	os.chdir("/Volumes/TOSHIBA EXT/3d_hst")
	
	
	
	dataz_a = ascii.read("aegis_3dhst.v4.1.zout")
	dataz_c = ascii.read("cosmos_3dhst.v4.1.zout")
	dataz_n = ascii.read("goodsn_3dhst.v4.1.zout")
	dataz_s = ascii.read("goodss_3dhst.v4.1.zout")
	dataz_u = ascii.read("uds_3dhst.v4.1.zout")
	
	for i in aegis_1:
		"dataz1_a" + str(i) = dataz_a[(dataz_a["id"] == i)]
	for i in aegis_2:
		"dataz2_a" + str(i) = dataz_a[(dataz_a["id"] == i)]
	for i in aegis_3:
		"dataz3_a" + str(i) = dataz_a[(dataz_a["id"] == i)]
	for i in cosmos_1:
		"dataz1_c" + str(i) = dataz_c[(dataz_c["id"] == i)]
	for i in cosmos_2:
		"dataz2_c" + str(i) = dataz_c[(dataz_c["id"] == i)]
	for i in cosmos_3:
		"dataz3_c" + str(i) = dataz_c[(dataz_c["id"] == i)]
	for i in goodsn_1:
		"dataz1_n" + str(i) = dataz_n[(dataz_n["id"] == i)]
	for i in goodsn_2:
		"dataz2_n" + str(i) = dataz_n[(dataz_n["id"] == i)]
	for i in goodsn_3:
		"dataz3_n" + str(i) = dataz_n[(dataz_n["id"] == i)]
	for i in goodss_1:
		"dataz1_s" + str(i) = dataz_s[(dataz_s["id"] == i)]
	for i in goodss_2:
		"dataz2_s" + str(i) = dataz_s[(dataz_s["id"] == i)]
	for i in goodss_3:
		"dataz3_s" + str(i) = dataz_s[(dataz_s["id"] == i)]
	for i in uds_1:
		"dataz1_u" + str(i) = dataz_u[(dataz_u["id"] == i)]
	for i in uds_2:
		"dataz2_u" + str(i) = dataz_u[(dataz_u["id"] == i)]
	for i in uds_3:
		"dataz3_u" + str(i) = dataz_u[(dataz_u["id"] == i)]
	
	
	for i in aegis_1:
		"lam_1_a" + str(i) = ("lam1_a" + str(i))/(1+ ("dataz1_a" + str(i))["z_peak"])
	for i in aegis_2:
		"lam_2_a" + str(i) = ("lam2_a" + str(i))/(1+ ("dataz2_a" + str(i))["z_peak"])
	for i in aegis_3:
		"lam_3_a" + str(i) = ("lam3_a" + str(i))/(1+ ("dataz3_a" + str(i))["z_peak"])
	for i in cosmos_1:
		"lam_1_c" + str(i) = ("lam1_c" + str(i))/(1+ ("dataz1_c" + str(i))["z_peak"])
	for i in cosmos_2:
		"lam_2_c" + str(i) = ("lam2_c" + str(i))/(1+ ("dataz2_c" + str(i))["z_peak"])
	for i in cosmos_3:
		"lam_3_c" + str(i) = ("lam3_c" + str(i))/(1+ ("dataz3_c" + str(i))["z_peak"])
	for i in goodsn_1:
		"lam_1_n" + str(i) = ("lam1_n" + str(i))/(1+ ("dataz1_n" + str(i))["z_peak"])
	for i in goodsn_2:
		"lam_2_n" + str(i) = ("lam2_n" + str(i))/(1+ ("dataz2_n" + str(i))["z_peak"])
	for i in goodsn_3:
		"lam_3_n" + str(i) = ("lam3_n" + str(i))/(1+ ("dataz3_n" + str(i))["z_peak"])
	for i in goodss_1:
		"lam_1_s" + str(i) = ("lam1_s" + str(i))/(1+ ("dataz1_s" + str(i))["z_peak"])
	for i in goodss_2:
		"lam_2_s" + str(i) = ("lam2_s" + str(i))/(1+ ("dataz2_s" + str(i))["z_peak"])
	for i in goodss_3:
		"lam_3_s" + str(i) = ("lam3_s" + str(i))/(1+ ("dataz3_s" + str(i))["z_peak"])
	for i in uds_1:
		"lam_1_u" + str(i) = ("lam1_u" + str(i))/(1+ ("dataz1_u" + str(i))["z_peak"])
	for i in uds_2:
		"lam_2_u" + str(i) = ("lam2_u" + str(i))/(1+ ("dataz2_u" + str(i))["z_peak"])
	for i in uds_3:
		"lam_3_u" + str(i) = ("lam3_u" + str(i))/(1+ ("dataz3_u" + str(i))["z_peak"])
	
	
	
	
	
	
	fig,pylab.axes = pylab.subplots(3, 1, sharex=True)

	a1 = pylab.axes[0]
	a2 = pylab.axes[1]
	a3 = pylab.axes[2]
	
	
	for i in aegis_1:
		a1.plot(("lam_1_a" + str(i)), ("flux1_a" + str(i)), color="g", alpha=0.7, markeredgecolor="none", linestyle="none", marker="o")
	for i in aegis_2:
		a2.plot(("lam_2_a" + str(i)), ("flux2_a" + str(i)), color="b", alpha=0.7, markeredgecolor="none", linestyle="none", marker="o")
	for i in aegis_3:
		a3.plot(("lam_3_a" + str(i)), ("flux3_a" + str(i)), color="purple", alpha=0.7, markeredgecolor="none", linestyle="none", marker="o")
	
	for i in cosmos_1:
		a1.plot(("lam_1_c" + str(i)), ("flux1_c" + str(i)), color="g", alpha=0.7, markeredgecolor="none", linestyle="none", marker="o")
	for i in cosmos_2:
		a2.plot(("lam_2_c" + str(i)), ("flux2_c" + str(i)), color="b", alpha=0.7, markeredgecolor="none", linestyle="none", marker="o")
	for i in cosmos_3:
		a3.plot(("lam_3_c" + str(i)), ("flux3_c" + str(i)), color="purple", alpha=0.7, markeredgecolor="none", linestyle="none", marker="o")
	
	for i in goodsn_1:
		a1.plot(("lam_1_n" + str(i)), ("flux1_n" + str(i)), color="g", alpha=0.7, markeredgecolor="none", linestyle="none", marker="o")
	for i in goodsn_2:
		a2.plot(("lam_2_n" + str(i)), ("flux2_n" + str(i)), color="b", alpha=0.7, markeredgecolor="none", linestyle="none", marker="o")
	for i in goodsn_3:
		a3.plot(("lam_3_n" + str(i)), ("flux3_n" + str(i)), color="purple", alpha=0.7, markeredgecolor="none", linestyle="none", marker="o")
	
	for i in goodss_1:
		a1.plot(("lam_1_s" + str(i)), ("flux1_s" + str(i)), color="g", alpha=0.7, markeredgecolor="none", linestyle="none", marker="o")
	for i in goodss_2:
		a2.plot(("lam_2_s" + str(i)), ("flux2_s" + str(i)), color="b", alpha=0.7, markeredgecolor="none", linestyle="none", marker="o")
	for i in goodss_3:
		a3.plot(("lam_3_s" + str(i)), ("flux3_s" + str(i)), color="purple", alpha=0.7, markeredgecolor="none", linestyle="none", marker="o")
	
	for i in uds_1:
		a1.plot(("lam_1_u" + str(i)), ("flux1_u" + str(i)), color="g", alpha=0.7, markeredgecolor="none", linestyle="none", marker="o")
	for i in uds_2:
		a2.plot(("lam_2_u" + str(i)), ("flux2_u" + str(i)), color="b", alpha=0.7, markeredgecolor="none", linestyle="none", marker="o")
	for i in uds_3:
		a3.plot(("lam_3_u" + str(i)), ("flux3_u" + str(i)), color="purple", alpha=0.7, markeredgecolor="none", linestyle="none", marker="o")
	
	
	
	
	os.chdir("/Volumes/TOSHIBA EXT/3d_hst/noah_massive/aegis_massive")
	
	for i in aegis_1:
		"_data1_a" + str(i) = ascii.read("%s..temp_sed" % (i))
		"_flux_1_a" + str(i) = ("_data1_a" + str(i))["tempflux"]
		"_lam1_a" + str(i) = ("_data1_a" + str(i))["lambda"]
		"_lam_1_a" + str(i) = ("_lam1_a" + str(i))/(1+ ("dataz1_a" + str(i))["z_peak"])
		"_flux1_a" + str(i) = (("_flux_1_a" + str(i))*(("_lam_1_a" + str(i))**-2.0))*factor
		
	for i in aegis_2:
		"_data2_a" + str(i) = ascii.read("%s..temp_sed" % (i))
		"_flux_2_a" + str(i) = ("_data2_a" + str(i))["tempflux"]
		"_lam2_a" + str(i) = ("_data2_a" + str(i))["lambda"]
		"_lam_2_a" + str(i) = ("_lam2_a" + str(i))/(1+ ("dataz2_a" + str(i))["z_peak"])
		"_flux2_a" + str(i) = (("_flux_2_a" + str(i))*(("_lam_2_a" + str(i))**-2.0))*factor
		
		
	for i in aegis_3:
		"_data3_a" + str(i) = ascii.read("%s..temp_sed" % (i))
		"_flux_3_a" + str(i) = ("_data3_a" + str(i))["tempflux"]
		"_lam3_a" + str(i) = ("_data3_a" + str(i))["lambda"]
		"_lam_3_a" + str(i) = ("_lam3_a" + str(i))/(1+ ("dataz3_a" + str(i))["z_peak"])
		"_flux3_a" + str(i) = (("_flux_3_a" + str(i))*(("_lam_3_a" + str(i))**-2.0))*factor
		
		
	os.chdir("/Volumes/TOSHIBA EXT/3d_hst/noah_massive/cosmos_massive")	
		
	for i in cosmos_1:
		"_data1_c" + str(i) = ascii.read("%s..temp_sed" % (i))
		"_flux_1_c" + str(i) = ("_data1_c" + str(i))["tempflux"]
		"_lam1_c" + str(i) = ("_data1_c" + str(i))["lambda"]
		"_lam_1_c" + str(i) = ("_lam1_c" + str(i))/(1+ ("dataz1_c" + str(i))["z_peak"])
		"_flux1_c" + str(i) = (("_flux_1_c" + str(i))*(("_lam_1_c" + str(i))**-2.0))*factor
		
	for i in cosmos_2:
		"_data2_c" + str(i) = ascii.read("%s..temp_sed" % (i))
		"_flux_2_c" + str(i) = ("_data2_c" + str(i))["tempflux"]
		"_lam2_c" + str(i) = ("_data2_c" + str(i))["lambda"]
		"_lam_2_c" + str(i) = ("_lam2_c" + str(i))/(1+ ("dataz2_c" + str(i))["z_peak"])
		"_flux2_c" + str(i) = (("_flux_2_c" + str(i))*(("_lam_2_c" + str(i))**-2.0))*factor
		
		
	for i in cosmos_3:
		"_data3_c" + str(i) = ascii.read("%s..temp_sed" % (i))
		"_flux_3_c" + str(i) = ("_data3_c" + str(i))["tempflux"]
		"_lam3_c" + str(i) = ("_data3_c" + str(i))["lambda"]
		"_lam_3_c" + str(i) = ("_lam3_c" + str(i))/(1+ ("dataz3_c" + str(i))["z_peak"])
		"_flux3_c" + str(i) = (("_flux_3_c" + str(i))*(("_lam_3_c" + str(i))**-2.0))*factor
		
	
	os.chdir("/Volumes/TOSHIBA EXT/3d_hst/noah_massive/goodsn_massive")
		
	for i in goodsn_1:
		"_data1_n" + str(i) = ascii.read("%s.temp_sed" % (i))
		"_flux_1_n" + str(i) = ("_data1_n" + str(i))["tempflux"]
		"_lam1_n" + str(i) = ("_data1_n" + str(i))["lambda"]
		"_lam_1_n" + str(i) = ("_lam1_n" + str(i))/(1+ ("dataz1_n" + str(i))["z_peak"])
		"_flux1_n" + str(i) = (("_flux_1_n" + str(i))*(("_lam_1_n" + str(i))**-2.0))*factor
		
	for i in goodsn_2:
		"_data2_n" + str(i) = ascii.read("%s.temp_sed" % (i))
		"_flux_2_n" + str(i) = ("_data2_n" + str(i))["tempflux"]
		"_lam2_n" + str(i) = ("_data2_n" + str(i))["lambda"]
		"_lam_2_n" + str(i) = ("_lam2_n" + str(i))/(1+ ("dataz2_n" + str(i))["z_peak"])
		"_flux2_n" + str(i) = (("_flux_2_n" + str(i))*(("_lam_2_n" + str(i))**-2.0))*factor
		
		
	for i in goodsn_3:
		"_data3_n" + str(i) = ascii.read("%s.temp_sed" % (i))
		"_flux_3_n" + str(i) = ("_data3_n" + str(i))["tempflux"]
		"_lam3_n" + str(i) = ("_data3_n" + str(i))["lambda"]
		"_lam_3_n" + str(i) = ("_lam3_n" + str(i))/(1+ ("dataz3_n" + str(i))["z_peak"])
		"_flux3_n" + str(i) = (("_flux_3_n" + str(i))*(("_lam_3_n" + str(i))**-2.0))*factor
		
	
	os.chdir("/Volumes/TOSHIBA EXT/3d_hst/noah_massive/goodss_massive")
		
	for i in goodss_1:
		"_data1_s" + str(i) = ascii.read("%s.temp_sed" % (i))
		"_flux_1_s" + str(i) = ("_data1_s" + str(i))["tempflux"]
		"_lam1_s" + str(i) = ("_data1_s" + str(i))["lambda"]
		"_lam_1_s" + str(i) = ("_lam1_s" + str(i))/(1+ ("dataz1_s" + str(i))["z_peak"])
		"_flux1_s" + str(i) = (("_flux_1_s" + str(i))*(("_lam_1_s" + str(i))**-2.0))*factor
		
	for i in goodss_2:
		"_data2_s" + str(i) = ascii.read("%s.temp_sed" % (i))
		"_flux_2_s" + str(i) = ("_data2_s" + str(i))["tempflux"]
		"_lam2_s" + str(i) = ("_data2_s" + str(i))["lambda"]
		"_lam_2_s" + str(i) = ("_lam2_s" + str(i))/(1+ ("dataz2_s" + str(i))["z_peak"])
		"_flux2_s" + str(i) = (("_flux_2_s" + str(i))*(("_lam_2_s" + str(i))**-2.0))*factor
		
		
	for i in goodss_3:
		"_data3_s" + str(i) = ascii.read("%s.temp_sed" % (i))
		"_flux_3_s" + str(i) = ("_data3_s" + str(i))["tempflux"]
		"_lam3_s" + str(i) = ("_data3_s" + str(i))["lambda"]
		"_lam_3_s" + str(i) = ("_lam3_s" + str(i))/(1+ ("dataz3_s" + str(i))["z_peak"])
		"_flux3_s" + str(i) = (("_flux_3_s" + str(i))*(("_lam_3_s" + str(i))**-2.0))*factor
		
	
	os.chdir("/Volumes/TOSHIBA EXT/3d_hst/noah_massive/uds_massive")
		
	for i in uds_1:
		"_data1_u" + str(i) = ascii.read("%s.temp_sed" % (i))
		"_flux_1_u" + str(i) = ("_data1_u" + str(i))["tempflux"]
		"_lam1_u" + str(i) = ("_data1_u" + str(i))["lambda"]
		"_lam_1_u" + str(i) = ("_lam1_u" + str(i))/(1+ ("dataz1_u" + str(i))["z_peak"])
		"_flux1_u" + str(i) = (("_flux_1_u" + str(i))*(("_lam_1_u" + str(i))**-2.0))*factor
		
	for i in uds_2:
		"_data2_u" + str(i) = ascii.read("%s.temp_sed" % (i))
		"_flux_2_u" + str(i) = ("_data2_u" + str(i))["tempflux"]
		"_lam2_u" + str(i) = ("_data2_u" + str(i))["lambda"]
		"_lam_2_u" + str(i) = ("_lam2_u" + str(i))/(1+ ("dataz2_u" + str(i))["z_peak"])
		"_flux2_u" + str(i) = (("_flux_2_u" + str(i))*(("_lam_2_u" + str(i))**-2.0))*factor
		
		
	for i in uds_3:
		"_data3_u" + str(i) = ascii.read("%s.temp_sed" % (i))
		"_flux_3_u" + str(i) = ("_data3_u" + str(i))["tempflux"]
		"_lam3_u" + str(i) = ("_data3_u" + str(i))["lambda"]
		"_lam_3_u" + str(i) = ("_lam3_u" + str(i))/(1+ ("dataz3_u" + str(i))["z_peak"])
		"_flux3_u" + str(i) = (("_flux_3_u" + str(i))*(("_lam_3_u" + str(i))**-2.0))*factor
		
	
	chunk1 = zip(_flux1_a2250, _flux1_a2481, _flux1_a6115, _flux1_a6234, _flux1_a6424, _flux1_a6691, _flux1_a7771, _flux1_a8838, _flux1_a11022, _flux1_a13524, _flux1_a14012, _flux1_a14056, _flux1_a14499, _flux1_a14609, _flux1_a16238, _flux1_a17297, _flux1_a17570, _flux1_a17775, _flux1_a18045, _flux1_a18257, _flux1_a19970, _flux1_a20106, _flux1_a20250, _flux1_a20794, _flux1_a21028, _flux1_a21357, _flux1_a22126, _flux1_a22678, _flux1_a23011, _flux1_a23089, _flux1_a24456, _flux1_a26850, _flux1_a26884, _flux1_a27177, _flux1_a28328, _flux1_a29863, _flux1_a30393, _flux1_a30421, _flux1_a30735, _flux1_a30920, _flux1_a31326, _flux1_a32114, _flux1_a32425, _flux1_a33028, _flux1_a33158, _flux1_a34141, _flux1_a34254, _flux1_a34722, _flux1_a35604, _flux1_a37919, _flux1_a38065, _flux1_a38130, _flux1_a38167, _flux1_c796, _flux1_c2348, _flux1_c5238, _flux1_c9111, _flux1_c10703, _flux1_c11783, _flux1_c11871, _flux1_c12699, _flux1_c12767, _flux1_c13206, _flux1_c13890, _flux1_c15066, _flux1_c17263, _flux1_c18575, _flux1_c25627, _flux1_c27769, _flux1_c31555, _flux1_c32549, _flux1_n57, _flux1_n128, _flux1_n576, _flux1_n1749, _flux1_n2265, _flux1_n2868, _flux1_n3133, _flux1_n4711, _flux1_n7372, _flux1_n9056, _flux1_n10280, _flux1_n10606, _flux1_n11706, _flux1_n11826, _flux1_n12342, _flux1_n12561, _flux1_n13971, _flux1_n17270, _flux1_n21156, _flux1_n23564, _flux1_n25216, _flux1_n25813, _flux1_n32162, _flux1_n35090, _flux1_n35299, _flux1_s1523, _flux1_s1924, _flux1_s2707, _flux1_s4210, _flux1_s6098, _flux1_s6106, _flux1_s7444, _flux1_s7503, _flux1_s19186, _flux1_s27442, _flux1_s29928, _flux1_s30394, _flux1_s30997, _flux1_s33163, _flux1_s33164, _flux1_s38111, _flux1_s39170, _flux1_s43042, _flux1_s45775, _flux1_s46392, _flux1_s46846, _flux1_s47742, _flux1_s47873, _flux1_s48631, _flux1_u1123, _flux1_u2393, _flux1_u2394, _flux1_u5126, _flux1_u5924, _flux1_u6299, _flux1_u6852, _flux1_u7071, _flux1_u7783, _flux1_u9261, _flux1_u10758, _flux1_u11533, _flux1_u13482, _flux1_u14152, _flux1_u14854, _flux1_u15063, _flux1_u16239, _flux1_u17879, _flux1_u19765, _flux1_u19954, _flux1_u21513, _flux1_u23590, _flux1_u24953, _flux1_u26552, _flux1_u26875, _flux1_u28773, _flux1_u30057, _flux1_u30192, _flux1_u30255, _flux1_u32077, _flux1_u32256, _flux1_u32691, _flux1_u32777, _flux1_u32921, _flux1_u32986, _flux1_u34353, _flux1_u34641, _flux1_u35071, _flux1_u35356, _flux1_u36013, _flux1_u40631, _flux1_u41412, _flux1_u41671, _flux1_u41835)
	chunk2 = zip(_flux2_a195, _flux2_a766, _flux2_a1420, _flux2_a1821, _flux2_a1925, _flux2_a2016, _flux2_a2289, _flux2_a2427, _flux2_a3106, _flux2_a3311, _flux2_a4826, _flux2_a5016, _flux2_a6262, _flux2_a6310, _flux2_a6880, _flux2_a7601, _flux2_a8646, _flux2_a10858, _flux2_a11465, _flux2_a11730, _flux2_a12049, _flux2_a12219, _flux2_a14481, _flux2_a14495, _flux2_a15069, _flux2_a15088, _flux2_a16491, _flux2_a17691, _flux2_a18280, _flux2_a18922, _flux2_a19743, _flux2_a21156, _flux2_a22423, _flux2_a22713, _flux2_a22741, _flux2_a22887, _flux2_a23027, _flux2_a23045, _flux2_a23234, _flux2_a24059, _flux2_a24333, _flux2_a24448, _flux2_a25346, _flux2_a25969, _flux2_a26079, _flux2_a26649, _flux2_a27315, _flux2_a28310, _flux2_a28843, _flux2_a28864, _flux2_a29144, _flux2_a29329, _flux2_a29399, _flux2_a31169, _flux2_a31680, _flux2_a31839, _flux2_a33196, _flux2_a33551, _flux2_a33770, _flux2_a34185, _flux2_a35002, _flux2_a35677, _flux2_a36347, _flux2_a37556, _flux2_a37592, _flux2_a38174, _flux2_a38290, _flux2_a39028, _flux2_a39239, _flux2_a40470, _flux2_c312, _flux2_c363, _flux2_c728, _flux2_c2616, _flux2_c2816, _flux2_c3200, _flux2_c4536, _flux2_c7216, _flux2_c7411, _flux2_c9667, _flux2_c10128, _flux2_c10592, _flux2_c10989, _flux2_c11973, _flux2_c17089, _flux2_c17406, _flux2_c18688, _flux2_c20983, _flux2_c21723, _flux2_c24462, _flux2_c25534, _flux2_c25581, _flux2_c28492, _flux2_c28523, _flux2_c29222, _flux2_c31090, _flux2_n1616, _flux2_n3186, _flux2_n4117, _flux2_n5932, _flux2_n9692, _flux2_n10311, _flux2_n14140, _flux2_n14532, _flux2_n16827, _flux2_n18633, _flux2_n19913, _flux2_n20709, _flux2_n23187, _flux2_n23548, _flux2_n25265, _flux2_n29464, _flux2_n32842, _flux2_n33780, _flux2_n35292, _flux2_n36582, _flux2_n37738, _flux2_s2383, _flux2_s4505, _flux2_s7457, _flux2_s8422, _flux2_s9704, _flux2_s10436, _flux2_s13369, _flux2_s13628, _flux2_s14152, _flux2_s14335, _flux2_s14747, _flux2_s15214, _flux2_s16769, _flux2_s16814, _flux2_s24308, _flux2_s26139, _flux2_s27881, _flux2_s29407, _flux2_s29652, _flux2_s29900, _flux2_s31397, _flux2_s32048, _flux2_s32783, _flux2_s33912, _flux2_s34491, _flux2_s34519, _flux2_s34567, _flux2_s35444, _flux2_s36095, _flux2_s39012, _flux2_s39208, _flux2_s39364, _flux2_s40889, _flux2_s41148, _flux2_s42113, _flux2_s42501, _flux2_s42705, _flux2_s42957, _flux2_s43114, _flux2_s44042, _flux2_s44157, _flux2_s48312, _flux2_u922, _flux2_u1513, _flux2_u1831, _flux2_u1854, _flux2_u2294, _flux2_u3445, _flux2_u4721, _flux2_u6590, _flux2_u6764, _flux2_u7258, _flux2_u9073, _flux2_u10237, _flux2_u10604, _flux2_u12441, _flux2_u12778, _flux2_u14723, _flux2_u15270, _flux2_u18803, _flux2_u19572, _flux2_u19703, _flux2_u19708, _flux2_u19850, _flux2_u20529, _flux2_u20917, _flux2_u20941, _flux2_u21031, _flux2_u21267, _flux2_u21665, _flux2_u22480, _flux2_u25206, _flux2_u25394, _flux2_u25630, _flux2_u27672, _flux2_u28791, _flux2_u30133, _flux2_u30737, _flux2_u31684, _flux2_u32468, _flux2_u32707, _flux2_u33422, _flux2_u33527, _flux2_u35616, _flux2_u35829, _flux2_u36010, _flux2_u36685, _flux2_u37182, _flux2_u37775, _flux2_u38246, _flux2_u38288, _flux2_u38631, _flux2_u39349, _flux2_u39487, _flux2_u40420, _flux2_u40472, _flux2_u40849, _flux2_u41302, _flux2_u41456, _flux2_u42319, _flux2_u43367)
	chunk3 = zip(_flux3_a531, _flux3_a1606, _flux3_a2578, _flux3_a2918, _flux3_a2957, _flux3_a8635, _flux3_a9128, _flux3_a9870, _flux3_a10755, _flux3_a10893, _flux3_a11416, _flux3_a11773, _flux3_a12227, _flux3_a12479, _flux3_a15709, _flux3_a15871, _flux3_a16065, _flux3_a17754, _flux3_a23040, _flux3_a23645, _flux3_a25300, _flux3_a26508, _flux3_a26952, _flux3_a27802, _flux3_a29106, _flux3_a29178, _flux3_a29861, _flux3_a29987, _flux3_a30967, _flux3_a31353, _flux3_a32014, _flux3_a32686, _flux3_a33863, _flux3_a33925, _flux3_a34685, _flux3_a34918, _flux3_a36104, _flux3_a36574, _flux3_a37853, _flux3_a38187, _flux3_c490, _flux3_c1769, _flux2_c2049, _flux2_c3182, _flux3_c3206, _flux3_c5473, _flux3_c5530, _flux3_c6159, _flux3_c7884, _flux3_c7951, _flux3_c9871, _flux3_c11314, _flux3_c11337, _flux3_c11363, _flux3_c11494, _flux3_c12020, _flux3_c12995, _flux3_c13083, _flux3_c13174, _flux3_c16419, _flux3_c19090, _flux3_c19153, _flux3_c20668, _flux3_c22995, _flux3_c23021, _flux3_c23673, _flux3_c25515, _flux3_c26039, _flux3_c26338, _flux3_c26957, _flux3_c27289, _flux3_c28344, _flux3_c28565, _flux3_c31922, _flux3_c33199, _flux3_n338, _flux3_n764, _flux2_n774, _flux2_n1678, _flux3_n2295, _flux3_n3776, _flux3_n4854, _flux3_n4927, _flux3_n5346, _flux3_n5371, _flux3_n5507, _flux3_n5677, _flux3_n6215, _flux3_n6789, _flux3_n6877, _flux3_n9122, _flux3_n10125, _flux3_n10657, _flux3_n11064, _flux3_n12066, _flux3_n12302, _flux3_n16129, _flux3_n16346, _flux3_n16879, _flux3_n19082, _flux3_n20052, _flux3_n20317, _flux3_n21738, _flux3_n23018, _flux3_n25942, _flux3_n26529, _flux3_n26888, _flux3_n28810, _flux3_n28826, _flux3_n30283, _flux3_n32002, _flux3_n32033, _flux3_n36988, _flux3_s1725, _flux3_s2467, _flux2_s4583, _flux2_s6341, _flux3_s7686, _flux3_s8683, _flux3_s11016, _flux3_s14813, _flux3_s15847, _flux3_s16888, _flux3_s22079, _flux3_s22825, _flux3_s28604, _flux3_s29288, _flux3_s30274, _flux3_s30534, _flux3_s39568, _flux3_s40185, _flux3_s41181, _flux3_s42607, _flux3_s43901, _flux3_s45475, _flux3_s48464, _flux3_s49285, _flux3_s49834, _flux3_u190, _flux3_u394, _flux2_u1620, _flux2_u2166, _flux3_u2211, _flux3_u4059, _flux3_u4128, _flux3_u4701, _flux3_u4706, _flux3_u5155, _flux3_u7516, _flux3_u9207, _flux3_u9367, _flux3_u11558, _flux3_u12010, _flux3_u13108, _flux3_u14409, _flux3_u14996, _flux3_u15598, _flux3_u16022, _flux3_u16709, _flux3_u17838, _flux3_u19068, _flux3_u19126, _flux3_u20694, _flux3_u20704, _flux3_u20770, _flux3_u21998, _flux3_u22227, _flux3_u22416, _flux3_u22685, _flux3_u23692, _flux3_u26581, _flux3_u28087, _flux3_u29179, _flux3_u29461, _flux3_u30196, _flux3_u30916, _flux3_u31615, _flux3_u32147, _flux3_u32351, _flux3_u32947, _flux3_u34150, _flux3_u34817, _flux3_u36247, _flux3_u38640, _flux3_u39126, _flux3_u39624, _flux3_u42529, _flux3_u42571, _flux3_u42812, _flux3_u43667)
	
	chunked1 = np.median(chunk1, axis=1)
	chunked2 = np.median(chunk2, axis=1)
	chunked3 = np.median(chunk3, axis=1)
	
	a1.plot(_lam_1_a1, chunked1, color="black", label="median", lw=2)
	a2.plot(_lam_2_a1, chunked2, color="black", label="median", lw=2)
	a3.plot(_lam_3_a1, chunked3, color="black", label="median", lw=2)
	
	
	
	
	
	
	
	a1.set_xlim([800,40000])
	a1.set_ylim([0,2.4])
	a1.set_xscale("log")
	a2.set_xlim([800,40000])
	a2.set_ylim([0,2.4])
	a2.set_xscale("log")
	a3.set_xlim([800,40000])
	a3.set_ylim([0,2.4])
	a3.set_xscale("log")
	pylab.suptitle("Wavelength vs Flux", fontsize=19)
	fig.text(0.5, 0.01, "Wavelength", fontsize=18)
	fig.text(0.01, 0.5, "Flux", rotation = "vertical", fontsize=18)
	fig.subplots_adjust(hspace=0, wspace=0)
	
	os.chdir("/Volumes/TOSHIBA EXT/3d_hst")
	
	pylab.ion()
	pylab.show()
		
def test():
	aegis_1 = [2250, 2481, 6115, 6234, 6424, 6691, 7771, 8838, 11022, 13524, 14012, 14056, 14499, 14609, 16238, 17297, 17570, 17775, 18045, 18257, 19970, 20106, 20250, 20794, 21028, 21357, 22126, 22678, 23011, 23089, 24456, 26850, 26884, 27177, 28328, 29863, 30393, 30421, 30735, 30920, 31326, 32114, 32425, 33028, 33158, 34141, 34254, 34722, 35604, 37919, 38065, 38130, 38167]
	aegis_2 = [195, 766, 1420, 1821, 1925, 2016, 2289, 2427, 3106, 3311, 4826, 5016, 6262, 6310, 6880, 7601, 8646, 10858, 11465, 11730, 12049, 12219, 14481, 14495, 15069, 15088, 16491, 17691, 18280, 18922, 19743, 21156, 22423, 22713, 22741, 22887, 23027, 23045, 23234, 24059, 24333, 24448, 25346, 25969, 26079, 26649, 27315, 28310, 28843, 28864, 29144, 29329, 29399, 31169, 31680, 31839, 33196, 33551, 33770, 34185, 35002, 35677, 36347, 37556, 37592, 38174, 38290, 39028, 39239, 40470]
	aegis_3 = [531, 1606, 2578, 2918, 2957, 8635, 9128, 9870, 10755, 10893, 11416, 11773, 12227, 12479, 15709, 15871, 16065, 17754, 23040, 23645, 25300, 26508, 26952, 27802, 29106, 29178, 29861, 29987, 30967, 31353, 32014, 32686, 33863, 33925, 34685, 34918, 36104, 36574, 37853, 38187]
	cosmos_1 = [796, 2348, 5238, 9111, 10703, 11783, 11871, 12699, 12767, 13206, 13890, 15066, 17263, 18575, 25627, 27769, 31555, 32549]
	cosmos_2 = [312, 363, 728, 2616, 2816, 3200, 4536, 7216, 7411, 9667, 10128, 10592, 10989, 11973, 17089, 17406, 18688, 20983, 21723, 24462, 25534, 25581, 28492, 28523, 29222, 31090]
	cosmos_3 = [490, 1769, 2049, 3182, 3206, 5473, 5530, 6159, 7884, 7951, 9871, 11314, 11337, 11363, 11494, 12020, 12995, 13083, 13174, 16419, 19090, 19153, 20668, 22995, 23021, 23673, 25515, 26039, 26338, 26957, 27289, 28344, 28565, 31922, 33199]
	goodsn_1 = [57, 128, 576, 1749, 2265, 2868, 3133, 4711, 7372, 9056, 10280, 10606, 11706, 11826, 12342, 12561, 13971, 17270, 21156, 23564, 25216, 25813, 32162, 35090, 35299]
	goodsn_2 = [1616, 3186, 4117, 5932, 9692, 10311, 14140, 14532, 16827, 18633, 19913, 20709, 23187, 23548, 25265, 29464, 32842, 33780, 35292, 36582, 37738]
	goodsn_3 = [338, 764, 774, 1678, 2295, 3776, 4854, 4927, 5346, 5371, 5507, 5677, 6215, 6789, 6877, 9122, 10125, 10657, 11064, 12066, 12302, 16129, 16346, 16879, 19082, 20052, 20317, 21738, 23018, 25942, 26529, 26888, 28810, 28826, 30283, 32002, 32033, 36988]
	goodss_1 = [1523, 1924, 2707, 4210, 6098, 6106, 7444, 7503, 19186, 27442, 29928, 30394, 30997, 33163, 33164, 38111, 39170, 43042, 45775, 46392, 46846, 47742, 47873, 48631]
	goodss_2 = [2383, 4505, 7457, 8422, 9704, 10436, 13369, 13628, 14152, 14335, 14747, 15214, 16769, 16814, 24308, 26139, 27881, 29407, 29652, 29900, 31397, 32048, 32783, 33912, 34491, 34519, 34567, 35444, 36095, 39012, 39208, 39364, 40889, 41148, 42113, 42501, 42705, 42957, 43114, 44042, 44157, 48312]
	goodss_3 = [1725, 2467, 4583, 6341, 7686, 8683, 11016, 14813, 15847, 16888, 22079, 22825, 28604, 29288, 30274, 30534, 39568, 40185, 41181, 42607, 43901, 45475, 48464, 49285, 49834]
	uds_1 = [1123, 2393, 2394, 5126, 5924, 6299, 6852, 7071, 7783, 9261, 10758, 11533, 13482, 14152, 14854, 15063, 16239, 17879, 19765, 19954, 21513, 23590, 24953, 26552, 26875, 28773, 30057, 30192, 30255, 32077, 32256, 32691, 32777, 32921, 32986, 34353, 34641, 35071, 35356, 36013, 40631, 41412, 41671, 41835]
	uds_2 = [922, 1513, 1831, 1854, 2294, 3445, 4721, 6590, 6764, 7258, 9073, 10237, 10604, 12441, 12778, 14723, 15270, 18803, 19572, 19703, 19708, 19850, 20529, 20917, 20941, 21031, 21267, 21665, 22480, 25206, 25394, 25630, 27672, 28791, 30133, 30737, 31684, 32468, 32707, 33422, 33527, 35616, 35829, 36010, 36685, 37182, 37775, 38246, 38288, 38631, 39349, 39487, 40420, 40472, 40849, 41302, 41456, 42319, 43367]
	uds_3 = [190, 394, 1620, 2166, 2211, 4059, 4128, 4701, 4706, 5155, 7516, 9207, 9367, 11558, 12010, 13108, 14409, 14996, 15598, 16022, 16709, 17838, 19068, 19126, 20694, 20704, 20770, 21998, 22227, 22416, 22685, 23692, 26581, 28087, 29179, 29461, 30196, 30916, 31615, 32147, 32351, 32947, 34150, 34817, 36247, 38640, 39126, 39624, 42529, 42571, 42812, 43667]
	
	for i in aegis_1:
		print "a" + str(i)



print "end"