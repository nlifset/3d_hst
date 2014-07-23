print "start"
import os
os.chdir("/Volumes/TOSHIBA EXT/3d_hst/noah")
from astropy.io import ascii
import numpy as np
import astropy.units as u
from astropy.io import fits
import astropy.constants
import matplotlib.pyplot as plt
import pylab
import operator 

#chunk1#

data18935 = ascii.read("18935..obs_sed")
data3740 = ascii.read("3740..obs_sed")
data37194 = ascii.read("37194..obs_sed")
data33990 = ascii.read("33990..obs_sed")
data7002 = ascii.read("7002..obs_sed")
data11221 = ascii.read("11221..obs_sed")
data25242 = ascii.read("25242..obs_sed")
data10509 = ascii.read("10509..obs_sed")
data33712 = ascii.read("33712..obs_sed")
data8664 = ascii.read("8664..obs_sed")

data_18935 = ascii.read("18935..temp_sed")
data_3740 = ascii.read("3740..temp_sed")
data_37194 = ascii.read("37194..temp_sed")
data_33990 = ascii.read("33990..temp_sed")
data_7002 = ascii.read("7002..temp_sed")
data_11221 = ascii.read("11221..temp_sed")
data_25242 = ascii.read("25242..temp_sed")
data_10509 = ascii.read("10509..temp_sed")
data_33712 = ascii.read("33712..temp_sed")
data_8664 = ascii.read("8664..temp_sed")

lam18935 = data18935["lambda"]
lam3740 = data3740["lambda"]
lam37194 = data37194["lambda"]
lam33990 = data33990["lambda"]
lam7002 = data7002["lambda"]
lam11221 = data11221["lambda"]
lam25242 = data25242["lambda"]
lam10509 = data10509["lambda"]
lam33712 = data33712["lambda"]
lam8664 = data8664["lambda"]

lam_18935 = data_18935["lambda"]
lam_3740 = data_3740["lambda"]
lam_37194 = data_37194["lambda"]
lam_33990 = data_33990["lambda"]
lam_7002 = data_7002["lambda"]
lam_11221 = data_11221["lambda"]
lam_25242 = data_25242["lambda"]
lam_10509 = data_10509["lambda"]
lam_33712 = data_33712["lambda"]
lam_8664 = data_8664["lambda"]

flux18935 = data18935["flux_cat"]
flux3740 = data3740["flux_cat"]
flux37194 = data37194["flux_cat"]
flux33990 = data33990["flux_cat"]
flux7002 = data7002["flux_cat"]
flux11221 = data11221["flux_cat"]
flux25242 = data25242["flux_cat"]
flux10509 = data10509["flux_cat"]
flux33712 = data33712["flux_cat"]
flux8664 = data8664["flux_cat"]

flux_18935 = data_18935["tempflux"]
flux_3740 = data_3740["tempflux"]
flux_37194 = data_37194["tempflux"]
flux_33990 = data_33990["tempflux"]
flux_7002 = data_7002["tempflux"]
flux_11221 = data_11221["tempflux"]
flux_25242 = data_25242["tempflux"]
flux_10509 = data_10509["tempflux"]
flux_33712 = data_33712["tempflux"]
flux_8664 = data_8664["tempflux"]

#chunk2#

data25354 = ascii.read("25354..obs_sed")
data16238 = ascii.read("16238..obs_sed")
data9849 = ascii.read("9849..obs_sed")
data26850 = ascii.read("26850..obs_sed")
data6424 = ascii.read("6424..obs_sed")
data17527 = ascii.read("17527..obs_sed")
data22126 = ascii.read("22126..obs_sed")
data18045 = ascii.read("18045..obs_sed")
data30421 = ascii.read("30421..obs_sed")
data38065 = ascii.read("38065..obs_sed")

data_25354 = ascii.read("25354..temp_sed")
data_16238 = ascii.read("16238..temp_sed")
data_9849 = ascii.read("9849..temp_sed")
data_26850 = ascii.read("26850..temp_sed")
data_6424 = ascii.read("6424..temp_sed")
data_17527 = ascii.read("17527..temp_sed")
data_22126 = ascii.read("22126..temp_sed")
data_18045 = ascii.read("18045..temp_sed")
data_30421 = ascii.read("30421..temp_sed")
data_38065 = ascii.read("38065..temp_sed")

lam25354 = data25354["lambda"]
lam16238 = data16238["lambda"]
lam9849 = data9849["lambda"]
lam26850 = data26850["lambda"]
lam6424 = data6424["lambda"]
lam17527 = data17527["lambda"]
lam22126 = data22126["lambda"]
lam18045 = data18045["lambda"]
lam30421 = data30421["lambda"]
lam38065 = data38065["lambda"]

lam_25354 = data_25354["lambda"]
lam_16238 = data_16238["lambda"]
lam_9849 = data_9849["lambda"]
lam_26850 = data_26850["lambda"]
lam_6424 = data_6424["lambda"]
lam_17527 = data_17527["lambda"]
lam_22126 = data_22126["lambda"]
lam_18045 = data_18045["lambda"]
lam_30421 = data_30421["lambda"]
lam_38065 = data_38065["lambda"]

flux25354 = data25354["flux_cat"]
flux16238 = data16238["flux_cat"]
flux9849 = data9849["flux_cat"]
flux26850 = data26850["flux_cat"]
flux6424 = data6424["flux_cat"]
flux17527 = data17527["flux_cat"]
flux22126 = data22126["flux_cat"]
flux18045 = data18045["flux_cat"]
flux30421 = data30421["flux_cat"]
flux38065 = data38065["flux_cat"]

flux_25354 = data_25354["tempflux"]
flux_16238 = data_16238["tempflux"]
flux_9849 = data_9849["tempflux"]
flux_26850 = data_26850["tempflux"]
flux_6424 = data_6424["tempflux"]
flux_17527 = data_17527["tempflux"]
flux_22126 = data_22126["tempflux"]
flux_18045 = data_18045["tempflux"]
flux_30421 = data_30421["tempflux"]
flux_38065 = data_38065["tempflux"]

#chunk3#

data2918 = ascii.read("2918..obs_sed")
data38187 = ascii.read("38187..obs_sed")
data24333 = ascii.read("24333..obs_sed")
data18922 = ascii.read("18922..obs_sed")
data33863 = ascii.read("33863..obs_sed")
data31169 = ascii.read("31169..obs_sed")
data23234 = ascii.read("23234..obs_sed")
data15709 = ascii.read("15709..obs_sed")
data28310 = ascii.read("28310..obs_sed")
data15088 = ascii.read("15088..obs_sed")

data_2918 = ascii.read("2918..temp_sed")
data_38187 = ascii.read("38187..temp_sed")
data_24333 = ascii.read("24333..temp_sed")
data_18922 = ascii.read("18922..temp_sed")
data_33863 = ascii.read("33863..temp_sed")
data_31169 = ascii.read("31169..temp_sed")
data_23234 = ascii.read("23234..temp_sed")
data_15709 = ascii.read("15709..temp_sed")
data_28310 = ascii.read("28310..temp_sed")
data_15088 = ascii.read("15088..temp_sed")

lam2918 = data2918["lambda"]
lam38187 = data38187["lambda"]
lam24333 = data24333["lambda"]
lam18922 = data18922["lambda"]
lam33863 = data33863["lambda"]
lam31169 = data31169["lambda"]
lam23234 = data23234["lambda"]
lam15709 = data15709["lambda"]
lam28310 = data28310["lambda"]
lam15088 = data15088["lambda"]

lam_2918 = data_2918["lambda"]
lam_38187 = data_38187["lambda"]
lam_24333 = data_24333["lambda"]
lam_18922 = data_18922["lambda"]
lam_33863 = data_33863["lambda"]
lam_31169 = data_31169["lambda"]
lam_23234 = data_23234["lambda"]
lam_15709 = data_15709["lambda"]
lam_28310 = data_28310["lambda"]
lam_15088 = data_15088["lambda"]

flux2918 = data2918["flux_cat"]
flux38187 = data38187["flux_cat"]
flux24333 = data24333["flux_cat"]
flux18922 = data18922["flux_cat"]
flux33863 = data33863["flux_cat"]
flux31169 = data31169["flux_cat"]
flux23234 = data23234["flux_cat"]
flux15709 = data15709["flux_cat"]
flux28310 = data28310["flux_cat"]
flux15088 = data15088["flux_cat"]

flux_2918 = data_2918["tempflux"]
flux_38187 = data_38187["tempflux"]
flux_24333 = data_24333["tempflux"]
flux_18922 = data_18922["tempflux"]
flux_33863 = data_33863["tempflux"]
flux_31169 = data_31169["tempflux"]
flux_23234 = data_23234["tempflux"]
flux_15709 = data_15709["tempflux"]
flux_28310 = data_28310["tempflux"]
flux_15088 = data_15088["tempflux"]

#chunk4#

data179 = ascii.read("179..obs_sed")
data8790 = ascii.read("8790..obs_sed")
data41033 = ascii.read("41033..obs_sed")
data21306 = ascii.read("21306..obs_sed")
data3048 = ascii.read("3048..obs_sed")
data1849 = ascii.read("1849..obs_sed")
data27584 = ascii.read("27584..obs_sed")
data19028 = ascii.read("19028..obs_sed")
data38731 = ascii.read("38731..obs_sed")
data37799 = ascii.read("37799..obs_sed")

data_179 = ascii.read("179..temp_sed")
data_8790 = ascii.read("8790..temp_sed")
data_41033 = ascii.read("41033..temp_sed")
data_21306 = ascii.read("21306..temp_sed")
data_3048 = ascii.read("3048..temp_sed")
data_1849 = ascii.read("1849..temp_sed")
data_27584 = ascii.read("27584..temp_sed")
data_19028 = ascii.read("19028..temp_sed")
data_38731 = ascii.read("38731..temp_sed")
data_37799 = ascii.read("37799..temp_sed")

lam179 = data179["lambda"]
lam8790 = data8790["lambda"]
lam41033 = data41033["lambda"]
lam21306 = data21306["lambda"]
lam3048 = data3048["lambda"]
lam1849 = data1849["lambda"]
lam27584 = data27584["lambda"]
lam19028 = data19028["lambda"]
lam38731 = data38731["lambda"]
lam37799 = data37799["lambda"]

lam_179 = data_179["lambda"]
lam_8790 = data_8790["lambda"]
lam_41033 = data_41033["lambda"]
lam_21306 = data_21306["lambda"]
lam_3048 = data_3048["lambda"]
lam_1849 = data_1849["lambda"]
lam_27584 = data_27584["lambda"]
lam_19028 = data_19028["lambda"]
lam_38731 = data_38731["lambda"]
lam_37799 = data_37799["lambda"]

flux179 = data179["flux_cat"]
flux8790 = data8790["flux_cat"]
flux41033 = data41033["flux_cat"]
flux21306 = data21306["flux_cat"]
flux3048 = data3048["flux_cat"]
flux1849 = data1849["flux_cat"]
flux27584 = data27584["flux_cat"]
flux19028 = data19028["flux_cat"]
flux38731 = data38731["flux_cat"]
flux37799 = data37799["flux_cat"]

flux_179 = data_179["tempflux"]
flux_8790 = data_8790["tempflux"]
flux_41033 = data_41033["tempflux"]
flux_21306 = data_21306["tempflux"]
flux_3048 = data_3048["tempflux"]
flux_1849 = data_1849["tempflux"]
flux_27584 = data_27584["tempflux"]
flux_19028 = data_19028["tempflux"]
flux_38731 = data_38731["tempflux"]
flux_37799 = data_37799["tempflux"]



#plot (outlier not included: id37194, chunk1)#

pylab.plot(lam_3740, flux_3740, color="r", label="chunk1, z < 0.5")
pylab.plot(lam_33990, flux_33990, color="r")
pylab.plot(lam_7002, flux_7002, color="r")
pylab.plot(lam_11221, flux_11221, color="r")
pylab.plot(lam_25242, flux_25242, color="r")
pylab.plot(lam_10509, flux_10509, color="r")
pylab.plot(lam_33712, flux_33712, color="r")
pylab.plot(lam_8664, flux_8664, color="r")
pylab.plot(lam_18935, flux_18935, color="r")

pylab.plot(lam_25354, flux_25354, color="g", label="chunk2, 0.5 < z < 1.5")
pylab.plot(lam_16238, flux_16238, color="g")
pylab.plot(lam_9849, flux_9849, color="g")
pylab.plot(lam_26850, flux_26850, color="g")
pylab.plot(lam_6424, flux_6424, color="g")
pylab.plot(lam_17527, flux_17527, color="g")
pylab.plot(lam_22126, flux_22126, color="g")
pylab.plot(lam_18045, flux_18045, color="g")
pylab.plot(lam_30421, flux_30421, color="g")
pylab.plot(lam_38065, flux_38065, color="g")

pylab.plot(lam_2918, flux_2918, color="b", label="chunk3, 1.5 < z < 2.5")
pylab.plot(lam_38187, flux_38187, color="b")
pylab.plot(lam_24333, flux_24333, color="b")
pylab.plot(lam_18922, flux_18922, color="b")
pylab.plot(lam_33863, flux_33863, color="b")
pylab.plot(lam_31169, flux_31169, color="b")
pylab.plot(lam_23234, flux_23234, color="b")
pylab.plot(lam_15709, flux_15709, color="b")
pylab.plot(lam_28310, flux_28310, color="b")
pylab.plot(lam_15088, flux_15088, color="b")

pylab.plot(lam_179, flux_179, color="k", label="chunk4, z > 2.5")
pylab.plot(lam_8790, flux_8790, color="k")
pylab.plot(lam_41033, flux_41033, color="k")
pylab.plot(lam_21306, flux_21306, color="k")
pylab.plot(lam_3048, flux_3048, color="k")
pylab.plot(lam_1849, flux_1849, color="k")
pylab.plot(lam_27584, flux_27584, color="k")
pylab.plot(lam_19028, flux_19028, color="k")
pylab.plot(lam_38731, flux_38731, color="k")
pylab.plot(lam_37799, flux_37799, color="k")



pylab.xlim([-1000,200000])
pylab.xlabel("Lambda", fontsize=16)
pylab.ylabel("Flux", fontsize=16)
pylab.suptitle("Ten most massive galaxies of each z-bin", fontsize=18)
pylab.legend()
pylab.ion()
pylab.show()
pylab.savefig("wavelength_flux_without")


os.chdir("/Volumes/TOSHIBA EXT/3d_hst")
print "end"