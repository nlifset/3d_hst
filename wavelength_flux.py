print "start"
import os
os.chdir("/Users/noah/Documents/3d_hst/noah")
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

pylab.scatter(lam3740, flux3740, color="r", label="id 3740")
pylab.plot(lam_3740, flux_3740, color="r")
pylab.scatter(lam37194, flux37194, color="g", label="id 37194")
pylab.plot(lam_37194, flux_37194, color="g")
pylab.scatter(lam33990, flux33990, color="b", label="id 33990")
pylab.plot(lam_33990, flux_33990, color="b")
pylab.scatter(lam7002, flux7002, color="k", label="id 7002")
pylab.plot(lam_7002, flux_7002, color="k")
pylab.scatter(lam11221, flux11221, color="y", label="id 11221")
pylab.plot(lam_11221, flux_11221, color="y")
pylab.scatter(lam25242, flux25242, color="c", label="id 25242")
pylab.plot(lam_25242, flux_25242, color="c")
pylab.scatter(lam10509, flux10509, color="0.7", label="id 10509")
pylab.plot(lam_10509, flux_10509, color="0.7")
pylab.scatter(lam33712, flux33712, color="m", label="id 33712")
pylab.plot(lam_33712, flux_33712, color="m")
pylab.scatter(lam8664, flux8664, color="purple", label="id 8664")
pylab.plot(lam_8664, flux_8664, color="purple")
pylab.scatter(lam18935, flux18935, color="orange", label="id 18935")
pylab.plot(lam_18935, flux_18935, color="orange")
pylab.xlim([-1000,100000])
pylab.ylim([-500, 12000])
pylab.legend(bbox_to_anchor=(1.1,1.1))
pylab.xlabel("lambda", fontsize=16)
pylab.ylabel("flux", fontsize=16)
pylab.suptitle("chunk1, z < 0.5", fontsize=20)
pylab.ion()
pylab.show()


os.chdir("/Users/noah/Documents/3d_hst")
print "end"