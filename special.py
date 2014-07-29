print "start"
import os
os.chdir("/Volumes/TOSHIBA EXT/3d_hst")
from astropy.io import ascii
import numpy as np
import astropy.units as u
from astropy.io import fits
import astropy.constants
import matplotlib.pyplot as plt
import pylab
import math

def make_plot(id = 123, field = 'aegis'):
	dataz = ascii.read("%s_3dhst.v4.0.sfr" % (field))
	dataz_id = dataz[dataz["id"] == id]
	z = dataz_id["z"]
	os.chdir("/Volumes/TOSHIBA EXT/3d_hst/noah/%s" % (field))
	data1 = ascii.read("%s..obs_sed" % (id))
	data2 = ascii.read("%s..temp_sed" % (id))
	lam1 = data1["lambda"]
	flux_1 = data1["flux_cat"]
	lam2 = data2["lambda"]
	flux_2 = data2["tempflux"]
	
	factor = 3.0*(10.0**5.56)
	flux1 = (flux_1*(lam1**-2.0))*factor
	flux2 = (flux_2*(lam2**-2.0))*factor
	
	error_ = data1["err_full"]
	error = (error_*(lam1**-2.0))*factor
	pylab.plot(lam2, flux2)
	pylab.scatter(lam1, flux1)
	pylab.errorbar(lam1, flux1, yerr=error, linestyle="none")
	pylab.xlim([2000,140000])
	pylab.ylim(ymin=0)
	pylab.xscale("log")
	pylab.xlabel("Wavelength")
	pylab.ylabel("Flux")
	pylab.suptitle("id %s of field %s" % (id, field))
	pylab.ion()
	pylab.show()
	os.chdir("/Volumes/TOSHIBA EXT/3d_hst")

def make_zplot(id = 123, field = 'aegis'):
	os.chdir("/Volumes/TOSHIBA EXT/3d_hst/noah/%s" % (field))
	data = ascii.read("%s..pz" % (id))
	z = data["z"]
	pz = data["pz"]
	pylab.plot(z, pz)
	pylab.xlabel("z")
	pylab.ylabel("pz")
	pylab.yscale("log")
	pylab.ylim([0.01,10])
	pylab.suptitle("id %s of field %s" % (id, field))
	os.chdir("/Volumes/TOSHIBA EXT/3d_hst")


print "end"