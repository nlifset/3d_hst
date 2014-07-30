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
import matplotlib.image as mpimg
from matplotlib.gridspec import GridSpec

def make_plot(id = 123, field = 'aegis'):
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
	

	data = ascii.read("%s..pz" % (id))
	z = data["z"]
	pz = data["pz"]
	
	field = field.upper()
	id = str(id)
	
	os.chdir("/Volumes/TOSHIBA EXT/3d_hst/%s-RGB_v4.0" % (field))
	
	if len(str(id)) == 5:
		img=mpimg.imread("/Volumes/TOSHIBA EXT/3d_hst/%s-RGB_v4.0/%s_%s_vJH_6.png" % (field, field, id))
	elif len(str(id)) == 4:
		id = "0" + id
		img=mpimg.imread("/Volumes/TOSHIBA EXT/3d_hst/%s-RGB_v4.0/%s_%s_vJH_6.png" % (field, field, id))
	else:
		id = "00" + id
		img=mpimg.imread("/Volumes/TOSHIBA EXT/3d_hst/%s-RGB_v4.0/%s_%s_vJH_6.png" % (field, field, id))
	
	gs=GridSpec(2,2, hspace=0.5)
	a1=pylab.subplot(gs[0,:])
	a2=pylab.subplot(gs[1,0])
	a3=pylab.subplot(gs[1,1])
	
	a1.plot(lam2, flux2)
	a1.scatter(lam1, flux1)
	a1.errorbar(lam1, flux1, yerr=error, linestyle="none")
	a1.set_xlim([2000,140000])
	a1.set_ylim(ymin=0)
	a1.set_xscale("log")
	a1.set_xlabel("Wavelength")
	a1.set_ylabel("Flux")
	
	a2.plot(z, pz)
	a2.set_xlabel("z")
	a2.set_ylabel("pz")
	a2.set_yscale("log")
	a2.set_ylim([0.01,10])
	a2.set_title("P of Z")
	
	a3.imshow(img)
	a3.axis("off")
	a3.set_title("Image of galaxy")
	
	pylab.suptitle("id %s of field %s" % (id, field), fontsize=19)
	pylab.text(0, -110, "Wavelength vs Flux", fontsize=16)
	
	os.chdir("/Volumes/TOSHIBA EXT/3d_hst")
	
def make_plot2(id = 123, field = 'aegis'):
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
	

	data = ascii.read("%s..pz" % (id))
	z = data["z"]
	pz = data["pz"]
	
	field = field.upper()
	id = str(id)
	
	os.chdir("/Volumes/TOSHIBA EXT/3d_hst/%s-RGB_v4.0" % (field))
	
	fig,pylab.axes = pylab.subplots(1, 2)

	a1 = pylab.axes[0]
	a2 = pylab.axes[1]
	
	a1.plot(lam2, flux2)
	a1.scatter(lam1, flux1)
	a1.errorbar(lam1, flux1, yerr=error, linestyle="none")
	a1.set_xlim([2000,140000])
	a1.set_ylim(ymin=0)
	a1.set_xscale("log")
	a1.set_xlabel("Wavelength")
	a1.set_ylabel("Flux")
	
	a2.plot(z, pz)
	a2.set_xlabel("z")
	a2.set_ylabel("pz")
	a2.set_yscale("log")
	a2.set_ylim([0.01,10])
	a2.set_title("P of Z")
	
	pylab.suptitle("id %s of field %s" % (id, field), fontsize=19)
	pylab.text(0, -110, "Wavelength vs Flux", fontsize=16)
	
	os.chdir("/Volumes/TOSHIBA EXT/3d_hst")


print "end"