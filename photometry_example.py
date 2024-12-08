# Example: this example creates three color-magnitude diagrams from an image of the Pleiades(M45) I took and stacked. The zipped image is included in the repository

# importing the image
import numpy as np
from astropy.io import fits
from photutils.detection import IRAFStarFinder as iraf
import matplotlib.pyplot as plt
from photometry import *

fits_data, header = fits.getdata('M45_final.FTS', header=True)
sectionR = fits_data[0,:,:]     # can adjust these to get certain regions of the image
sectionG = fits_data[1,:,:]
sectionB = fits_data[2,:,:]

# Finding the stars
iraffind = iraf(fwhm = 3.0, threshold = 1e-5, roundlo=-1.0, roundhi=1.0)
sources = iraffind(sectionR-np.median(sectionR))

# Establishing the magnitudes, calculating magnitude error and finding 3 sigma upper limit 
I_mag, uncI, I_sul = photometry(section=sectionR, sources=sources, r_ap=20.0,r_in=28,r_out=38,exptime=header['exptime'],zeropoint=0,color='red')
V_mag, uncV, V_sul = photometry(section=sectionG, sources=sources, r_ap=20.0,r_in=28,r_out=38,exptime=header['exptime'],zeropoint=0,color='green')
B_mag, uncB, B_sul = photometry(section=sectionB, sources=sources, r_ap=20.0,r_in=28,r_out=38,exptime=header['exptime'],zeropoint=0,color='blue')

print('Three sigma upper limit I = ', I_sul, 'Three sigma upper limit V = ', V_sul, 'Three sigma upper limit B = ', B_sul)

BV = np.subtract(B_mag,V_mag)
diagram, (I_mags,V_mags,B_mags) = plt.subplots(1,3,sharey=True)

def subplot(name,xlabel,ylabel,x,y,xerr,yerr,errcolor,plotcolor):   # defines a function to create the subplots
    name.set_xlabel(xlabel)
    name.set_ylabel(ylabel)
    name.invert_yaxis()
    name.errorbar(x,y,xerr=xerr,yerr=yerr,ls='none', color=errcolor)
    name.plot(x,y,'o',color=plotcolor)

subplot(name=I_mags,xlabel='B-V',ylabel='I magnitude',x=BV,y=I_mag,xerr=(np.array(uncB)+np.array(uncV)),yerr=uncI,errcolor='purple',plotcolor='Red')
subplot(name=V_mags,xlabel='B-V',ylabel='V magnitude',x=BV,y=V_mag,xerr=(np.array(uncB)+np.array(uncV)),yerr=uncV,errcolor='purple',plotcolor='Green')
subplot(name=B_mags,xlabel='B-V',ylabel='B magnitude',x=BV,y=B_mag,xerr=(np.array(uncB)+np.array(uncV)),yerr=uncB,errcolor='purple',plotcolor='Blue')
plt.show()
