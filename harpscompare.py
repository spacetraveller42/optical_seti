# Code to download HARPS spectral and CCD data, and plot it side-by-side.  Steals code and ideas from
# Benjamin Fields (fields_benjamin@wheatoncollege.edu)
#
# Author: Jason Goodman (goodman_jason@wheatoncollege.edu)

import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astroquery.eso import Eso

# Download spectrum from ESO HARPS archive, given "Arcfile" name (leave out the ".fits")
def download_spectrum(specfilearc,replace_underscores=True):
    eso = Eso()
    eso.login("goodmanj")
    if replace_underscores:
        specfilearc = specfilearc.replace("_",":") # Replace underscores with colons
    specfilename = eso.retrieve_data(specfilearc)
    return specfilename

# Extract spectral data from spectral FITS file, given filename.
def read_spectrum(specfilename):
    specfits = fits.open(specfilename)
    return (specfits[1].data[0][0], specfits[1].data[0][1])

# Look in FITS header to find arcfile name of raw image file associated with a given spectral file, and download it.
def download_associated_raw(specfilename,decompress=True):
    specfits = fits.open(specfilename)
    outfilename = download_raw(specfits[0].header['PROV1'][:-5],decompress=decompress)
    specfits.close()
    return outfilename

# Download raw image file, given its "Arcfile" name.  
# Set decompress=true to use unlzw3 module to decompress (slow)
def download_raw(rawfilearc,decompress=True):
    import unlzw3
    eso = Eso()
    eso.login("goodmanj")
    print("Downloading "+rawfilearc)
    rawfilename = eso.retrieve_data(rawfilearc,unzip=False)
    if decompress:
        # Uncompress using unlzw3
        uncompressedfilename = rawfilename[:-2]
        print("Uncompressing "+rawfilename+" to "+uncompressedfilename+"...")
        fin = open(rawfilename,"rb")
        fout = open(uncompressedfilename,"wb")
        fout.write(unlzw3.unlzw(fin.read(-1)))
        fout.close()
        fin.close()
        return uncompressedfilename
    else:
        return rawfilename

# Read spectral layout information (location of spectral orders on CCD chip and what wavelengths they correspond to)
start_wavelength = [float(x.split('\t')[7]) for x in open("harps_spectralpositioning.txt").readlines()]
end_wavelength = [float(x.split('\t')[8]) for x in open("harps_spectralpositioning.txt").readlines()]
y_spectral_locations = [int(x.split('\t')[2]) for x in open("harps_spectralpositioning.txt").readlines()]

# Figure out which spectral order (where on the CCD) a given wavelength is.

def find_order(lamb):  # lamb in nanometers
    # Scan through wavelength list, find the index of the first start/end pair that straddles lamb
    ix = np.where([lamb > start and lamb < end for (start,end) in zip(start_wavelength,end_wavelength)])[0][0]
    if (lamb < 531): # Use "Linda", fits image #1
        ccd=1
        starty = y_spectral_locations[ix]+25  # Orders are offset from what's claimed on 
        endy = y_spectral_locations[ix]+125   # https://www.eso.org/sci/facilities/lasilla/instruments/harps/inst/spec_form.html.
    else: # use "Jasmin", fits image #2
        ccd=2
        starty = y_spectral_locations[ix]  # Different offsets for the two chips
        endy = y_spectral_locations[ix]+100        
    return(ccd,starty,endy, start_wavelength[ix],end_wavelength[ix])

# Plot spectrum and CCD image file on the same figure.  specfilename and rawfilename should be files already
# downloaded to the local computer.  This function plots one full spectral order containing wavelength
# "lamb" (in nanometers)  Countmin and countmax are optional arguments describing the scale of the spectral graph, and the color
# range of the CCD image.
def compare_spec_to_raw(specfilename,rawfilename,lamb,countmin=-1,countmax=-1):
# Open spectral and raw image files
    specfits = fits.open(specfilename)
    rawfits = fits.open(rawfilename)
    wave = specfits[1].data[0][0]
    arr1 = specfits[1].data[0][1]
# Find the spectral order for the desired wavelength so we can look at the right part of the CCD
    (ccd,starty,endy,lam_start,lam_end) = find_order(lamb)
    im = rawfits[ccd].data
# Subplot: Upper panel shows spectrum, lower panel shows CCD image
    plt.subplot(2,1,1)
    if (countmin==-1): #if no value for countmin and countmax were provided
        countmin = np.mean(arr1)-2*np.std(arr1)  # Try to make plot axis limits reasonable based on range of data
        countmax = np.mean(arr1)+6*np.std(arr1)
    plt.plot(wave,arr1)
    plt.plot([lamb*10, lamb*10],[countmin, countmax],':')  # Plot a dotted line showing the location of our spike
    plt.xlim(lam_start*10+1,lam_end*10)  # nanometers to angstroms
    plt.xlabel("Wavelength (Angstroms)")
    plt.ylim([countmin,countmax]) 
    plt.subplot(2,1,2)
    plt.imshow(im.transpose(),vmin=countmin+np.min(im),vmax=countmax+np.min(im),aspect="auto")
    plt.xlim(160,4096)
    plt.xlabel("CCD pixels")
    plt.gca().invert_xaxis()  # CCD image has wavelength in reverse order
    plt.ylim(starty,endy)
    plt.show()