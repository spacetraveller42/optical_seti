# Code to download HARPS spectral and CCD data, and plot it side-by-side.  Steals code and ideas from
# Benjamin Fields (fields_benjamin@wheatoncollege.edu)
#
# Author: Jason Goodman (goodman_jason@wheatoncollege.edu)

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.image import NonUniformImage
from astropy.io import fits
from astroquery.eso import Eso
from pathlib import Path

# Configure these lines for your local setup
eso_login = "goodmanj"
opticalseti_lib_path = Path("C:/Users/goodm/OneDrive/Documents/optical_seti")

# Download spectrum from ESO HARPS archive, given "Arcfile" name (leave out the ".fits")
def download_spectrum(specfilearc,replace_underscores=True):
    eso = Eso()
    eso.login(eso_login)
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
    eso.login(eso_login)
    print("Downloading "+rawfilearc)
    rawfilename = eso.retrieve_data(rawfilearc,unzip=False)
    print(rawfilename)
    if decompress:
        # Uncompress using unlzw3
        uncompressedfilename = str(rawfilename)[:-2]
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
start_wavelength = [float(x.split('\t')[7]) for x in open(opticalseti_lib_path / "harps_spectralpositioning.txt").readlines()]
end_wavelength = [float(x.split('\t')[8]) for x in open(opticalseti_lib_path / "harps_spectralpositioning.txt").readlines()]
mid_wavelength = [float(x.split('\t')[1]) for x in open(opticalseti_lib_path / "harps_spectralpositioning.txt").readlines()]
y_spectral_locations = [int(x.split('\t')[2]) for x in open(opticalseti_lib_path / "harps_spectralpositioning.txt").readlines()]

# Figure out which spectral order (where on the CCD) a given wavelength is.
# lamb: wavelength in *nanometers*
# optional order_shift: Offset the calculated spectral order (integer)
def find_order(lamb,order_shift=0):  # lamb in nanometers
    # Scan through wavelength list, find the index of the first start/end pair that straddles lamb[0][0]
    ix = np.where([lamb > start and lamb < end for (start,end) in zip(start_wavelength,end_wavelength)])[0][0] + order_shift
    if (lamb < 531): # Use "Linda", fits image #1
        ccd=1
        starty = y_spectral_locations[ix]  # Orders are offset from what's claimed on 
        endy = y_spectral_locations[ix]+100   # https://www.eso.org/sci/facilities/lasilla/instruments/harps/inst/spec_form.html.
    else: # use "Jasmin", fits image #2
        ccd=2
        starty = y_spectral_locations[ix]  # Different offsets for the two chips
        endy = y_spectral_locations[ix]+100        
    return(ccd,starty,endy, start_wavelength[ix],end_wavelength[ix],mid_wavelength[ix])

# Plot spectrum and CCD image file on the same figure.  Plots one full spectral order containing a given wavelength.
# specfilename: file name of downloaded HARPS spectrum file
# rawfilename: file name of downloaded raw CCD image
# lamb: target wavelength, in *nanometers*
# lamb_range: Range of wavelengths to plot, in *nanometers*
# spec_countmin, spec_countmax: minimum and maximum range of spectral plot.  Optional: if left blank or set to -1,
#   "reasonable" values will be autoselected.
# raw_countmin, raw_countmax: minimum and maximum range of ccd image.  Optional: if left blank or set to -1,
#   will use the same range as spec_countmin, spec_countmax above.
def compare_spec_to_raw(specfilename,rawfilename,lamb,lamb_range = [],
                        spec_countmin=-1,spec_countmax=-1,raw_countmin=-1,raw_countmax=-1,
                        order_shift=0):
    # Open spectral and raw image files
    specfits = fits.open(specfilename)
    rawfits = fits.open(rawfilename)
    wave = specfits[1].data[0][0]
    arr1 = specfits[1].data[0][1]
    # Find the spectral order for the desired wavelength so we can look at the right part of the CCD
    (ccd,starty,endy,lam_start,lam_end,lam_mid) = find_order(lamb,order_shift)
    im = rawfits[ccd].data
    # Subplot: Upper panel shows spectrum, lower panel shows CCD image
    plt.subplot(2,1,1)
    if (spec_countmin==-1):
        spec_countmin = np.mean(arr1)-2*np.std(arr1)  # Try to make plot axis limits reasonable based on range of data
        spec_countmax = np.mean(arr1)+6*np.std(arr1)
    if (raw_countmin==-1):
        raw_countmin = spec_countmin
        raw_countmax = spec_countmax
    plt.plot(wave,arr1)
    plt.plot([lamb*10, lamb*10],[spec_countmin, spec_countmax],':')  # Plot a dotted line showing the location of our spike
    plt.xlim(lam_start*10+1,lam_end*10)  # nanometers to angstroms
    plt.ylim([spec_countmin,spec_countmax]) 
    plt.title(specfits[0].header['ARCFILE'])
    specax = plt.gca()
    plt.subplot(2,1,2,sharex=specax) # Upper and lower panel zoom together
    rows,cols = im.shape
    # CCD data not regullarly spaced: fit CCD wavelength data to a 2nd-degree polynomioal and do a nonuniform image
    y = np.linspace(1,cols,cols)
    p = np.polynomial.Polynomial.fit([1,4096/2,4096],[lam_start*10,lam_mid*10,lam_end*10],2)
    lam_pix,lam_grid = p.linspace(n=4096)
    dispim = NonUniformImage(plt.gca())    
    dispim.set_data(lam_grid, y, np.flip(im.transpose(),axis=1))
    dispim.set_clim(vmin=raw_countmin+np.min(im),vmax=raw_countmax+np.min(im))
    plt.gca().add_image(dispim)
    plt.gca().invert_yaxis()  # Higher orders up
    plt.ylim(starty,endy)
    if (lamb_range):
        plt.xlim(lamb_range[0]*10,lamb_range[1]*10)
    else:
        plt.xlim(lam_start*10,lam_end*10)
    plt.title(Path(rawfilename).name)
    plt.xlabel("Wavelength (Angstroms)")
    plt.show()
    

def plot_raw_image(specfilename,rawfilename,lamb,lamb_range = [],countmin=-1,countmax=-1):
    specfits = fits.open(specfilename)
    rawfits = fits.open(rawfilename)
    arr1 = specfits[1].data[0][1]
    (ccd,starty,endy,lam_start,lam_end,lam_mid) = find_order(lamb)
    if (countmin==-1):
        countmin = np.mean(arr1)-2*np.std(arr1)  # Try to make plot axis limits reasonable based on range of data
        countmax = np.mean(arr1)+6*np.std(arr1)
    im = rawfits[ccd].data
    rows,cols = im.shape
    # CCD data not regullarly spaced: fit CCD wavelength data to a 2nd-degree polynomioal and do a nonuniform image
    y = np.linspace(1,cols,cols)
    p = np.polynomial.Polynomial.fit([1,4096/2,4096],[lam_start*10,lam_mid*10,lam_end*10],2)
    lam_pix,lam_grid = p.linspace(n=4096)
    dispim = NonUniformImage(plt.gca())    
    dispim.set_data(lam_grid, y, np.flip(im.transpose(),axis=1))
    dispim.set_clim(vmin=countmin+np.min(im),vmax=countmax+np.min(im))
    plt.gca().add_image(dispim)
    plt.gca().invert_yaxis()  # Higher orders up
    plt.ylim(starty,endy)
    if (lamb_range):
        plt.xlim(lamb_range[0]*10,lamb_range[1]*10)
    else:
        plt.xlim(lam_start*10,lam_end*10)
    plt.title(Path(rawfilename).name)
    plt.xlabel("Wavelength (Angstroms)")
