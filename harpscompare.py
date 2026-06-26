# Code to download HARPS spectral and CCD data, and plot it side-by-side.  Steals code and ideas from
# Benjamin Fields (fields_benjamin@wheatoncollege.edu)
#
# Author: Jason Goodman (goodman_jason@wheatoncollege.edu)

import os
import unlzw3
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.image import NonUniformImage
from astropy.io import fits
from astroquery.eso import Eso
from pathlib import Path

# Configure these lines for your local setup
eso_login = "goodmanj"

# Module-level ESO session — login happens once, then gets reused.
_eso = None

def get_eso():
    """Return a logged-in Eso instance, creating and logging in only on the first call."""
    global _eso
    if _eso is None:
        _eso = Eso()
        _eso.login(username=eso_login)
    return _eso


# Path to spectral positioning file - update this for your local setup or set
# the HARPS_SPECPOS_FILE environment variable to override.
_default_specpos = os.path.join(os.path.dirname(os.path.abspath(__file__)), "harps_spectralpositioning.txt")
harps_spectralpositioning_file = os.environ.get("HARPS_SPECPOS_FILE", _default_specpos)

# Download spectrum from ESO HARPS archive, given "Arcfile" name (leave out the ".fits")
def download_spectrum(specfilearc, replace_underscores=True):
    eso = get_eso()
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
# Set decompress=True to use unlzw3 to uncompress the files.
def download_raw(rawfilearc, decompress=True):
    eso = get_eso()
    print("Downloading " + rawfilearc)
    result = eso.retrieve_data(rawfilearc, unzip=False)

    if not result:
        raise RuntimeError(f"Download failed for {rawfilearc}: retrieve_data returned an empty list")

    # If retrieve_data returned a list, extract the single path.
    if isinstance(result,list):
        rawfilename = Path(result[0])
    elif isinstance(result,str):
        rawfilename = Path(result)
    else:
        raise RuntimeError(f"Unrecognized response for {rawfilearc}: {result}")

    if decompress:
        uncompressedfilename = rawfilename.with_suffix("")
        # If the .Z file is gone but the decompressed file already exists (e.g. from
        # a previous run), skip decompression and return the existing file.
        if not rawfilename.exists() and uncompressedfilename.exists():
            print(f"Already decompressed: {uncompressedfilename}")
            return str(uncompressedfilename)
        print(f"Downloaded: {rawfilename}")
        print(f"Uncompressing {rawfilename} to {uncompressedfilename}...")
        # Uncompress using unlzw3 for cross-platform compatibility
        fin = open(rawfilename,"rb")
        fout = open(uncompressedfilename,"wb")
        fout.write(unlzw3.unlzw(fin.read(-1)))
        fout.close()
        fin.close()
        return str(uncompressedfilename)
    else:
        return str(rawfilename)

# Read spectral layout information (location of spectral orders on CCD chip and what wavelengths they correspond to)
# Column indices in harps_spectralpositioning.txt: 1=mid wavelength, 2=y pixel location,
#   7=start wavelength, 8=end wavelength
_specpos_lines = open(harps_spectralpositioning_file).readlines()
start_wavelength = [float(x.split('\t')[7]) for x in _specpos_lines]
end_wavelength = [float(x.split('\t')[8]) for x in _specpos_lines]
mid_wavelength = [float(x.split('\t')[1]) for x in _specpos_lines]
y_spectral_locations = [int(x.split('\t')[2]) for x in _specpos_lines]

# Figure out which spectral orders (where on the CCD) a given wavelength is.
# lamb: wavelength in *nanometers*
# optional order_shift: Offset the calculated spectral order (integer)
def find_order(lamb,order_shift=0):  # lamb in nanometers
    if (lamb < 531): # Use "Linda", fits image #1
        ccd=1
    else: # use "Jasmin", fits image #2
        ccd=2
    # Scan through wavelength list, find the index of the first start/end pair that straddles lamb[0][0]
    ixes = np.where([lamb > start and lamb < end for (start,end) in zip(start_wavelength,end_wavelength)])[0] + order_shift
    order_data = [(ccd,y_spectral_locations[ix],y_spectral_locations[ix]+100,start_wavelength[ix],end_wavelength[ix],mid_wavelength[ix]) for ix in ixes]
    return order_data

# Doppler shift function
def doppler(wave,v):
    c = 2.998e5 # km/s
    beta = v/c
    return wave*np.sqrt((1-beta)/(1+beta))

# Plot spectrum and CCD image file on the same figure.  Plots one full spectral order containing a given wavelength.
# specfilename: file name of downloaded HARPS spectrum file
# rawfilename: file name of downloaded raw CCD image
# lamb: target wavelength, in *nanometers*
# lamb_range: Range of wavelengths to plot, in *nanometers*
# spec_countmin, spec_countmax: minimum and maximum range of spectral plot.  Optional: if left blank or set to -1,
#   "reasonable" values will be autoselected.
# raw_countmin, raw_countmax: minimum and maximum range of ccd image.  Optional: if left blank or set to -1,
#   will use "reasonable" percentile-based color scaling.
# multi_order: plot both spectral orders from the raw image that include the wavelength, if possible
# doppler_shift: doppler-shift CCD image wavelengths to stellar frame

def compare_spec_to_raw(specfilename,rawfilename,lamb,lamb_range = [],
                        spec_countmin=-1,spec_countmax=-1,raw_countmin=-1,raw_countmax=-1,
                        order_shift=0,multi_order=False,doppler_shift=False):
    # Open spectral and raw image files
    specfits = fits.open(specfilename)
    rawfits = fits.open(rawfilename)
    wave = specfits[1].data[0][0]
    arr1 = specfits[1].data[0][1]
    # Find the spectral order(s) for the desired wavelength so we can look at the right part of the CCD.  There may be one or two.
    order_info = find_order(lamb,order_shift)
    n_orders = len(order_info) # Number of orders in which this wavelength appears
    if not multi_order:
        n_orders = 1           # Only show the first order

    # Subplot: Upper panel shows spectrum, lower panel(s) show CCD image
    # Show spectrum
    plt.subplot(n_orders+1,1,1)
    if (spec_countmin==-1):
        spec_countmin = np.mean(arr1)-2*np.std(arr1)  # Try to make plot axis limits reasonable based on range of data
        spec_countmax = np.mean(arr1)+6*np.std(arr1)
    plt.plot(wave,arr1)
    plt.plot([lamb*10, lamb*10],[spec_countmin, spec_countmax],':')  # Plot a dotted line showing the location of our spike
    plt.ylim([spec_countmin,spec_countmax]) 
    plt.title(specfits[0].header['ARCFILE'])
    specax = plt.gca()

    # Show raw image(s)
    for order in range(n_orders):
        (ccd,starty,endy,lam_start,lam_end,lam_mid) = order_info[order]
        # Doppler-shift CCD order wavelengths into star's reference frame (-berv)
        if (doppler_shift):
            berv = specfits[0].header['HIERARCH ESO DRS BERV']
            lam_start = doppler(lam_start,-berv)
            lam_end = doppler(lam_end,-berv)
            lam_mid = doppler(lam_mid,-berv)
        if len(rawfits)<2:  # 2003 FITS files are in a different format
            if ccd==1:
                im = rawfits[0].data[:,0:2148]
            else:
                im = rawfits[0].data[:,2149:]
        else:
            im = rawfits[ccd].data
        rows,cols = im.shape
        if (raw_countmin==-1):
            # Set image color scale limits based on percentiles.
            raw_countmin = np.percentile(im[:,starty:endy],20)
            raw_countmax = np.percentile(im[:,starty:endy],99)
        plt.subplot(n_orders+1,1,order+2,sharex=specax) # All panels zoom together
        # CCD data not regullarly spaced: fit CCD wavelength data to a 2nd-degree polynomioal and do a nonuniform image
        y = np.linspace(1,cols,cols)
        p = np.polynomial.Polynomial.fit([1,4096/2,4096],[lam_start*10,lam_mid*10,lam_end*10],2)
        lam_pix,lam_grid = p.linspace(n=4096)
        dispim = NonUniformImage(plt.gca())    
        dispim.set_data(lam_grid, y, np.flip(im.transpose(),axis=1))
        dispim.set_clim(vmin=raw_countmin,vmax=raw_countmax)
        plt.gca().add_image(dispim)
        plt.gca().invert_yaxis()  # Higher orders up
        plt.ylim(starty,endy)
        if (lamb_range):
            plt.xlim(lamb_range[0]*10,lamb_range[1]*10)
        else:
            plt.xlim(lam_start*10,lam_end*10)
        plt.title(Path(rawfilename).name)
        if (order == n_orders-1):
            plt.xlabel("Wavelength (Angstroms)")
    plt.subplots_adjust(hspace=0.5)
    plt.show()
    

def plot_raw_image(specfilename,rawfilename,lamb,lamb_range = [],countmin=-1,countmax=-1):
    specfits = fits.open(specfilename)
    rawfits = fits.open(rawfilename)
    arr1 = specfits[1].data[0][1]
    (ccd,starty,endy,lam_start,lam_end,lam_mid) = find_order(lamb)[0]
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
