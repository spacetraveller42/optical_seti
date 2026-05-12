# ###################################
# #                                 #
# #    OPTICAL_SETI_FUNCTIONS.PY    #
# #                                 #
# ###################################

# Functions to support Optical SETI searches using HARPS data.
#
# HARPS: https://www.eso.org/sci/facilities/lasilla/instruments/harps.html
#
# This file contains the main search algorithm and support routines.
#
# Authors: Benjamin Fields, Jason Goodman
# Repo: https://github.com/spacetraveller42/optical_seti
#

# general imports
import numpy as np
from matplotlib import pyplot as plt

import astropy
from astropy.io import fits

# curve fitting imports
from astropy.modeling import models
from astropy import units as u
from specutils.spectra import Spectrum1D
from specutils.fitting import fit_lines


from pathlib import Path
eso_cache_path = Path(astropy.config.get_cache_dir()) / "astroquery" / "Eso"

# Speed of light in km/s
_SPEED_OF_LIGHT_KM_S = 299792.458

# ##### 1.  STATISTICS FUNCTIONS

# Calculate running median of data, using window size.
# Output size is arr1 size - window_size
# Use sliding window view for speed, see
# https://numpy.org/devdocs/reference/generated/numpy.lib.stride_tricks.sliding_window_view.html
def running_median(data, window_size):
  return np.median(np.lib.stride_tricks.sliding_window_view(data,window_size),1)

# Calculate running median of arr1, using window size.
# Uses a for loop instead of numpy tricks: much slower.  For testing only.
def running_median_old(arr1,window_size):
    fluxes = np.array(arr1)
    median_iterations = len(arr1) - window_size + 1
    running_median = []
    for i in range(0, median_iterations):
        start = i 
        end = i + window_size
        running_median.append(np.median(fluxes[start:end])) 
    return(running_median)

                  
# Calculate running mean of data, using window size.
# Output size is arr1 size - window_size
# Use sliding window view for speed, see
# https://numpy.org/devdocs/reference/generated/numpy.lib.stride_tricks.sliding_window_view.html
def running_mean(data, window_size):
  return np.mean(np.lib.stride_tricks.sliding_window_view(data,window_size),1)

# Calculate running mean of arr1, using window size.
# Uses a for loop instead of numpy tricks: much slower.  For testing only.
def running_mean_old(arr1, window_size):
    fluxes = np.array(arr1)
    mean_iterations = len(fluxes) - window_size + 1
    running_mean = []
    for i in range(0, mean_iterations):
        start = i 
        end = i + window_size
        running_mean.append(np.mean(fluxes[start:end])) 
    return(running_mean)

# Calculate running standard dev of data, using window size.
# Output size is arr1 size - window_size
# Use sliding window view for speed, see
# https://numpy.org/devdocs/reference/generated/numpy.lib.stride_tricks.sliding_window_view.html
def running_standarddev(data, stwindow):
   return np.std(np.lib.stride_tricks.sliding_window_view(data,stwindow),1)

# Calculate running stdev of arr1, using window size.
# Uses a for loop instead of numpy tricks: much slower.  For testing only.
def running_standarddev_old(arr1, stwindow):
    import numpy as np
    fluxes = np.array(arr1)
    standard_iterations = len(arr1) - stwindow + 1
    running_stdeviation = []
    for i in range(0, standard_iterations):
        start = i 
        end = i + stwindow
        running_stdeviation.append(np.std(fluxes[start:end])) 
    return(running_stdeviation)

# ##### 2.  WAVELENGTH UTILITIES

# Apply a Doppler shift to a wavelength given a radial velocity.
# Input arguments:
#   wavelength: rest-frame wavelength (any unit, e.g. Angstroms or nm)
#   velocity_km_s: radial velocity in km/s (positive = receding, negative = approaching)
#                  For HARPS, pass the BERV (Barycentric Earth Radial Velocity) to
#                  convert from the stellar/barycentric frame to the observer frame.
# Output arguments:
#   shifted wavelength in the same units as the input wavelength

def doppler(wavelength, velocity_km_s):
    return wavelength * (1.0 + velocity_km_s / _SPEED_OF_LIGHT_KM_S)

# Find the index of the element in array that is closest to value.
# Input arguments:
#   value: target value
#   array: 1-D numpy array to search
# Output arguments:
#   index of the closest element

def find_closest_index_numpy(value, array):
    return np.argmin(np.abs(np.asarray(array) - value))

# ##### 3.  FILE LOADING

# Read a HARPS data file.
# Input arguments:
#   file: filename
# Output arguments:
#   wave: wavelength data (nm)
#   arr1: spectral brightness data

def read_harps_file(file):
    filepath = Path(file)    # Check for full path name
    if not filepath.exists():
        filepath = eso_cache_path / file # Check for file in eso cache
        if not filepath.exists():
            print("File "+str(file)+" not found.")
            return [],[]

    with fits.open(filepath) as fits_file:
        spectral_data = fits_file[1].data
        wave = spectral_data[0][0]
        arr1 = spectral_data[0][1]
    return wave, arr1

# ##### 4.  SEARCH ALGORITHM

# The main optical seti search routine.  Identifies spikes that rise more than
# threshold_multiplier times the local standard deviation above the local median.
# Spikes must have at least min_count and no more than max_count pixels in a row
# above the threshold.
# 
# Input arguments:
#   arr1: spectral intensity
#   min_cont: minimum number of bright pixels in a row to count as a spike
#   max_count: maximum number of bright pixels in a row to count as a spike
#   threshold_multiplier: how many standard deviations above median must pixel be to count as a spike
#   cosmic_ray_threshold: unused
#   window_size: running median/stdev window size
# Output arguments: 
#   hits_start: list of indices of start of identified spikes
#   hits_end: list of indices of end of identified spikes
#   count: length of last hit

def seti_spike_analyzer(arr1, min_count = 4, max_count = 8, threshold_multiplier = 3.5, window_size = 101):
    half_window_size = round((window_size-1)/2)
    continuum = running_median(arr1, window_size)
    count = 0 #reset bright pixel count variable
#    import numpy as np # JCG: Not needed
    flux_threshold = np.array(running_standarddev(arr1, window_size)) * (threshold_multiplier)
#    cosmic_ray_threshold = np.array(running_standarddev(arr1, stwindow)) * (cosmic_ray_threshold) # JCG: Not used
    hits_start = []  # List of starting wavelength indices for spikes
    hits_end = []    # List of ending wavelength indices for spikes
    # List of known airglow lines, this needs to be clearer
    prohibited_wavelengths = list(range(179450, 179650)) + list(range(251750, 251950)) + list(range(210750, 210950)) + list(range(258150, 258350)) + list(range(141550, 141750)) + list(range(141750, 141950)) + list(range(211350, 211550)) # JCG: This is damn ugly
    # Loop over all wavelengths where "continuum" has been calculated.
    for i in range(50,len(continuum) - half_window_size):
            if arr1[i] >= continuum[i - half_window_size] + flux_threshold[i - half_window_size]:  # If pixel is above threshold
                  count += 1                                           # increment count
            else:                                                      # if pixel falls below threshold
                if (count >= min_count) and (count <= max_count):      # if spike isn't too wide or narrow
                    if i not in prohibited_wavelengths:                # if it's not in our list of airglow lines
                        hits_start.append(i-count)                     # Add to the list of spikes found
                        hits_end.append(i)
                count = 0
    print(hits_start, hits_end)
    return hits_start, hits_end, count                                 # Return list of hits found, and length of last hit

# ##### 5.  PLOTTING

# Plot spectral data between index1 and index2.
def original_spectrum_plot(wave, arr1, index1, index2):
    plt.plot(wave[index1:index2], arr1[index1:index2], '.-')

# Plot spectral data with continuum and flux threshold lines.
# Save figures as png.
#
# Input:
#   file: HARPS spectrum filename
#   window_size: window size for running medians and stdev
#   threshold_multiplier: how many standard deviations above median to plot threshold line
#   center_index: Index of center of graph

def spike_plotter(file, window_size = 101, threshold_multiplier = 3.5, center_index = 1000): 
    wave, arr1 = read_harps_file(file)
    continuum = running_median(arr1,window_size) 
    start_index = (center_index - 1000)
    end_index = (center_index + 1000)
    threshold = continuum + np.array(running_standarddev(arr1, window_size)) * threshold_multiplier
    x = wave[start_index:end_index]
    plt.plot(x, arr1[start_index:end_index],'.-', x, continuum[(start_index - 50):(end_index - 50)], x, threshold[(start_index - 50):(end_index - 50)])
    plt.savefig(str(file[48:]) + "zoom_out" + ".png")
    start_index = (center_index - 100)
    end_index = (center_index + 100)
    x = wave[start_index:end_index]
    plt.plot(x, arr1[start_index:end_index],'.-', x, continuum[(start_index - 50):(end_index - 50)], x, threshold[(start_index - 50):(end_index - 50)])
    plt.savefig(str(file[48:]) + "zoom_in" + ".png")

# Fit a Gaussian curve to a spectral line found at hits_start to hits_end, plot both, return the width of the Gaussian fit.
#
# Inputs:
#   file: file name to fit
#   hits_start: starting index of best-guess peak location
#   hits_end: ending index of best-guess peak location
def gaussian_curve_fit(file,hits_start,hits_end):
    wave,arr1 = read_harps_file(file)
    windowpoint1 = hits_start - 100
    windowpoint2 = hits_end + 100
    subtracted = arr1 - np.mean(arr1[windowpoint1:windowpoint2])
    peak_guess = np.max(subtracted[hits_start:hits_end]) #makes a highly "educated guess" for the fitted curve's peak by taking the actual maximum from the subtracted continuum
    mean_guess = np.mean(wave[hits_start:hits_end])
    st_deviation_guess_wide = (wave[hits_end] - wave[hits_start]) * 10
    st_deviation_guess_narrow = (wave[hits_end] - wave[hits_start]) * 2
    spectrum = Spectrum1D(flux=subtracted[windowpoint1:windowpoint2]*u.dimensionless_unscaled, spectral_axis=wave[windowpoint1:windowpoint2]*u.AA)
    g_init = models.Gaussian1D(amplitude=peak_guess*u.dimensionless_unscaled, mean=mean_guess*u.AA, stddev=st_deviation_guess_wide*u.AA)
    g_fit = fit_lines(spectrum, g_init, window=(wave[windowpoint1]*u.AA, wave[windowpoint2]*u.AA))
    standard_deviation = g_fit.stddev.value
    if standard_deviation < 1e-12: # Didn't find peak, try a narrower guess
        alternate_guess = models.Gaussian1D(amplitude=peak_guess*u.dimensionless_unscaled, mean=mean_guess*u.AA, stddev=st_deviation_guess_narrow*u.AA)
        g_fit = fit_lines(spectrum, alternate_guess, window=(wave[windowpoint1]*u.AA, wave[windowpoint2]*u.AA))
        standard_deviation= g_fit.stddev.value
    fwhm = standard_deviation * 2.35
    y_fit = g_fit(wave[windowpoint1:windowpoint2]*u.AA)
    plt.plot(spectrum.spectral_axis, spectrum.flux) 
    plt.plot(wave[windowpoint1:windowpoint2], y_fit)
    print(fwhm) 
    #Cosmic ray hits_start: 4082, green auroral emission line: 179487, 179487
    #GJ551 HITS_STARTS: 13812, 39670 (something's odd), 43403, 45101, 59377, 64514, 67948, 76729, 80165
    plt.show()
    return fwhm