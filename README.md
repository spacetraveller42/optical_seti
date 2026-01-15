# HARPS Optical SETI

This project is an attempt to find artificial narrow laser signals in the spectra of stars, using data from ESO's [HARPS](https://www.eso.org/sci/facilities/lasilla/instruments/harps.html) spectrometer.  It allows the user to download spectra from the [ESO data archive](http://archive.eso.org/eso/eso_archive_main.html) and search for narrowband spikes in the spectra.  It attempts to reject known terrestrial airglow lines, emission lines from active stars, and other false positives.  The user can compare the processed spectral data to the raw CCD image data, which can help to identify false positives from cosmic rays.

It also includes other utility functions, such as fitting a simple Gaussian line profile to a line to estimate its width, or plotting spectral data.

## Dependencies

This project runs in Python 3.  To use it, the user must install the following packages using `pip`:
- `numpy`
- `astropy`
- `astroquery`
- `matplotlib`
- `specutils`

You must also have an [account set up](https://www.eso.org/sso/login?service=https%3a%2f%2farchive.eso.org%2fwdb%2fforms%2fcas%2feso_archive_main.html) at the ESO Archive, and have configured `astroquery` to [store your password](https://astroquery.readthedocs.io/en/latest/eso/eso.html#automatic-password) in its keyring.

## Modules and Primary Functions

### `optical_seti_functions.py`
 This file contains the main search algorithm and support routines.

- `seti_spike_analyzer()`
  
  The main optical seti search routine.  Identifies spikes that rise more than `threshold_multiplier` times the local standard deviation above the local median.  Spikes must have at least `min_count` and no more than `max_count` pixels in a row above the threshold.
  
        hits_start, hits_end, count = seti_spike_analyzer(arr1, min_count = 4, max_count = 8, threshold_multiplier = 3.5, window_size = 101)

- `spike_plotter()`
  
  Plot spectral data with continuum and flux threshold lines.  Save figures as png.
  
        spike_plotter(file, window_size = 101, threshold_multiplier = 3.5, center_index = 1000)

- `gaussian_curve_fit()` 
  
  Fit a Gaussian curve to a spectral line found at `hits_start` to `hits_end`, plot both, return the width of the Gaussian fit.
  
        fwhm = gaussian_curve_fit("filename.fits",hits_start,hits_end)

- `generate_gaussian()`
  
  Generate a Gaussian curve with specified Full Width Half Maximum (FWHM). You can specify either the peak amplitude OR the total area under the curve. This is useful for creating synthetic test data or adding artificial signals to spectra.
  
        # Using amplitude (traditional method)
        gaussian_array = generate_gaussian(fwhm=10.0, amplitude=100.0, center=50.0, array_length=100)
        
        # Using area under the curve (new feature)
        gaussian_array = generate_gaussian(fwhm=10.0, area=1000.0, center=50.0, array_length=100)

- `add_gaussian_to_array()`
  
  Add a Gaussian curve to an existing data array with broadcasting support for multi-dimensional arrays. The Gaussian can have a different length than the data array (will be zero-padded or truncated as needed). You can specify either amplitude or area.
  
        # 1D array with amplitude (default behavior)
        data_with_gaussian = add_gaussian_to_array(data, fwhm=10.0, amplitude=100.0, center=50.0)
        
        # 1D array with area (controls total signal strength)
        data_with_gaussian = add_gaussian_to_array(data, fwhm=10.0, area=1000.0, center=50.0)
        
        # Different Gaussian length (shorter Gaussians are zero-padded, longer ones truncated)
        data_with_gaussian = add_gaussian_to_array(data, fwhm=10.0, amplitude=100.0, center=50.0, array_length=150)
        
        # 2D array - add Gaussian along specific axis (broadcasting)
        data_2d = np.ones((50, 100)) * 500
        result_2d = add_gaussian_to_array(data_2d, fwhm=10.0, amplitude=100.0, center=50.0, axis=1)  # broadcasts across rows

### `seti_catalog_functions.py`
These functions operate on entire catalogs of stars (CSV/TSV text files).

- `do_search()`
  
  Search for spikes in all stars in a catalog by calling `optical_seti_functions.seti_spike_analyzer()` on each one.  
  
  `categorize`: If True, two results files will be created, one for stars with 3 spikes or fewer, one with 4 or more (possible active stars).  

  `predownloader_format`: If `True`, assume the input star catalog is a CSV file in the format created by `predowloader()` below.  If `False`, assume it's in the TSV format used by OSETI_targets.txt.
    
        do_search(withlogin="your_eso_login",withstarlist="starlist.txt",withresults="results.txt",categorize=False,predownloader_format=False)
- `predownloader()`
  
   Predownload one spectrum from each of the stars in the catalog.

        predownload(withlogin="your_eso_login",withstarlist="starlist.txt",withresults="results.txt")

- `bulk_predownloader()`
  
  Predownload one or more spectra from each of the stars in the catalog. Group up files to download to reduce the number of requests we make to ESO servers.

        bulk_predownload(withlogin="your_eso_login",withstarlist="starlist.txt",withresults="results.txt",obs_per_star=1)

### `harpscompare.py`

  Download HARPS spectral and CCD data and plot it side-by-side.  This can allow the user to identify cosmic rays as false positives.

- `compare_spec_to_raw()`
  Plot spectrum and CCD image file on the same figure. `countmin` and `countmax` are optional arguments describing the scale of the spectral graph and the color range of the CCD image.

        compare_spec_to_raw(specfilename,rawfilename,lamb,countmin,countmax)

### `test/`

Test scripts for major functions.  See the code here for usage examples and verification.

### `OSETI_targets.txt`

A list of 2800+ stars that have been observed by HARPS.  (Reference source URL?)

### `harps_spectralpositioning.txt`

Data file describing the position of the spectral orders of HARPS's echelle spectrograph on its CCD chips.  Used by `harpscompare.py`.  (Reference source URL?)