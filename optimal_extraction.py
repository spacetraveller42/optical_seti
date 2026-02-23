# #####################################
# #                                   #
# #    OPTIMAL_EXTRACTION.PY          #
# #                                   #
# #####################################
#
# Standalone module for optimal extraction of HARPS spectra using PyReduce.
#
# This module can be copied into any project as a self-contained file.
# It filters out anything that is not actual light entering the instrument's
# optics, including bias, dark current, bad/hot pixels, scattered light,
# and cosmic rays.
#
# Dependencies (install via pip):
#   numpy, astropy, scipy, pyreduce-astro
#
# Usage:
#   from optimal_extraction import optimal_extraction_harps
#   spectrum, uncertainties, slitfunction = optimal_extraction_harps(raw_image)
#
# Authors: Benjamin Fields, Jason Goodman
# Repo: https://github.com/spacetraveller42/optical_seti

import numpy as np
from astropy.io import fits
from pathlib import Path


def optimal_extraction_harps(raw_image, order_traces=None, extraction_width=5,
                             column_range=None, bias=None, dark=None,
                             bad_pixel_mask=None, scatter=None,
                             gain=1.0, readnoise=0.0,
                             lambda_sf=0.1, lambda_sp=0, osample=1,
                             swath_width=None, maxiter=20, plot=False):
    """Perform optimal extraction on a raw HARPS CCD image using PyReduce.

    This function extracts spectra while filtering out anything that is not
    actual light entering the instrument's optics.  It handles:

    - Bias: electronic offset from the readout electronics
    - Dark current: thermal electrons that accumulate during exposure
    - Bad/hot pixels: detector defects that produce spurious signal
    - Scattered light: stray light inside the spectrograph
    - Cosmic rays: high-energy particle hits affecting one or two pixels

    After calibration, the Horne (1986) optimal extraction algorithm uses
    the known spatial profile (slit function) of each spectral order to
    weight pixels.  Any remaining single-pixel artifacts that do not match
    the expected spatial profile are rejected.  Real light spread by the
    telescope optics matches the slit function and is preserved.

    Parameters
    ----------
    raw_image : array-like or str or Path
        2D numpy array of the raw CCD image, or a FITS filename to load.
    order_traces : array-like, optional
        Polynomial coefficients for order tracing, shape ``[nord, degree]``.
        Each row contains polynomial coefficients (highest degree first)
        that describe the y-position of the order center as a function of
        x-pixel.  If ``None``, orders are traced automatically from the image.
    extraction_width : scalar or array-like, optional
        Extraction width in pixels (below, above), shape ``[nord, 2]``,
        or a scalar applied symmetrically to all orders.  Default is 5.
    column_range : array-like, optional
        Array ``[nord, 2]`` specifying the valid column range for each order.
        ``None`` uses the full image width.
    bias : array-like or scalar, optional
        Bias frame (or scalar) to subtract.  Removes the electronic offset
        added by the CCD readout electronics.
    dark : array-like, optional
        Dark frame to subtract.  Removes thermally generated electrons
        that are not from photons.
    bad_pixel_mask : array-like of bool, optional
        2D boolean array where ``True`` marks bad/hot pixels to exclude.
        Masked pixels are replaced with the median of their neighbors.
    scatter : array-like, optional
        Background scattered light model, array ``[nord, 4, ncol]`` of
        polynomial coefficients for inter-order background.
    gain : float, optional
        CCD gain in electrons/ADU.  Default is 1.0.
    readnoise : float, optional
        CCD read noise in electrons.  Default is 0.0.
    lambda_sf : float, optional
        Slit function smoothing parameter.  Default is 0.1.
    lambda_sp : float, optional
        Spectrum smoothing parameter.  Default is 0.
    osample : int, optional
        Oversampling factor for the slit function.  Default is 1.
    swath_width : int, optional
        Width of swaths for extraction.  ``None`` auto-determines the width.
    maxiter : int, optional
        Maximum iterations for the slit function decomposition.  Default is 20.
    plot : bool, optional
        Whether to plot extraction progress.  Default is ``False``.

    Returns
    -------
    spectrum : ndarray, shape ``[nord, ncol]``
        Extracted spectra with non-optical artifacts filtered out.
    uncertainties : ndarray, shape ``[nord, ncol]``
        Uncertainties on the extracted spectra.
    slitfunction : list of ndarray
        Recovered slit functions for each order.
    """
    from pyreduce.extract import optimal_extraction, fix_extraction_width, fix_column_range
    from pyreduce.trace_orders import mark_orders
    from scipy.ndimage import median_filter

    # Load image from file if a filename is given
    if isinstance(raw_image, (str, Path)):
        with fits.open(raw_image) as hdul:
            img_data = hdul[0].data
            if img_data is None:
                # Try the first extension with image data
                for ext in hdul:
                    if ext.data is not None and ext.data.ndim == 2:
                        img_data = ext.data
                        break
            if img_data is None:
                raise ValueError("No 2D image data found in FITS file")
            raw_image = img_data.astype(float)
    raw_image = np.asarray(raw_image, dtype=float)

    if raw_image.ndim != 2:
        raise ValueError("raw_image must be a 2D array, got shape {}".format(raw_image.shape))

    # --- Calibration: remove non-optical artifacts before extraction ---
    # Work on a copy so the caller's original array is not modified.
    if bias is not None or dark is not None or bad_pixel_mask is not None:
        raw_image = raw_image.copy()

    # Subtract bias (electronic readout offset)
    if bias is not None:
        raw_image -= np.asarray(bias, dtype=float)

    # Subtract dark current (thermal electrons)
    if dark is not None:
        raw_image -= np.asarray(dark, dtype=float)

    # Replace bad/hot pixels with local median so they don't corrupt extraction.
    # Only the flagged pixel positions are overwritten; real signal is untouched.
    if bad_pixel_mask is not None:
        bad_pixel_mask = np.asarray(bad_pixel_mask, dtype=bool)
        filtered = median_filter(raw_image, size=5)
        raw_image[bad_pixel_mask] = filtered[bad_pixel_mask]

    nrow, ncol = raw_image.shape

    # Trace orders automatically if not provided
    if order_traces is None:
        order_traces = mark_orders(raw_image, manual=False, plot=plot)

    order_traces = np.asarray(order_traces)
    nord = len(order_traces)

    # Set up column range
    if column_range is None:
        column_range = np.tile([0, ncol], (nord, 1))
    column_range = np.asarray(column_range, dtype=int)

    # Set up extraction width
    if np.isscalar(extraction_width):
        extraction_width = np.full((nord, 2), extraction_width, dtype=float)
    extraction_width = np.asarray(extraction_width, dtype=float)

    # Fix extraction width and column range to stay within image bounds
    extraction_width = fix_extraction_width(extraction_width, order_traces, column_range, ncol)
    column_range, order_traces = fix_column_range(column_range, order_traces, extraction_width, nrow, ncol)

    # No tilt or shear correction by default (straight orders).
    # optimal_extraction expects per-order arrays (list of arrays or Nones).
    tilt = [None] * nord
    shear = [None] * nord

    # Perform optimal extraction.  The slit-function decomposition rejects any
    # remaining pixel-level artifacts (cosmic rays, residual hot pixels, etc.)
    # that do not match the expected spatial profile, while preserving real light.
    spectrum, slitfunction, uncertainties = optimal_extraction(
        raw_image, order_traces, extraction_width, column_range,
        tilt=tilt, shear=shear,
        scatter=scatter,
        gain=gain, readnoise=readnoise,
        lambda_sf=lambda_sf, lambda_sp=lambda_sp,
        osample=osample, swath_width=swath_width,
        maxiter=maxiter, plot=plot
    )

    return spectrum, uncertainties, slitfunction
