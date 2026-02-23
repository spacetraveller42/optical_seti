#!/usr/bin/env python
"""
Gaussian Generation with FWHM - Complete Standalone Module

This module provides functions to generate Gaussian curves with specified
Full Width Half Maximum (FWHM) and add them to existing arrays.

FEATURES:
- Generate Gaussians with configurable FWHM
- Specify amplitude OR area under curve
- Index-based OR wavelength-based center positioning
- Wavelength interpolation (supports wavelengths between array elements)
- Broadcasting support for multi-dimensional arrays
- Different array lengths (zero-padding/truncation)
- Comprehensive input validation

Author: Generated for optical_seti project
Repository: https://github.com/spacetraveller42/optical_seti

QUICK START:
    import numpy as np
    from gaussian_generator import generate_gaussian, add_gaussian_to_array
    
    # Method 1: Index-based (traditional)
    gaussian = generate_gaussian(fwhm=10.0, amplitude=100, center=50, array_length=100)
    
    # Method 2: Wavelength-based (NEW!)
    wavelengths = np.linspace(400, 700, 100)  # nm
    gaussian = generate_gaussian(fwhm=10.0, amplitude=100, 
                                 wavelength=550,  # nm
                                 array_length=100, wavelengths=wavelengths)
    
    # Add to existing data
    data = np.ones(100) * 500
    result = add_gaussian_to_array(data, fwhm=5.0, amplitude=200, 
                                   wavelength=532, wavelengths=wavelengths)
"""

import numpy as np
import matplotlib.pyplot as plt

# ============================================================================
# CONSTANTS
# ============================================================================

# Conversion factor from FWHM to standard deviation for a Gaussian
# FWHM = 2 * sqrt(2 * ln(2)) * sigma ≈ 2.35482 * sigma
FWHM_TO_SIGMA_FACTOR = 2 * np.sqrt(2 * np.log(2))

# ============================================================================
# CORE FUNCTIONS
# ============================================================================

def generate_gaussian(fwhm, amplitude=None, center=None, array_length=None, area=None,
                     wavelength=None, wavelengths=None):
    """
    Generate a Gaussian curve with specified Full Width Half Maximum (FWHM).
    
    Supports both index-based and wavelength-based center specification.
    Wavelength-based positioning works even for wavelengths between array elements
    using linear interpolation.
    
    Parameters
    ----------
    fwhm : float
        Full Width Half Maximum of the Gaussian (in array index units).
        Must be positive.
    amplitude : float, optional
        Peak amplitude of the Gaussian. Specify either amplitude or area, not both.
    center : float, optional
        Center position as array index (0 to array_length-1).
        Specify either center or wavelength, not both.
    array_length : int
        Length of the output array. Must be positive.
    area : float, optional
        Total area under the Gaussian curve. Specify either amplitude or area, not both.
        When area is specified, the amplitude is automatically calculated.
    wavelength : float, optional
        Target wavelength for Gaussian center (requires wavelengths array).
        Specify either center or wavelength, not both. Can be between array elements.
    wavelengths : array-like, optional
        Array of wavelength values corresponding to array indices.
        Required if wavelength parameter is used. Must be monotonic.
    
    Returns
    -------
    numpy.ndarray
        Array containing the Gaussian curve.
    
    Examples
    --------
    >>> # Index-based center
    >>> g = generate_gaussian(fwhm=10.0, amplitude=100.0, center=50.0, array_length=100)
    
    >>> # Wavelength-based center
    >>> wavelengths = np.linspace(400, 700, 100)
    >>> g = generate_gaussian(fwhm=10.0, amplitude=100.0, wavelength=550.0, 
    ...                       array_length=100, wavelengths=wavelengths)
    
    >>> # Area instead of amplitude
    >>> g = generate_gaussian(fwhm=5.0, area=1000.0, wavelength=525.5,
    ...                       array_length=100, wavelengths=wavelengths)
    """
    # Input validation
    if fwhm <= 0:
        raise ValueError(f"FWHM must be positive, got {fwhm}")
    if array_length <= 0:
        raise ValueError(f"array_length must be positive, got {array_length}")
    
    # Check that either center or wavelength is specified, but not both
    if center is None and wavelength is None:
        raise ValueError("Either 'center' or 'wavelength' must be specified")
    if center is not None and wavelength is not None:
        raise ValueError("Cannot specify both 'center' and 'wavelength' - choose one")
    
    # Handle wavelength-based center specification
    if wavelength is not None:
        if wavelengths is None:
            raise ValueError("'wavelengths' array must be provided when using 'wavelength' parameter")
        
        wavelengths = np.asarray(wavelengths)
        
        # Validate wavelengths array
        if len(wavelengths) != array_length:
            raise ValueError(f"wavelengths array length ({len(wavelengths)}) must match array_length ({array_length})")
        
        # Check if wavelengths array is monotonic
        if not (np.all(np.diff(wavelengths) > 0) or np.all(np.diff(wavelengths) < 0)):
            raise ValueError("wavelengths array must be monotonically increasing or decreasing")
        
        # Check if target wavelength is within range
        wave_min, wave_max = np.min(wavelengths), np.max(wavelengths)
        if wavelength < wave_min or wavelength > wave_max:
            raise ValueError(f"wavelength {wavelength} is outside wavelengths range [{wave_min}, {wave_max}]")
        
        # Convert wavelength to fractional array index using interpolation
        # This allows wavelengths between array elements
        indices = np.arange(array_length)
        center = float(np.interp(wavelength, wavelengths, indices))
    
    # Validate center position
    if center < 0 or center >= array_length:
        raise ValueError(f"center must be within [0, {array_length}), got {center}")
    
    # Warn if center looks like it might be a wavelength value instead of an index
    if center < 10 and wavelengths is not None:
        wavelengths_arr = np.asarray(wavelengths)
        wave_min, wave_max = np.min(wavelengths_arr), np.max(wavelengths_arr)
        if wave_min <= center <= wave_max:
            import warnings
            warnings.warn(
                f"\n⚠️  WARNING: center={center} looks like a WAVELENGTH value, not an array INDEX!\n"
                f"\nYou specified: center={center}\n"
                f"But 'center' expects an array INDEX (0 to {array_length-1}),\n"
                f"and {center} is within your wavelengths range [{wave_min:.3f}, {wave_max:.3f}].\n"
                f"\nDid you mean to use:\n"
                f"  wavelength={center}, wavelengths=<your_wavelength_array>\n",
                UserWarning, stacklevel=2
            )
    
    # Check that either amplitude or area is specified, but not both
    if amplitude is None and area is None:
        raise ValueError("Either 'amplitude' or 'area' must be specified")
    if amplitude is not None and area is not None:
        raise ValueError("Cannot specify both 'amplitude' and 'area' - choose one")
    
    # Convert FWHM to standard deviation
    sigma = fwhm / FWHM_TO_SIGMA_FACTOR
    
    # If area is specified, calculate the required amplitude
    # For a Gaussian: Area = amplitude * sigma * sqrt(2*pi)
    if area is not None:
        if area <= 0:
            raise ValueError(f"area must be positive, got {area}")
        amplitude = area / (sigma * np.sqrt(2 * np.pi))
    
    # Create x-axis array
    x = np.arange(array_length)
    
    # Generate Gaussian
    gaussian_array = amplitude * np.exp(-((x - center) ** 2) / (2 * sigma ** 2))
    
    return gaussian_array


def add_gaussian_to_array(data, fwhm, amplitude=None, center=None, array_length=None, axis=-1, area=None,
                         wavelength=None, wavelengths=None):
    """
    Add a Gaussian curve to an existing array with broadcasting support.
    
    Supports multi-dimensional arrays, different array lengths, and both
    index-based and wavelength-based center specification.
    
    Parameters
    ----------
    data : numpy.ndarray
        Array to which the Gaussian will be added (1D, 2D, or higher).
    fwhm : float
        Full Width Half Maximum of the Gaussian.
    amplitude : float, optional
        Peak amplitude. Specify either amplitude or area, not both.
    center : float, optional
        Center position as array index. Specify either center or wavelength.
    array_length : int, optional
        Length of Gaussian to generate. If None, uses data length along axis.
    axis : int, optional
        Axis for multi-dimensional arrays. Default is -1 (last axis).
    area : float, optional
        Total area under curve. Specify either amplitude or area, not both.
    wavelength : float, optional
        Target wavelength for center (requires wavelengths array).
    wavelengths : array-like, optional
        Array of wavelength values. Required if wavelength parameter used.
    
    Returns
    -------
    numpy.ndarray
        Array with the Gaussian added.
    
    Examples
    --------
    >>> # 1D with index-based center
    >>> data = np.ones(100) * 500
    >>> result = add_gaussian_to_array(data, fwhm=10.0, amplitude=100.0, center=50.0)
    
    >>> # 1D with wavelength-based center
    >>> wavelengths = np.linspace(400, 700, 100)
    >>> result = add_gaussian_to_array(data, fwhm=10.0, amplitude=100.0, 
    ...                                wavelength=550.0, wavelengths=wavelengths)
    
    >>> # 2D array
    >>> data_2d = np.ones((50, 100)) * 500
    >>> result = add_gaussian_to_array(data_2d, fwhm=10.0, amplitude=100.0, 
    ...                                center=50.0, axis=1)
    """
    data = np.asarray(data)
    
    if data.size == 0:
        raise ValueError("Input data array is empty")
    
    # Handle multi-dimensional arrays
    if data.ndim > 1:
        # Normalize axis to positive index
        if axis < 0:
            axis = data.ndim + axis
        
        if axis < 0 or axis >= data.ndim:
            raise ValueError(f"axis {axis} is out of bounds for array of dimension {data.ndim}")
        
        if array_length is None:
            array_length = data.shape[axis]
        
        if array_length <= 0:
            raise ValueError(f"array_length must be positive, got {array_length}")
        
        if center is not None and (center < 0 or center >= array_length):
            raise ValueError(f"center {center} is out of bounds for array_length {array_length}")
        
        # Generate Gaussian
        gaussian = generate_gaussian(fwhm, amplitude=amplitude, center=center, 
                                    array_length=array_length, area=area,
                                    wavelength=wavelength, wavelengths=wavelengths)
        
        # Handle different lengths
        data_len = data.shape[axis]
        gauss_len = len(gaussian)
        
        if gauss_len < data_len:
            pad_gaussian = np.zeros(data_len)
            pad_gaussian[:gauss_len] = gaussian
            gaussian = pad_gaussian
        elif gauss_len > data_len:
            gaussian = gaussian[:data_len]
        
        # Reshape for broadcasting
        broadcast_shape = [1] * data.ndim
        broadcast_shape[axis] = len(gaussian)
        gaussian = gaussian.reshape(broadcast_shape)
        
        result = data + gaussian
    else:
        # 1D array
        if array_length is None:
            array_length = len(data)
        
        if array_length <= 0:
            raise ValueError(f"array_length must be positive, got {array_length}")
        
        if center is not None and (center < 0 or center >= array_length):
            raise ValueError(f"center {center} is out of bounds for array_length {array_length}")
        
        # Generate Gaussian
        gaussian = generate_gaussian(fwhm, amplitude=amplitude, center=center, 
                                    array_length=array_length, area=area,
                                    wavelength=wavelength, wavelengths=wavelengths)
        
        # Handle different lengths
        data_len = len(data)
        gauss_len = len(gaussian)
        
        if gauss_len == data_len:
            result = data + gaussian
        elif gauss_len < data_len:
            padded_gaussian = np.zeros(data_len)
            padded_gaussian[:gauss_len] = gaussian
            result = data + padded_gaussian
        else:
            result = data + gaussian[:data_len]
    
    return result


def plot_gaussian_comparison(data, fwhm, amplitude, center, title="Gaussian Addition Comparison"):
    """
    Plot original data and data with added Gaussian side by side.
    
    Parameters
    ----------
    data : numpy.ndarray
        Original data array.
    fwhm : float
        Full Width Half Maximum of the Gaussian.
    amplitude : float
        Peak amplitude of the Gaussian.
    center : float
        Center position of the Gaussian.
    title : str, optional
        Title for the plot.
    
    Returns
    -------
    tuple
        (figure, axes) - matplotlib figure and axes objects.
    """
    gaussian = generate_gaussian(fwhm, amplitude, center, len(data))
    result = data + gaussian
    
    fig, axes = plt.subplots(3, 1, figsize=(12, 10))
    
    # Plot 1: Pure Gaussian
    axes[0].plot(gaussian, 'b-', linewidth=2)
    axes[0].axhline(y=amplitude/2, color='r', linestyle='--', alpha=0.5, 
                    label=f'Half Max (FWHM={fwhm})')
    axes[0].axvline(x=center, color='g', linestyle=':', alpha=0.5, label='Center')
    axes[0].set_title(f'Gaussian (FWHM={fwhm}, Amplitude={amplitude})', 
                      fontsize=12, fontweight='bold')
    axes[0].set_xlabel('Array Index')
    axes[0].set_ylabel('Value')
    axes[0].legend()
    axes[0].grid(True, alpha=0.3)
    
    # Plot 2: Original data
    axes[1].plot(data, 'gray', linewidth=1.5)
    axes[1].set_title('Original Data', fontsize=12, fontweight='bold')
    axes[1].set_xlabel('Array Index')
    axes[1].set_ylabel('Value')
    axes[1].grid(True, alpha=0.3)
    
    # Plot 3: Combined
    axes[2].plot(data, 'gray', alpha=0.4, linewidth=1, label='Original')
    axes[2].plot(result, 'b-', linewidth=2, label='With Gaussian')
    axes[2].axvline(x=center, color='r', linestyle=':', alpha=0.5)
    axes[2].set_title(title, fontsize=12, fontweight='bold')
    axes[2].set_xlabel('Array Index')
    axes[2].set_ylabel('Value')
    axes[2].legend()
    axes[2].grid(True, alpha=0.3)
    
    plt.tight_layout()
    return fig, axes


# ============================================================================
# DEMONSTRATION AND TESTING
# ============================================================================

def run_demo():
    """Run comprehensive demonstration of functionality."""
    print("=" * 70)
    print("GAUSSIAN GENERATION WITH FWHM - DEMONSTRATION")
    print("=" * 70)
    
    # Demo 1: Simple generation
    print("\n1. Generate a simple Gaussian:")
    gaussian = generate_gaussian(fwhm=15.0, amplitude=200.0, center=100.0, array_length=200)
    print(f"   ✓ Generated Gaussian with FWHM=15.0, amplitude=200.0")
    print(f"   ✓ Array shape: {gaussian.shape}")
    print(f"   ✓ Peak value: {np.max(gaussian):.2f}")
    
    # Demo 2: Add to data
    print("\n2. Add Gaussian to existing spectral data:")
    spectral_data = np.random.randn(500) * 10 + 1000
    laser_signal = add_gaussian_to_array(spectral_data, fwhm=5.0, amplitude=500, center=250)
    print(f"   ✓ Original data range: [{np.min(spectral_data):.2f}, {np.max(spectral_data):.2f}]")
    print(f"   ✓ With signal range: [{np.min(laser_signal):.2f}, {np.max(laser_signal):.2f}]")
    
    # Demo 3: Multiple Gaussians
    print("\n3. Add multiple Gaussians:")
    baseline = np.random.randn(1000) * 5 + 500
    result = add_gaussian_to_array(baseline, fwhm=4.0, amplitude=300, center=200)
    result = add_gaussian_to_array(result, fwhm=12.0, amplitude=200, center=500)
    result = add_gaussian_to_array(result, fwhm=30.0, amplitude=150, center=800)
    print(f"   ✓ Added 3 Gaussians with different FWHM values")
    
    print("\n" + "=" * 70)
    print("DEMONSTRATION COMPLETE")
    print("=" * 70)


def run_tests():
    """Run basic tests to verify functionality."""
    print("\n" + "=" * 70)
    print("RUNNING TESTS")
    print("=" * 70)
    
    all_passed = True
    
    # Test 1: Basic generation
    print("\nTest 1: Basic Gaussian generation")
    try:
        gaussian = generate_gaussian(fwhm=10.0, amplitude=100.0, center=50.0, array_length=100)
        assert np.argmax(gaussian) == 50
        assert abs(np.max(gaussian) - 100.0) < 0.01
        print("   ✓ PASSED")
    except Exception as e:
        print(f"   ✗ FAILED: {e}")
        all_passed = False
    
    # Test 2: Array addition
    print("\nTest 2: Array addition")
    try:
        baseline = np.ones(100) * 50.0
        result = add_gaussian_to_array(baseline, fwhm=10.0, amplitude=100.0, center=50.0)
        assert abs(np.max(result) - 150.0) < 0.01
        print("   ✓ PASSED")
    except Exception as e:
        print(f"   ✗ FAILED: {e}")
        all_passed = False
    
    # Test 3: Input validation
    print("\nTest 3: Input validation")
    try:
        try:
            generate_gaussian(fwhm=-5.0, amplitude=100.0, center=50.0, array_length=100)
            all_passed = False
        except ValueError:
            print("   ✓ Correctly rejected negative FWHM")
        
        try:
            generate_gaussian(fwhm=10.0, amplitude=100.0, center=150.0, array_length=100)
            all_passed = False
        except ValueError:
            print("   ✓ Correctly rejected out-of-bounds center")
    except Exception as e:
        print(f"   ✗ FAILED: {e}")
        all_passed = False
    
    print("\n" + "=" * 70)
    if all_passed:
        print("ALL TESTS PASSED ✓")
    else:
        print("SOME TESTS FAILED ✗")
    print("=" * 70)
    
    return all_passed


# ============================================================================
# MAIN EXECUTION
# ============================================================================

if __name__ == "__main__":
    """When run as a script, execute tests and demonstration."""
    import sys
    
    tests_passed = run_tests()
    
    if tests_passed:
        print("\n")
        run_demo()
    else:
        print("\nSkipping demo due to test failures.")
        sys.exit(1)
