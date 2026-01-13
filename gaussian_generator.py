#!/usr/bin/env python
"""
Gaussian Generation with FWHM - Standalone Module

This module provides functions to generate Gaussian curves with specified
Full Width Half Maximum (FWHM) and add them to existing arrays for visualization
and analysis.

Author: Generated for optical_seti project
Repository: https://github.com/spacetraveller42/optical_seti

Usage Example:
    import numpy as np
    import matplotlib.pyplot as plt
    from gaussian_generator import generate_gaussian, add_gaussian_to_array
    
    # Generate a Gaussian
    gaussian = generate_gaussian(fwhm=10.0, amplitude=100, center=50, array_length=100)
    
    # Add to existing data
    data = np.random.randn(100) * 5 + 500
    result = data + gaussian
    
    # Plot
    plt.plot(result)
    plt.show()
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

def generate_gaussian(fwhm, amplitude, center, array_length):
    """
    Generate a Gaussian curve with specified Full Width Half Maximum (FWHM).
    
    The Gaussian is generated on an x-axis that spans indices 0 to array_length-1.
    Returns a numpy array that can be directly added to existing arrays.
    
    Parameters
    ----------
    fwhm : float
        Full Width Half Maximum of the Gaussian (in array index units).
        Must be positive.
    amplitude : float
        Peak amplitude of the Gaussian.
    center : float
        Center position of the Gaussian (in array index units).
        Must be within [0, array_length).
    array_length : int
        Length of the output array. Must be positive.
    
    Returns
    -------
    numpy.ndarray
        Array containing the Gaussian curve.
    
    Raises
    ------
    ValueError
        If fwhm <= 0, array_length <= 0, or center is out of bounds.
    
    Examples
    --------
    >>> gaussian = generate_gaussian(fwhm=10.0, amplitude=100.0, center=50.0, array_length=100)
    >>> print(gaussian.shape)
    (100,)
    >>> print(np.max(gaussian))
    100.0
    
    >>> # Add to existing data
    >>> data = np.ones(100) * 500
    >>> result = data + gaussian
    >>> print(np.max(result))
    600.0
    """
    # Input validation
    if fwhm <= 0:
        raise ValueError(f"FWHM must be positive, got {fwhm}")
    if array_length <= 0:
        raise ValueError(f"array_length must be positive, got {array_length}")
    if center < 0 or center >= array_length:
        raise ValueError(f"center must be within [0, {array_length}), got {center}")
    
    # Convert FWHM to standard deviation: FWHM = 2 * sqrt(2 * ln(2)) * sigma
    sigma = fwhm / FWHM_TO_SIGMA_FACTOR
    
    # Create x-axis array
    x = np.arange(array_length)
    
    # Generate Gaussian: amplitude * exp(-(x-center)^2 / (2*sigma^2))
    gaussian_array = amplitude * np.exp(-((x - center) ** 2) / (2 * sigma ** 2))
    
    return gaussian_array


def add_gaussian_to_array(data, fwhm, amplitude, center, array_length=None):
    """
    Add a Gaussian curve to an existing array.
    
    This function generates a Gaussian with the specified FWHM and adds it
    to the input array. The Gaussian can have a different length than the data
    array - if shorter, it will be zero-padded; if longer, it will be truncated.
    
    Parameters
    ----------
    data : numpy.ndarray
        Numpy array to which the Gaussian will be added.
    fwhm : float
        Full Width Half Maximum of the Gaussian (in array index units).
        Must be positive.
    amplitude : float
        Peak amplitude of the Gaussian.
    center : float
        Center position of the Gaussian (in array index units).
    array_length : int, optional
        Length of Gaussian to generate. If None (default), uses len(data).
        If different from data length, Gaussian will be zero-padded or truncated.
    
    Returns
    -------
    numpy.ndarray
        Array with the Gaussian added to the input data.
    
    Raises
    ------
    ValueError
        If fwhm <= 0 or center is out of bounds.
    
    Examples
    --------
    >>> # Same length (default behavior)
    >>> data = np.ones(100) * 500
    >>> result = add_gaussian_to_array(data, fwhm=10.0, amplitude=100.0, center=50.0)
    >>> print(np.max(result))
    600.0
    
    >>> # Add Gaussian of different length (will be zero-padded/truncated)
    >>> data = np.ones(200) * 500
    >>> result = add_gaussian_to_array(data, fwhm=10.0, amplitude=100.0, center=50.0, array_length=100)
    """
    # If array_length not specified, match the data length
    if array_length is None:
        array_length = len(data)
    
    # Generate Gaussian with specified length
    gaussian = generate_gaussian(fwhm, amplitude, center, array_length)
    
    # Handle different lengths
    data_len = len(data)
    gauss_len = len(gaussian)
    
    if gauss_len == data_len:
        # Same length - direct addition
        result = data + gaussian
    elif gauss_len < data_len:
        # Gaussian is shorter - zero-pad it
        padded_gaussian = np.zeros(data_len)
        padded_gaussian[:gauss_len] = gaussian
        result = data + padded_gaussian
    else:
        # Gaussian is longer - truncate it
        result = data + gaussian[:data_len]
    
    return result


def plot_gaussian_comparison(data, fwhm, amplitude, center, title="Gaussian Addition Comparison"):
    """
    Plot original data and data with added Gaussian side by side.
    
    This is a utility function for visualizing the effect of adding a Gaussian
    to existing data.
    
    Parameters
    ----------
    data : numpy.ndarray
        Original data array.
    fwhm : float
        Full Width Half Maximum of the Gaussian to add.
    amplitude : float
        Peak amplitude of the Gaussian.
    center : float
        Center position of the Gaussian.
    title : str, optional
        Title for the plot (default: "Gaussian Addition Comparison").
    
    Returns
    -------
    tuple
        (figure, axes) - matplotlib figure and axes objects.
    
    Examples
    --------
    >>> data = np.random.randn(200) * 10 + 500
    >>> fig, axes = plot_gaussian_comparison(data, fwhm=12.0, amplitude=200.0, center=100.0)
    >>> plt.show()
    """
    # Generate Gaussian
    gaussian = generate_gaussian(fwhm, amplitude, center, len(data))
    result = data + gaussian
    
    # Create figure
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
    """
    Run a comprehensive demonstration of the Gaussian generation functionality.
    
    This function demonstrates:
    1. Simple Gaussian generation
    2. Adding Gaussian to existing data
    3. Multiple Gaussians with different FWHM values
    4. Visualization of results
    """
    print("=" * 70)
    print("GAUSSIAN GENERATION WITH FWHM - DEMONSTRATION")
    print("=" * 70)
    
    # Demo 1: Simple Gaussian generation
    print("\n1. Generate a simple Gaussian:")
    fwhm = 15.0
    amplitude = 200.0
    center = 100.0
    array_length = 200
    
    gaussian = generate_gaussian(fwhm, amplitude, center, array_length)
    print(f"   ✓ Generated Gaussian with FWHM={fwhm}, amplitude={amplitude}")
    print(f"   ✓ Array shape: {gaussian.shape}")
    print(f"   ✓ Peak value: {np.max(gaussian):.2f}")
    print(f"   ✓ Peak location: {np.argmax(gaussian)}")
    
    # Demo 2: Add to existing data
    print("\n2. Add Gaussian to existing spectral data:")
    spectral_data = np.random.randn(500) * 10 + 1000
    
    laser_signal = add_gaussian_to_array(
        spectral_data, 
        fwhm=5.0,
        amplitude=500,
        center=250
    )
    
    print(f"   ✓ Original data range: [{np.min(spectral_data):.2f}, {np.max(spectral_data):.2f}]")
    print(f"   ✓ With signal range: [{np.min(laser_signal):.2f}, {np.max(laser_signal):.2f}]")
    print(f"   ✓ Successfully added narrow Gaussian (FWHM=5.0)")
    
    # Demo 3: Multiple Gaussians
    print("\n3. Add multiple Gaussians with different FWHM values:")
    baseline = np.random.randn(1000) * 5 + 500
    
    result = add_gaussian_to_array(baseline, fwhm=4.0, amplitude=300, center=200)
    print(f"   ✓ Added narrow signal (FWHM=4.0) at pixel 200")
    
    result = add_gaussian_to_array(result, fwhm=12.0, amplitude=200, center=500)
    print(f"   ✓ Added medium signal (FWHM=12.0) at pixel 500")
    
    result = add_gaussian_to_array(result, fwhm=30.0, amplitude=150, center=800)
    print(f"   ✓ Added wide signal (FWHM=30.0) at pixel 800")
    
    # Visualization
    print("\n4. Creating visualizations...")
    
    # Figure 1: Simple comparison
    fig1, axes1 = plot_gaussian_comparison(
        spectral_data[:200], 
        fwhm=15.0, 
        amplitude=400.0, 
        center=100.0,
        title="Single Gaussian Addition"
    )
    plt.figure(fig1.number)
    plt.savefig('gaussian_demo_single.png', dpi=150, bbox_inches='tight')
    print("   ✓ Saved 'gaussian_demo_single.png'")
    
    # Figure 2: Multiple Gaussians
    fig2, axes2 = plt.subplots(figsize=(14, 6))
    axes2.plot(baseline, 'lightgray', alpha=0.5, linewidth=0.5, label='Baseline')
    axes2.plot(result, 'b-', linewidth=1.5, label='With 3 Gaussians')
    axes2.axvline(x=200, color='r', linestyle=':', alpha=0.5)
    axes2.axvline(x=500, color='r', linestyle=':', alpha=0.5)
    axes2.axvline(x=800, color='r', linestyle=':', alpha=0.5)
    axes2.set_title('Multiple Gaussians (FWHM: 4, 12, 30)', fontsize=14, fontweight='bold')
    axes2.set_xlabel('Pixel Index')
    axes2.set_ylabel('Flux')
    axes2.legend()
    axes2.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig('gaussian_demo_multiple.png', dpi=150, bbox_inches='tight')
    print("   ✓ Saved 'gaussian_demo_multiple.png'")
    
    print("\n" + "=" * 70)
    print("DEMONSTRATION COMPLETE")
    print("=" * 70)
    print("\n✓ You CAN set the FWHM directly")
    print("✓ You GET a numpy array from generate_gaussian()")
    print("✓ You CAN add it directly to existing arrays")
    print("✓ Gaussians are clearly visible in plots")
    print("=" * 70)
    
    plt.show()


def run_tests():
    """
    Run basic tests to verify functionality.
    
    Returns
    -------
    bool
        True if all tests pass, False otherwise.
    """
    print("\n" + "=" * 70)
    print("RUNNING TESTS")
    print("=" * 70)
    
    all_passed = True
    
    # Test 1: Basic generation
    print("\nTest 1: Basic Gaussian generation")
    try:
        gaussian = generate_gaussian(fwhm=10.0, amplitude=100.0, center=50.0, array_length=100)
        peak_index = np.argmax(gaussian)
        peak_value = np.max(gaussian)
        
        assert peak_index == 50, f"Peak at {peak_index}, expected 50"
        assert abs(peak_value - 100.0) < 0.01, f"Peak value {peak_value}, expected 100.0"
        print("   ✓ PASSED")
    except Exception as e:
        print(f"   ✗ FAILED: {e}")
        all_passed = False
    
    # Test 2: Array addition
    print("\nTest 2: Array addition")
    try:
        baseline = np.ones(100) * 50.0
        result = add_gaussian_to_array(baseline, fwhm=10.0, amplitude=100.0, center=50.0)
        expected_peak = 150.0
        actual_peak = np.max(result)
        
        assert abs(actual_peak - expected_peak) < 0.01, \
            f"Peak {actual_peak}, expected {expected_peak}"
        print("   ✓ PASSED")
    except Exception as e:
        print(f"   ✗ FAILED: {e}")
        all_passed = False
    
    # Test 3: Input validation
    print("\nTest 3: Input validation")
    try:
        # Should raise ValueError for negative FWHM
        try:
            generate_gaussian(fwhm=-5.0, amplitude=100.0, center=50.0, array_length=100)
            print("   ✗ FAILED: Should have raised ValueError for negative FWHM")
            all_passed = False
        except ValueError:
            print("   ✓ Correctly rejected negative FWHM")
        
        # Should raise ValueError for out-of-bounds center
        try:
            generate_gaussian(fwhm=10.0, amplitude=100.0, center=150.0, array_length=100)
            print("   ✗ FAILED: Should have raised ValueError for out-of-bounds center")
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
    """
    When run as a script, execute tests and demonstration.
    """
    import sys
    
    # Run tests
    tests_passed = run_tests()
    
    # Run demonstration if tests passed
    if tests_passed:
        print("\n")
        run_demo()
    else:
        print("\nSkipping demo due to test failures.")
        sys.exit(1)
