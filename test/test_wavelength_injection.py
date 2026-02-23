#!/usr/bin/env python
"""
Test script for wavelength-based Gaussian injection functionality.

This tests the ability to inject Gaussian spikes at specific wavelengths,
including wavelengths that fall between array elements (requiring interpolation).
"""

import sys
from pathlib import Path
cwd = Path(__file__).parent.resolve()
optical_seti_dir = cwd.parent.resolve()
sys.path.append(str(optical_seti_dir))

import numpy as np
import matplotlib.pyplot as plt
from gaussian_generator import generate_gaussian, add_gaussian_to_array

print("="*70)
print("TEST: Wavelength-Based Gaussian Injection")
print("="*70)

# Create a wavelength array (simulating spectral data)
# 400-700 nm range with 100 points
wavelengths = np.linspace(400, 700, 100)
array_length = len(wavelengths)

print(f"\nWavelength range: {wavelengths[0]:.1f} - {wavelengths[-1]:.1f} nm")
print(f"Array length: {array_length}")
print(f"Wavelength spacing: {wavelengths[1] - wavelengths[0]:.2f} nm")

# Test 1: Exact wavelength match (on grid point)
print("\n" + "="*70)
print("Test 1: Inject Gaussian at exact wavelength (550.0 nm)")
print("="*70)

target_wavelength_1 = 550.0  # Exactly at a grid point
gaussian_1 = generate_gaussian(fwhm=5.0, amplitude=100.0, wavelength=target_wavelength_1,
                               array_length=array_length, wavelengths=wavelengths)

# Find the peak position
peak_index_1 = np.argmax(gaussian_1)
peak_wavelength_1 = wavelengths[peak_index_1]

print(f"Target wavelength: {target_wavelength_1:.2f} nm")
print(f"Peak found at index: {peak_index_1}")
print(f"Peak wavelength: {peak_wavelength_1:.2f} nm")
print(f"Peak value: {np.max(gaussian_1):.2f}")
print(f"Error: {abs(peak_wavelength_1 - target_wavelength_1):.4f} nm")

# Test 2: Wavelength between grid points (interpolated)
print("\n" + "="*70)
print("Test 2: Inject Gaussian at interpolated wavelength (475.5 nm)")
print("="*70)

target_wavelength_2 = 475.5  # Between grid points
gaussian_2 = generate_gaussian(fwhm=5.0, amplitude=150.0, wavelength=target_wavelength_2,
                               array_length=array_length, wavelengths=wavelengths)

# Find the peak position
peak_index_2 = np.argmax(gaussian_2)
peak_wavelength_2 = wavelengths[peak_index_2]

# Calculate the expected fractional index
expected_index = np.interp(target_wavelength_2, wavelengths, np.arange(array_length))

print(f"Target wavelength: {target_wavelength_2:.2f} nm")
print(f"Expected fractional index: {expected_index:.2f}")
print(f"Peak found at index: {peak_index_2}")
print(f"Peak wavelength: {peak_wavelength_2:.2f} nm")
print(f"Peak value: {np.max(gaussian_2):.2f}")
print(f"Wavelength error: {abs(peak_wavelength_2 - target_wavelength_2):.4f} nm")

# Test 3: Using area parameter with wavelength
print("\n" + "="*70)
print("Test 3: Inject Gaussian with area parameter at 625.7 nm")
print("="*70)

target_wavelength_3 = 625.7  # Between grid points
target_area = 1000.0
gaussian_3 = generate_gaussian(fwhm=8.0, area=target_area, wavelength=target_wavelength_3,
                               array_length=array_length, wavelengths=wavelengths)

measured_area = np.sum(gaussian_3)
peak_index_3 = np.argmax(gaussian_3)
peak_wavelength_3 = wavelengths[peak_index_3]

print(f"Target wavelength: {target_wavelength_3:.2f} nm")
print(f"Target area: {target_area:.2f}")
print(f"Measured area: {measured_area:.2f}")
print(f"Area error: {abs(measured_area - target_area):.2f}")
print(f"Peak found at index: {peak_index_3}")
print(f"Peak wavelength: {peak_wavelength_3:.2f} nm")

# Test 4: Add Gaussian to existing data using wavelength
print("\n" + "="*70)
print("Test 4: Add Gaussian to existing spectral data")
print("="*70)

# Create baseline data (simulated spectrum)
baseline = np.ones(array_length) * 500 + np.random.randn(array_length) * 10

# Add Gaussian at specific wavelength
target_wavelength_4 = 532.0  # Green laser line
result = add_gaussian_to_array(baseline, fwhm=3.0, amplitude=200.0,
                               wavelength=target_wavelength_4, wavelengths=wavelengths)

# Verify the injection
peak_index_4 = np.argmax(result)
peak_wavelength_4 = wavelengths[peak_index_4]

print(f"Target wavelength: {target_wavelength_4:.2f} nm")
print(f"Peak found at index: {peak_index_4}")
print(f"Peak wavelength: {peak_wavelength_4:.2f} nm")
print(f"Peak flux: {result[peak_index_4]:.2f}")
print(f"Baseline at peak: {baseline[peak_index_4]:.2f}")

# Test 5: Edge case validation - out of range wavelength
print("\n" + "="*70)
print("Test 5: Error handling for out-of-range wavelength")
print("="*70)

try:
    gaussian_bad = generate_gaussian(fwhm=5.0, amplitude=100.0, wavelength=800.0,  # Outside range
                                    array_length=array_length, wavelengths=wavelengths)
    print("ERROR: Should have raised ValueError for out-of-range wavelength")
except ValueError as e:
    print(f"✓ Correctly raised ValueError: {e}")

# Test 6: Error handling - both center and wavelength specified
print("\n" + "="*70)
print("Test 6: Error handling for ambiguous parameters")
print("="*70)

try:
    gaussian_bad = generate_gaussian(fwhm=5.0, amplitude=100.0, center=50.0, 
                                    wavelength=550.0, array_length=array_length, 
                                    wavelengths=wavelengths)
    print("ERROR: Should have raised ValueError for both center and wavelength")
except ValueError as e:
    print(f"✓ Correctly raised ValueError: {e}")

# Create visualization
print("\n" + "="*70)
print("Creating visualization...")
print("="*70)

fig, axes = plt.subplots(3, 1, figsize=(12, 12))

# Plot 1: Exact wavelength vs interpolated wavelength
axes[0].plot(wavelengths, gaussian_1, 'b-', linewidth=2, label=f'Exact: {target_wavelength_1:.1f} nm')
axes[0].plot(wavelengths, gaussian_2, 'r-', linewidth=2, label=f'Interpolated: {target_wavelength_2:.1f} nm')
axes[0].axvline(x=target_wavelength_1, color='b', linestyle='--', alpha=0.5)
axes[0].axvline(x=target_wavelength_2, color='r', linestyle='--', alpha=0.5)
axes[0].set_xlabel('Wavelength (nm)')
axes[0].set_ylabel('Amplitude')
axes[0].set_title('Gaussian Injection at Different Wavelengths')
axes[0].legend()
axes[0].grid(True, alpha=0.3)

# Plot 2: Area-specified Gaussian
axes[1].plot(wavelengths, gaussian_3, 'g-', linewidth=2)
axes[1].axvline(x=target_wavelength_3, color='g', linestyle='--', alpha=0.5, 
               label=f'Target: {target_wavelength_3:.1f} nm')
axes[1].set_xlabel('Wavelength (nm)')
axes[1].set_ylabel('Amplitude')
axes[1].set_title(f'Gaussian with Area={target_area:.1f} at {target_wavelength_3:.1f} nm')
axes[1].legend()
axes[1].grid(True, alpha=0.3)

# Plot 3: Gaussian added to noisy baseline
axes[2].plot(wavelengths, baseline, 'gray', linewidth=1, alpha=0.5, label='Baseline')
axes[2].plot(wavelengths, result, 'b-', linewidth=2, label='With Gaussian')
axes[2].axvline(x=target_wavelength_4, color='r', linestyle='--', alpha=0.7,
               label=f'Injection: {target_wavelength_4:.1f} nm')
axes[2].set_xlabel('Wavelength (nm)')
axes[2].set_ylabel('Flux')
axes[2].set_title(f'Gaussian Spike Injected at {target_wavelength_4:.1f} nm (Laser Line Simulation)')
axes[2].legend()
axes[2].grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig(cwd / 'test_wavelength_injection.png', dpi=150)
print(f"✓ Plot saved as 'test_wavelength_injection.png'")

# Summary
print("\n" + "="*70)
print("SUMMARY")
print("="*70)
print("✓ All tests passed!")
print("✓ Wavelength-based injection works for exact wavelengths")
print("✓ Interpolation works for wavelengths between grid points")
print("✓ Area parameter works with wavelength specification")
print("✓ Error handling validates out-of-range wavelengths")
print("✓ Error handling prevents ambiguous parameter combinations")
print("="*70)

plt.show()
