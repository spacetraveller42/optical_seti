#!/usr/bin/env python
"""
Example: Using wavelength-based injection with HARPS spectral data

This example shows how to inject a Gaussian spike into actual HARPS spectral data
at a specific wavelength (e.g., simulating a laser reference line).
"""

import sys
from pathlib import Path
cwd = Path(__file__).parent.resolve()
optical_seti_dir = cwd.parent.resolve()
sys.path.append(str(optical_seti_dir))

import numpy as np
from gaussian_generator import add_gaussian_to_array

# Simulate HARPS data format
# In real usage, you would load this with: wave, flux = read_harps_file(filename)
wave = np.linspace(400, 700, 1000)  # Wavelength array in nm
flux = np.ones(1000) * 1000 + np.random.randn(1000) * 50  # Flux array

print("HARPS Data Simulation")
print("=" * 70)
print(f"Wavelength range: {wave[0]:.1f} - {wave[-1]:.1f} nm")
print(f"Number of pixels: {len(wave)}")
print(f"Wavelength spacing: ~{(wave[-1] - wave[0])/len(wave):.3f} nm/pixel")

# Inject a laser calibration line at 632.8 nm (He-Ne laser)
target_wavelength = 632.8  # nm
print(f"\nInjecting calibration line at {target_wavelength} nm...")

# Add the Gaussian spike
flux_with_laser = add_gaussian_to_array(
    flux,
    fwhm=3.0,              # Line width in pixels
    amplitude=2000.0,       # Peak amplitude above continuum
    wavelength=target_wavelength,  # Target wavelength in nm
    wavelengths=wave       # Wavelength array from HARPS data
)

# Verify the injection
peak_index = np.argmax(flux_with_laser)
peak_wavelength = wave[peak_index]
peak_flux = flux_with_laser[peak_index]

print("\nResults:")
print(f"✓ Target wavelength: {target_wavelength} nm")
print(f"✓ Peak found at index: {peak_index}")
print(f"✓ Peak wavelength: {peak_wavelength:.2f} nm")
print(f"✓ Peak flux: {peak_flux:.1f}")
print(f"✓ Wavelength error: {abs(peak_wavelength - target_wavelength):.3f} nm")

print("\n" + "=" * 70)
print("USAGE WITH REAL HARPS DATA:")
print("=" * 70)
print("""
# Load HARPS spectral data
wave, flux = read_harps_file('your_harps_file.fits')

# Inject Gaussian at specific wavelength
flux_with_spike = add_gaussian_to_array(
    flux,
    fwhm=5.0,
    amplitude=1000.0,
    wavelength=589.0,  # Sodium D-line
    wavelengths=wave
)

# The wavelength can be:
# - Exactly on a grid point
# - Between grid points (will interpolate)
# - Any value within the wavelength range
""")
print("=" * 70)
