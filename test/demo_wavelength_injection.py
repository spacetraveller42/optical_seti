#!/usr/bin/env python
"""
Demonstration: Inject Gaussian at specific wavelength (including between grid points)

This example shows how to inject a Gaussian spike at any wavelength,
including wavelengths that fall between array elements.
"""

import sys
from pathlib import Path
cwd = Path(__file__).parent.resolve()
optical_seti_dir = cwd.parent.resolve()
sys.path.append(str(optical_seti_dir))

import numpy as np
import matplotlib.pyplot as plt
from gaussian_generator import generate_gaussian, add_gaussian_to_array

# Create a wavelength array (simulating HARPS spectral data)
# 500-600 nm range with 50 points
wavelengths = np.linspace(500, 600, 50)
print(f"Wavelength range: {wavelengths[0]:.1f} - {wavelengths[-1]:.1f} nm")
print(f"Wavelength spacing: {wavelengths[1] - wavelengths[0]:.2f} nm/pixel")

# Create baseline spectral data (simulated continuum with noise)
baseline_flux = 1000 + np.random.randn(len(wavelengths)) * 20

# Example 1: Inject laser line at 532 nm (green laser)
# This wavelength falls between grid points
laser_wavelength = 532.0  # nm
print(f"\nInjecting laser line at {laser_wavelength} nm...")

# Add Gaussian spike at the laser wavelength
result = add_gaussian_to_array(
    baseline_flux,
    fwhm=2.0,           # Narrow line (2 pixels FWHM)
    amplitude=500.0,     # Peak amplitude
    wavelength=laser_wavelength,
    wavelengths=wavelengths
)

# Example 2: Inject multiple laser lines at different wavelengths
result_multi = np.copy(baseline_flux)

# Green laser (532 nm)
result_multi = add_gaussian_to_array(
    result_multi,
    fwhm=2.0,
    amplitude=500.0,
    wavelength=532.0,
    wavelengths=wavelengths
)

# Red laser (632.8 nm - He-Ne)
result_multi = add_gaussian_to_array(
    result_multi,
    fwhm=2.5,
    amplitude=400.0,
    wavelength=565.5,  # Between grid points
    wavelengths=wavelengths
)

# Create visualization
fig, axes = plt.subplots(2, 1, figsize=(14, 10))

# Plot 1: Single laser injection
axes[0].plot(wavelengths, baseline_flux, 'gray', linewidth=1.5, alpha=0.6, label='Baseline spectrum')
axes[0].plot(wavelengths, result, 'b-', linewidth=2, label='With injected laser line')
axes[0].axvline(x=laser_wavelength, color='red', linestyle='--', linewidth=2, alpha=0.7,
               label=f'Target: {laser_wavelength} nm')
axes[0].set_xlabel('Wavelength (nm)', fontsize=12)
axes[0].set_ylabel('Flux (counts)', fontsize=12)
axes[0].set_title('Gaussian Spike Injected at Specific Wavelength\n(Wavelength can be between grid points)', 
                 fontsize=14, fontweight='bold')
axes[0].legend(fontsize=11)
axes[0].grid(True, alpha=0.3)

# Zoom in on the laser line
zoom_range = 15  # nm
mask = (wavelengths > laser_wavelength - zoom_range) & (wavelengths < laser_wavelength + zoom_range)
axes[0].set_xlim(laser_wavelength - zoom_range, laser_wavelength + zoom_range)

# Plot 2: Multiple laser injections
axes[1].plot(wavelengths, baseline_flux, 'gray', linewidth=1.5, alpha=0.6, label='Baseline spectrum')
axes[1].plot(wavelengths, result_multi, 'b-', linewidth=2, label='With multiple laser lines')
axes[1].axvline(x=532.0, color='green', linestyle='--', linewidth=1.5, alpha=0.7, label='532.0 nm')
axes[1].axvline(x=565.5, color='orange', linestyle='--', linewidth=1.5, alpha=0.7, label='565.5 nm (interpolated)')
axes[1].set_xlabel('Wavelength (nm)', fontsize=12)
axes[1].set_ylabel('Flux (counts)', fontsize=12)
axes[1].set_title('Multiple Gaussian Spikes at Different Wavelengths', fontsize=14, fontweight='bold')
axes[1].legend(fontsize=11)
axes[1].grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig(cwd / 'demo_wavelength_injection.png', dpi=150, bbox_inches='tight')
print(f"\n✓ Visualization saved as 'demo_wavelength_injection.png'")

print("\n" + "="*70)
print("KEY FEATURES:")
print("="*70)
print("✓ Inject Gaussian at ANY wavelength (not limited to grid points)")
print("✓ Wavelength 532.0 nm is automatically mapped to correct array index")
print("✓ Interpolation handles wavelengths between grid points")
print("✓ Perfect for simulating laser lines, emission lines, etc.")
print("✓ Works with real HARPS spectral data")
print("="*70)

plt.show()
