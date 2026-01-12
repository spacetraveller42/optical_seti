#!/usr/bin/env python3
"""
Example script demonstrating how to use the Gaussian generation functions.

This script shows how to:
1. Generate synthetic Gaussian peaks with specified FWHM
2. Add Gaussians to existing spectral data
3. Visualize the results
"""

import sys
from pathlib import Path
cwd = Path(__file__).parent.resolve() 
optical_seti_dir = cwd.parent.resolve() 
sys.path.append(str(optical_seti_dir))

import numpy as np
import matplotlib.pyplot as plt
import optical_seti_functions

# Example 1: Simple Gaussian generation
print("=" * 60)
print("Example 1: Generate a simple Gaussian")
print("=" * 60)

# Generate a Gaussian with FWHM of 15 pixels, amplitude 200, centered at pixel 100
fwhm = 15.0
amplitude = 200.0
center = 100.0
array_length = 200

gaussian = optical_seti_functions.generate_gaussian(fwhm, amplitude, center, array_length)
print(f"Generated Gaussian with FWHM={fwhm}, amplitude={amplitude}, center={center}")
print(f"Peak value: {np.max(gaussian):.2f}")
print(f"Peak location: {np.argmax(gaussian)}")

# Example 2: Add Gaussian to simulated spectral data
print("\n" + "=" * 60)
print("Example 2: Add Gaussian to simulated spectral data")
print("=" * 60)

# Create simulated spectral data with noise and a baseline
spectral_data = np.random.randn(500) * 10 + 1000  # Noise with baseline of 1000

# Add a bright, narrow Gaussian (simulating a laser signal)
laser_signal = optical_seti_functions.add_gaussian_to_array(
    spectral_data, 
    fwhm=5.0,      # Very narrow - characteristic of laser
    amplitude=500, # Bright signal
    center=250     # Middle of array
)

print("Added a narrow Gaussian (FWHM=5.0) to simulated spectral data")
print(f"Original data range: [{np.min(spectral_data):.2f}, {np.max(spectral_data):.2f}]")
print(f"Data with signal range: [{np.min(laser_signal):.2f}, {np.max(laser_signal):.2f}]")

# Example 3: Multiple signals with different FWHM values
print("\n" + "=" * 60)
print("Example 3: Add multiple signals with different FWHM")
print("=" * 60)

# Start with a clean baseline
baseline = np.random.randn(1000) * 5 + 500

# Add narrow signal (laser-like)
data = optical_seti_functions.add_gaussian_to_array(baseline, fwhm=4.0, amplitude=300, center=200)
print("Added narrow signal (FWHM=4.0) at pixel 200")

# Add medium width signal (emission line)
data = optical_seti_functions.add_gaussian_to_array(data, fwhm=12.0, amplitude=200, center=500)
print("Added medium signal (FWHM=12.0) at pixel 500")

# Add wide signal (broad feature)
data = optical_seti_functions.add_gaussian_to_array(data, fwhm=30.0, amplitude=150, center=800)
print("Added wide signal (FWHM=30.0) at pixel 800")

# Visualization
print("\n" + "=" * 60)
print("Creating visualization...")
print("=" * 60)

fig, axes = plt.subplots(3, 1, figsize=(12, 10))

# Plot 1: Simple Gaussian
axes[0].plot(gaussian, 'b-', linewidth=2)
axes[0].axhline(y=amplitude/2, color='r', linestyle='--', alpha=0.5, label=f'Half Max (FWHM={fwhm})')
axes[0].set_title('Example 1: Simple Gaussian Generation', fontsize=14, fontweight='bold')
axes[0].set_xlabel('Array Index')
axes[0].set_ylabel('Amplitude')
axes[0].legend()
axes[0].grid(True, alpha=0.3)

# Plot 2: Single signal in noisy data
axes[1].plot(spectral_data, 'gray', alpha=0.4, linewidth=0.5, label='Original noisy data')
axes[1].plot(laser_signal, 'b-', linewidth=1.5, label='With narrow laser signal (FWHM=5)')
axes[1].set_title('Example 2: Narrow Gaussian Added to Spectral Data', fontsize=14, fontweight='bold')
axes[1].set_xlabel('Pixel Index')
axes[1].set_ylabel('Flux')
axes[1].legend()
axes[1].grid(True, alpha=0.3)

# Plot 3: Multiple signals
axes[2].plot(baseline, 'lightgray', alpha=0.5, linewidth=0.5, label='Baseline')
axes[2].plot(data, 'b-', linewidth=1.5, label='With 3 Gaussians (FWHM: 4, 12, 30)')
axes[2].axvline(x=200, color='r', linestyle=':', alpha=0.5)
axes[2].axvline(x=500, color='r', linestyle=':', alpha=0.5)
axes[2].axvline(x=800, color='r', linestyle=':', alpha=0.5)
axes[2].set_title('Example 3: Multiple Gaussians with Different FWHM', fontsize=14, fontweight='bold')
axes[2].set_xlabel('Pixel Index')
axes[2].set_ylabel('Flux')
axes[2].legend()
axes[2].grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig('gaussian_examples.png', dpi=150, bbox_inches='tight')
print("Visualization saved as 'gaussian_examples.png'")

print("\n" + "=" * 60)
print("All examples completed successfully!")
print("=" * 60)

# Show plot
plt.show()
