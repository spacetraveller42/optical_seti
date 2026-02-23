#!/usr/bin/env python
"""
Example usage of find_closest_wavelength_index() function.

This script demonstrates how to use the new wavelength index finder function.
"""

import numpy as np
import sys
import os

# Add parent directory to path
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from gaussian_generator import find_closest_wavelength_index, add_gaussian_to_array

print("="*70)
print("EXAMPLE: Using find_closest_wavelength_index()")
print("="*70)

# User's actual 151-point wavelength array (0.2-1.995 microns)
wavelengths = np.array([0.2, 0.2285714286, 0.2612244898, 0.2985422741, 
                       0.3411911703, 0.3899327661, 0.445637447, 0.5092999394, 
                       0.515, 0.5186785714, 0.5223834184, 0.5261147285, 
                       0.5298726908, 0.5336574958, 0.537469335, 0.5413084017, 
                       0.5451748903, 0.5490689967, 0.5529909181, 0.5569408532, 
                       0.5609190021, 0.5649255664, 0.5689607491, 0.5730247544, 
                       0.5771177884, 0.5812400583, 0.585391773, 0.5895731428, 
                       0.5937843795, 0.5980256965, 0.6022973086, 0.6065994323, 
                       0.6109322854, 0.6152960874, 0.6196910595, 0.6241174242, 
                       0.6285754058, 0.6330652301, 0.6375871246, 0.6421413183, 
                       0.646728042, 0.6513475281, 0.6560000104, 0.6606857248, 
                       0.6654049085, 0.6701578007, 0.6749446422, 0.6797656753, 
                       0.6846211444, 0.6895112954, 0.6944363761, 0.699396636, 
                       0.7043923262, 0.7094237, 0.7144910121, 0.7195945193, 
                       0.7247344802, 0.7299111551, 0.7351248062, 0.7403756976, 
                       0.7456640955, 0.7509902676, 0.7563544838, 0.7617570158, 
                       0.7671981374, 0.7726781241, 0.7781972535, 0.7837558053, 
                       0.7893540611, 0.7949923044, 0.8006708208, 0.8063898981, 
                       0.812149826, 0.8179508961, 0.8237934025, 0.8296776411, 
                       0.83560391, 0.8415725094, 0.8475837416, 0.8536379112, 
                       0.8597353248, 0.8658762914, 0.8720611221, 0.8782901301, 
                       0.884563631, 0.8908819427, 0.8972453851, 0.9036542807, 
                       0.9101089541, 0.9166097324, 0.9231569448, 0.9297509229, 
                       0.936392001, 0.9430805153, 0.9498168046, 0.9566012104, 
                       0.9634340762, 0.9703157482, 0.9772465749, 0.9842269076, 
                       0.9912570998, 0.9983375077, 1.01, 1.024428571, 
                       1.039063265, 1.053907026, 1.068962841, 1.084233739, 
                       1.099722792, 1.115433118, 1.131367876, 1.147530275, 
                       1.163923564, 1.180551044, 1.197416059, 1.214522002, 
                       1.231872317, 1.249470493, 1.267320071, 1.285424644, 
                       1.303787853, 1.322413394, 1.341305013, 1.360466514, 
                       1.37990175, 1.399614632, 1.419609126, 1.439889257, 
                       1.460459103, 1.481322805, 1.502484559, 1.523948624, 
                       1.545719319, 1.567801024, 1.590198181, 1.612915298, 
                       1.635956945, 1.659327758, 1.683032441, 1.707075761, 
                       1.731462558, 1.756197737, 1.781286276, 1.806733223, 
                       1.832543698, 1.858722894, 1.885276078, 1.912208593, 
                       1.939525859, 1.967233371, 1.995336705])

# Create some example flux data
flux = np.ones(len(wavelengths)) * 0.08  # Baseline flux

print(f"\nWavelength array: {len(wavelengths)} points")
print(f"Range: {wavelengths[0]:.3f} to {wavelengths[-1]:.3f} microns")
print(f"Flux array: {len(flux)} points, mean = {flux.mean():.3f}")

# ============================================================================
# Example 1: Find closest wavelength to green laser (0.532 microns)
# ============================================================================
print("\n" + "-"*70)
print("EXAMPLE 1: Find closest wavelength to 0.532 microns (green laser)")
print("-"*70)

target = 0.532
idx, wl, diff = find_closest_wavelength_index(target, wavelengths)

print(f"\nTarget wavelength: {target} microns")
print(f"Closest index: {idx}")
print(f"Closest wavelength: {wl:.6f} microns")
print(f"Difference: {diff:.6f} microns ({diff*1000:.2f} nm)")

# Extract flux at that wavelength
flux_at_wavelength = flux[idx]
print(f"Flux at that wavelength: {flux_at_wavelength:.6f}")

# ============================================================================
# Example 2: Find multiple wavelengths (laser calibration lines)
# ============================================================================
print("\n" + "-"*70)
print("EXAMPLE 2: Find multiple calibration wavelengths")
print("-"*70)

calibration_wavelengths = [0.532, 0.633, 0.780, 1.064]  # Common laser lines

print("\nCalibration wavelengths:")
for target_wl in calibration_wavelengths:
    idx, wl, diff = find_closest_wavelength_index(target_wl, wavelengths)
    print(f"  {target_wl:.3f} μm → index {idx:3d}, actual {wl:.6f} μm, diff {diff*1000:5.2f} nm")

# ============================================================================
# Example 3: Use with Gaussian injection
# ============================================================================
print("\n" + "-"*70)
print("EXAMPLE 3: Combine with Gaussian injection")
print("-"*70)

target = 0.532
print(f"\nTarget: {target} microns")

# Find the closest wavelength index
idx, wl, diff = find_closest_wavelength_index(target, wavelengths)
print(f"Closest index: {idx}, wavelength: {wl:.6f} microns")

# Inject Gaussian at that wavelength
result = add_gaussian_to_array(flux, fwhm=3.0, amplitude=0.05, 
                               wavelength=target, wavelengths=wavelengths)

# Verify the peak is at the expected location
peak_idx = np.argmax(result)
print(f"\nAfter injection:")
print(f"  Peak at index: {peak_idx}")
print(f"  Peak wavelength: {wavelengths[peak_idx]:.6f} microns")
print(f"  Peak flux: {result[peak_idx]:.6f}")
print(f"  Baseline flux: {flux.mean():.6f}")
print(f"  Increase: {result[peak_idx] - flux.mean():.6f}")

# ============================================================================
# Example 4: Extract spectral features
# ============================================================================
print("\n" + "-"*70)
print("EXAMPLE 4: Extract spectral features at specific wavelengths")
print("-"*70)

# Add some features to the spectrum
result = flux.copy()
result = add_gaussian_to_array(result, fwhm=3.0, amplitude=0.02, 
                              wavelength=0.589, wavelengths=wavelengths)  # Sodium D-line
result = add_gaussian_to_array(result, fwhm=2.0, amplitude=0.03, 
                              wavelength=0.656, wavelengths=wavelengths)  # H-alpha

# Now extract features at specific wavelengths
features_to_extract = {
    'Sodium D': 0.589,
    'H-alpha': 0.656,
    'Random point': 1.2
}

print("\nExtracting spectral features:")
for name, target_wl in features_to_extract.items():
    idx, wl, diff = find_closest_wavelength_index(target_wl, wavelengths)
    flux_value = result[idx]
    print(f"  {name:15s} ({target_wl:.3f} μm) → index {idx:3d}, "
          f"flux {flux_value:.6f}, diff {diff*1000:5.2f} nm")

# ============================================================================
# Example 5: Verify wavelength range before injection
# ============================================================================
print("\n" + "-"*70)
print("EXAMPLE 5: Check if wavelength is in range before injection")
print("-"*70)

test_wavelengths = [0.15, 0.532, 2.5]  # Below, in, above range

for target in test_wavelengths:
    print(f"\nTarget: {target} microns")
    
    # Check if in range
    if target < wavelengths.min() or target > wavelengths.max():
        print(f"  ⚠️  OUT OF RANGE! (valid: {wavelengths.min():.3f}-{wavelengths.max():.3f} μm)")
    else:
        idx, wl, diff = find_closest_wavelength_index(target, wavelengths)
        print(f"  ✓ In range: closest index {idx}, wavelength {wl:.6f} μm")

print("\n" + "="*70)
print("✓ All examples completed successfully!")
print("="*70)
