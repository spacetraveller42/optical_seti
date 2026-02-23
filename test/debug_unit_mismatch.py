#!/usr/bin/env python
"""
Test what happens with unit mismatch - using 532 nm when array is in microns.
This might be causing the "edge of spectrum" issue!
"""

import sys
import os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import numpy as np
from gaussian_generator import add_gaussian_to_array

# User's exact 151-point wavelength array (0.2 to 1.995 microns)
wavelengths = np.array([0.2, 0.2285714286, 0.2612244898, 0.2985422741, 0.3411911703, 
                       0.3899327661, 0.445637447, 0.5092999394, 0.515, 0.5186785714, 
                       0.5223834184, 0.5261147285, 0.5298726908, 0.5336574958, 0.537469335, 
                       0.5413084017, 0.5451748903, 0.5490689967, 0.5529909181, 0.5569408532, 
                       0.5609190021, 0.5649255664, 0.5689607491, 0.5730247544, 0.5771177884, 
                       0.5812400583, 0.585391773, 0.5895731428, 0.5937843795, 0.5980256965, 
                       0.6022973086, 0.6065994323, 0.6109322854, 0.6152960874, 0.6196910595, 
                       0.6241174242, 0.6285754058, 0.6330652301, 0.6375871246, 0.6421413183, 
                       0.646728042, 0.6513475281, 0.6560000104, 0.6606857248, 0.6654049085, 
                       0.6701578007, 0.6749446422, 0.6797656753, 0.6846211444, 0.6895112954, 
                       0.6944363761, 0.699396636, 0.7043923262, 0.7094237, 0.7144910121, 
                       0.7195945193, 0.7247344802, 0.7299111551, 0.7351248062, 0.7403756976, 
                       0.7456640955, 0.7509902676, 0.7563544838, 0.7617570158, 0.7671981374, 
                       0.7726781241, 0.7781972535, 0.7837558053, 0.7893540611, 0.7949923044, 
                       0.8006708208, 0.8063898981, 0.812149826, 0.8179508961, 0.8237934025, 
                       0.8296776411, 0.83560391, 0.8415725094, 0.8475837416, 0.8536379112, 
                       0.8597353248, 0.8658762914, 0.8720611221, 0.8782901301, 0.884563631, 
                       0.8908819427, 0.8972453851, 0.9036542807, 0.9101089541, 0.9166097324, 
                       0.9231569448, 0.9297509229, 0.936392001, 0.9430805153, 0.9498168046, 
                       0.9566012104, 0.9634340762, 0.9703157482, 0.9772465749, 0.9842269076, 
                       0.9912570998, 0.9983375077, 1.01, 1.024428571, 1.039063265, 
                       1.053907026, 1.068962841, 1.084233739, 1.099722792, 1.115433118, 
                       1.131367876, 1.147530275, 1.163923564, 1.180551044, 1.197416059, 
                       1.214522002, 1.231872317, 1.249470493, 1.267320071, 1.285424644, 
                       1.303787853, 1.322413394, 1.341305013, 1.360466514, 1.37990175, 
                       1.399614632, 1.419609126, 1.439889257, 1.460459103, 1.481322805, 
                       1.502484559, 1.523948624, 1.545719319, 1.567801024, 1.590198181, 
                       1.612915298, 1.635956945, 1.659327758, 1.683032441, 1.707075761, 
                       1.731462558, 1.756197737, 1.781286276, 1.806733223, 1.832543698, 
                       1.858722894, 1.885276078, 1.912208593, 1.939525859, 1.967233371, 
                       1.995336705])

print("=" * 70)
print("TESTING UNIT MISMATCH: wavelength=532 (nm) vs array in microns")
print("=" * 70)
print()

print(f"Wavelength array range: {wavelengths.min():.3f} to {wavelengths.max():.3f} microns")
print(f"Target: 532 (likely meant as 532 nm = 0.532 microns)")
print()

flux = np.ones(len(wavelengths)) * 0.08

# Test 1: Using 532 (what user might be doing - WRONG!)
print("Test 1: wavelength=532 (UNIT MISMATCH!)")
print("-" * 70)
target_wrong = 532  # nanometers interpreted as microns

print(f"Target wavelength: {target_wrong}")
print(f"Is in range [{wavelengths.min():.3f}, {wavelengths.max():.3f}]: ", end="")

if wavelengths.min() <= target_wrong <= wavelengths.max():
    print("YES - but shouldn't be!")
    print("⚠️  WARNING: 532 > 1.995, should fail but let's see what np.interp does...")
else:
    print("NO")
    print(f"✗ Out of range! {target_wrong} is way outside [{wavelengths.min()}, {wavelengths.max()}]")
    print()
    
    # np.interp extrapolates outside range!
    indices = np.arange(len(wavelengths))
    extrapolated_idx = np.interp(target_wrong, wavelengths, indices)
    print(f"np.interp still returns: {extrapolated_idx}")
    print(f"This is BEYOND the array (max index = {len(wavelengths)-1})")
    print()
    
    # What happens when we try to inject?
    print("Attempting injection with wavelength=532...")
    try:
        result = add_gaussian_to_array(
            flux,
            fwhm=3.0,
            amplitude=0.05,
            wavelength=target_wrong,
            wavelengths=wavelengths
        )
        peak_idx = np.argmax(result)
        peak_wave = wavelengths[peak_idx]
        print(f"✓ Injection succeeded (shouldn't have!)")
        print(f"  Peak at index: {peak_idx}")
        print(f"  Peak wavelength: {peak_wave:.6f} microns")
        
        if peak_idx > len(wavelengths) - 10:
            print(f"  ⚠️  FOUND IT! Peak is at RIGHT EDGE (index {peak_idx}/{len(wavelengths)-1})")
            print(f"  This is why user sees laser at edge of spectrum!")
        
    except ValueError as e:
        print(f"✗ Correctly rejected: {e}")

print()
print("=" * 70)
print("Test 2: wavelength=0.532 (CORRECT - microns)")
print("-" * 70)

target_correct = 0.532
print(f"Target wavelength: {target_correct} microns")
print(f"Is in range: {wavelengths.min() <= target_correct <= wavelengths.max()}")

result = add_gaussian_to_array(
    flux,
    fwhm=3.0,
    amplitude=0.05,
    wavelength=target_correct,
    wavelengths=wavelengths
)

peak_idx = np.argmax(result)
peak_wave = wavelengths[peak_idx]
print(f"✓ Injection successful!")
print(f"  Peak at index: {peak_idx}")
print(f"  Peak wavelength: {peak_wave:.6f} microns")
print(f"  Peak is in MIDDLE of spectrum (index {peak_idx}/{len(wavelengths)-1})")

print()
print("=" * 70)
print("CONCLUSION")
print("=" * 70)
print()
print("If user sees laser at edge of spectrum, they are likely using:")
print("  wavelength=532  (nanometers, but array is in microns)")
print()
print("Instead of:")
print("  wavelength=0.532  (microns, matching array units)")
print()
print("This is a UNIT MISMATCH error, not a code bug!")
print()
