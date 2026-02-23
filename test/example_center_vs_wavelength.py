#!/usr/bin/env python
"""
Clear Example: center vs wavelength Parameters

This script demonstrates the difference between using 'center' and 'wavelength'
parameters, and shows how to properly use the wavelength injection feature.
"""

import numpy as np
import sys
import os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from gaussian_generator import add_gaussian_to_array

# User's wavelength array (151 elements, 0.2-1.995 microns)
wavelengths = np.array([0.2, 0.2285714286, 0.2612244898, 0.2985422741, 
                        0.3411911703, 0.3899327661, 0.445637447, 
                        0.5092999394, 0.515, 0.5186785714, 0.5223834184,
                        0.5261147285, 0.5298726908, 0.5336574958, 
                        0.537469335, 0.5413084017, 0.5451748903,
                        0.5490689967, 0.5529909181, 0.5569408532,
                        0.5609190021, 0.5649255664, 0.5689607491,
                        0.5730247544, 0.5771177884, 0.5812400583,
                        0.585391773, 0.5895731428, 0.5937843795,
                        0.5980256965, 0.6022973086, 0.6065994323,
                        0.6109322854, 0.6152960874, 0.6196910595,
                        0.6241174242, 0.6285754058, 0.6330652301,
                        0.6375871246, 0.6421413183, 0.646728042,
                        0.6513475281, 0.6560000104, 0.6606857248,
                        0.6654049085, 0.6701578007, 0.6749446422,
                        0.6797656753, 0.6846211444, 0.6895112954,
                        0.6944363761, 0.699396636, 0.7043923262,
                        0.7094237, 0.7144910121, 0.7195945193,
                        0.7247344802, 0.7299111551, 0.7351248062,
                        0.7403756976, 0.7456640955, 0.7509902676,
                        0.7563544838, 0.7617570158, 0.7671981374,
                        0.7726781241, 0.7781972535, 0.7837558053,
                        0.7893540611, 0.7949923044, 0.8006708208,
                        0.8063898981, 0.812149826, 0.8179508961,
                        0.8237934025, 0.8296776411, 0.83560391,
                        0.8415725094, 0.8475837416, 0.8536379112,
                        0.8597353248, 0.8658762914, 0.8720611221,
                        0.8782901301, 0.884563631, 0.8908819427,
                        0.8972453851, 0.9036542807, 0.9101089541,
                        0.9166097324, 0.9231569448, 0.9297509229,
                        0.936392001, 0.9430805153, 0.9498168046,
                        0.9566012104, 0.9634340762, 0.9703157482,
                        0.9772465749, 0.9842269076, 0.9912570998,
                        0.9983375077, 1.01, 1.024428571, 1.039063265,
                        1.053907026, 1.068962841, 1.084233739,
                        1.099722792, 1.115433118, 1.131367876,
                        1.147530275, 1.163923564, 1.180551044,
                        1.197416059, 1.214522002, 1.231872317,
                        1.249470493, 1.267320071, 1.285424644,
                        1.303787853, 1.322413394, 1.341305013,
                        1.360466514, 1.37990175, 1.399614632,
                        1.419609126, 1.439889257, 1.460459103,
                        1.481322805, 1.502484559, 1.523948624,
                        1.545719319, 1.567801024, 1.590198181,
                        1.612915298, 1.635956945, 1.659327758,
                        1.683032441, 1.707075761, 1.731462558,
                        1.756197737, 1.781286276, 1.806733223,
                        1.832543698, 1.858722894, 1.885276078,
                        1.912208593, 1.939525859, 1.967233371,
                        1.995336705])

flux = np.ones(len(wavelengths)) * 0.08  # Baseline flux

print("=" * 70)
print("EXAMPLE: center vs wavelength Parameters")
print("=" * 70)
print()

# Example 1: Using 'center' parameter (WRONG for wavelength 0.532!)
print("Example 1: Using center=0.532 (TREATS 0.532 AS ARRAY INDEX)")
print("-" * 70)
try:
    result1 = add_gaussian_to_array(flux.copy(), fwhm=3.0, amplitude=0.05,
                                    center=0.532)
    peak_idx1 = np.argmax(result1)
    print(f"✓ Injection successful")
    print(f"  Peak at index: {peak_idx1}")
    print(f"  Peak wavelength: {wavelengths[peak_idx1]:.6f} microns")
    print(f"  Position: {peak_idx1}/{len(wavelengths)-1} ({100*peak_idx1/(len(wavelengths)-1):.1f}% through array)")
    print()
    if peak_idx1 < 10:
        print("  ⚠️  WARNING: Peak is near the LEFT EDGE!")
        print("  This is probably NOT what you wanted!")
except Exception as e:
    print(f"✗ Error: {e}")
print()

# Example 2: Using 'wavelength' parameter (CORRECT!)
print("Example 2: Using wavelength=0.532 (INJECTS AT 0.532 MICRONS)")
print("-" * 70)
try:
    result2 = add_gaussian_to_array(flux.copy(), fwhm=3.0, amplitude=0.05,
                                    wavelength=0.532,
                                    wavelengths=wavelengths)
    peak_idx2 = np.argmax(result2)
    print(f"✓ Injection successful")
    print(f"  Peak at index: {peak_idx2}")
    print(f"  Peak wavelength: {wavelengths[peak_idx2]:.6f} microns")
    print(f"  Target wavelength: 0.532000 microns")
    print(f"  Wavelength error: {abs(wavelengths[peak_idx2] - 0.532)*1000:.2f} nm")
    print(f"  Position: {peak_idx2}/{len(wavelengths)-1} ({100*peak_idx2/(len(wavelengths)-1):.1f}% through array)")
    print()
    if 20 < peak_idx2 < len(wavelengths) - 20:
        print("  ✓ Peak is in the MIDDLE of the spectrum - correct!")
except Exception as e:
    print(f"✗ Error: {e}")
print()

# Example 3: Show the wavelength injection feature works with interpolation
print("Example 3: Wavelength NOT in array (0.550 microns)")
print("-" * 70)
try:
    # 0.550 is NOT exactly in the wavelengths array
    print(f"Is 0.550 in wavelengths array? {0.550 in wavelengths}")
    print()
    
    result3 = add_gaussian_to_array(flux.copy(), fwhm=3.0, amplitude=0.05,
                                    wavelength=0.550,
                                    wavelengths=wavelengths)
    peak_idx3 = np.argmax(result3)
    print(f"✓ Injection successful (even though 0.550 not in array!)")
    print(f"  Peak at index: {peak_idx3}")
    print(f"  Peak wavelength: {wavelengths[peak_idx3]:.6f} microns")
    print(f"  Target wavelength: 0.550000 microns")
    print(f"  Wavelength error: {abs(wavelengths[peak_idx3] - 0.550)*1000:.2f} nm")
    print()
    print("  ✓ The feature DOES support wavelengths between array elements!")
except Exception as e:
    print(f"✗ Error: {e}")
print()

# Summary
print("=" * 70)
print("SUMMARY")
print("=" * 70)
print()
print("When you want to inject at a specific WAVELENGTH:")
print()
print("  ✓ USE: wavelength=0.532, wavelengths=wavelengths")
print("  ✗ DON'T USE: center=0.532")
print()
print("The 'center' parameter is for ARRAY INDICES (0, 1, 2, ...)")
print("The 'wavelength' parameter is for PHYSICAL WAVELENGTHS (0.2, 0.532, 1.5 microns, ...)")
print()
print("Both features exist and work correctly!")
print("=" * 70)
