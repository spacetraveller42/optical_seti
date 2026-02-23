#!/usr/bin/env python
"""
COMPREHENSIVE DIAGNOSTIC FOR EDGE INJECTION ISSUE

Run this script with your actual data to diagnose why the laser
is appearing at the edge of the spectrum.

Usage:
    python diagnose_my_issue.py

This script will:
1. Check your wavelength array properties
2. Test injection at 0.532 microns
3. Identify the exact cause of the issue
4. Provide specific fix recommendations
"""

import sys
import os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import numpy as np
from gaussian_generator import add_gaussian_to_array

def diagnose_injection_issue():
    """
    Interactive diagnostic tool for wavelength injection issues.
    """
    print("=" * 80)
    print("GAUSSIAN INJECTION DIAGNOSTIC TOOL")
    print("=" * 80)
    print()
    print("This tool will help identify why your laser injection isn't working.")
    print()
    
    # Step 1: Get wavelength array
    print("STEP 1: Wavelength Array")
    print("-" * 80)
    print()
    print("Please provide your wavelength array.")
    print()
    print("Options:")
    print("  1. Use the 151-element array you provided earlier (0.2 to 1.995 microns)")
    print("  2. Enter your own wavelength array")
    print()
    
    choice = input("Enter choice (1 or 2): ").strip()
    
    if choice == "1":
        # User's provided wavelength array
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
        print("✓ Using your 151-element wavelength array")
    else:
        print("Enter wavelengths as comma-separated values (e.g., 0.4, 0.5, 0.6, 0.7):")
        wave_str = input("> ")
        try:
            wavelengths = np.array([float(x.strip()) for x in wave_str.split(",")])
            print(f"✓ Loaded {len(wavelengths)} wavelengths")
        except:
            print("✗ Error parsing wavelengths. Using default array.")
            return
    
    print()
    print("Wavelength Array Properties:")
    print(f"  Length: {len(wavelengths)}")
    print(f"  Range: {wavelengths.min():.6f} to {wavelengths.max():.6f}")
    print(f"  Is monotonic: {np.all(np.diff(wavelengths) > 0) or np.all(np.diff(wavelengths) < 0)}")
    print()
    
    # Step 2: Get target wavelength
    print("STEP 2: Target Wavelength")
    print("-" * 80)
    print()
    print("What wavelength are you trying to inject?")
    print()
    print("Common laser wavelengths:")
    print("  532 nm = 0.532 microns (green)")
    print("  632.8 nm = 0.6328 microns (red, He-Ne)")
    print()
    
    target_str = input("Enter target wavelength: ").strip()
    try:
        target = float(target_str)
    except:
        print(f"✗ Could not parse '{target_str}' as number")
        return
    
    print()
    print(f"Target wavelength: {target}")
    print()
    
    # Step 3: Check units
    print("STEP 3: Unit Check")
    print("-" * 80)
    print()
    
    # Detect likely unit mismatch
    if target > 100:
        print(f"⚠️  WARNING: Target value {target} looks like NANOMETERS")
        print(f"   But your wavelength array range is {wavelengths.min():.3f} to {wavelengths.max():.3f}")
        print()
        print("   This suggests your array is in MICRONS, not nanometers!")
        print()
        print(f"   Convert: {target} nm = {target/1000:.6f} microns")
        print()
        
        convert = input("Would you like to convert to microns? (y/n): ").strip().lower()
        if convert == 'y':
            target = target / 1000
            print(f"✓ Converted to {target} microns")
            print()
    
    # Check if in range
    wave_min, wave_max = wavelengths.min(), wavelengths.max()
    in_range = wave_min <= target <= wave_max
    
    print(f"Is target in wavelength array range? {in_range}")
    if not in_range:
        print(f"✗ TARGET IS OUT OF RANGE!")
        print(f"   Target: {target}")
        print(f"   Array range: [{wave_min:.6f}, {wave_max:.6f}]")
        print()
        print("This will cause an error or edge injection.")
        print()
        print("SOLUTION: Check your units!")
        print("  - If target is in nm, convert to microns (divide by 1000)")
        print("  - If target is in microns, check wavelength array units")
        print()
        return
    else:
        print(f"✓ Target is in range [{wave_min:.6f}, {wave_max:.6f}]")
        print()
    
    # Step 4: Test injection
    print("STEP 4: Test Injection")
    print("-" * 80)
    print()
    
    flux = np.ones(len(wavelengths)) * 0.08  # Baseline like user's data
    
    try:
        result = add_gaussian_to_array(
            flux,
            fwhm=3.0,
            amplitude=0.05,
            wavelength=target,
            wavelengths=wavelengths
        )
        
        peak_idx = np.argmax(result)
        peak_wave = wavelengths[peak_idx]
        
        print("✓ Injection SUCCESSFUL!")
        print()
        print(f"Peak location:")
        print(f"  Index: {peak_idx} out of {len(wavelengths)-1}")
        print(f"  Wavelength: {peak_wave:.6f}")
        print(f"  Target was: {target}")
        print(f"  Error: {abs(peak_wave - target) * 1000:.2f} nm")
        print()
        
        # Check where peak is
        if peak_idx < 10:
            print("⚠️  PEAK IS AT LEFT EDGE!")
            print(f"   Index {peak_idx} is very close to start (0)")
            print()
            print("This is likely due to unit mismatch or incorrect wavelength array.")
        elif peak_idx > len(wavelengths) - 11:
            print("⚠️  PEAK IS AT RIGHT EDGE!")
            print(f"   Index {peak_idx} is very close to end ({len(wavelengths)-1})")
            print()
            print("This is likely due to unit mismatch or incorrect wavelength array.")
        else:
            print("✓ Peak is in middle of spectrum (correct!)")
            print(f"  Position: {peak_idx}/{len(wavelengths)-1} ({100*peak_idx/(len(wavelengths)-1):.1f}% through array)")
        
        print()
        print("Values around peak:")
        start = max(0, peak_idx - 3)
        end = min(len(wavelengths), peak_idx + 4)
        for i in range(start, end):
            marker = " ← PEAK" if i == peak_idx else ""
            print(f"  [{i:3d}] {wavelengths[i]:.6f}: {result[i]:.6f}{marker}")
        
    except Exception as e:
        print(f"✗ Injection FAILED!")
        print(f"  Error: {e}")
        print()
        import traceback
        traceback.print_exc()
    
    print()
    print("=" * 80)
    print("DIAGNOSTIC COMPLETE")
    print("=" * 80)
    print()

if __name__ == "__main__":
    diagnose_injection_issue()
