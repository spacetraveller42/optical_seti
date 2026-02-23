#!/usr/bin/env python
"""
Visual demonstration showing why laser appears at edge.
Creates a plot comparing center=0.532 vs wavelength=0.532
"""

import sys
import os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import numpy as np
import matplotlib.pyplot as plt
from gaussian_generator import add_gaussian_to_array

# User's wavelength array
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

# Baseline flux
flux = np.ones(len(wavelengths)) * 0.08

# Test 1: WRONG way (center=0.532)
result_wrong = add_gaussian_to_array(flux, fwhm=3.0, amplitude=0.10, center=0.532)

# Test 2: CORRECT way (wavelength=0.532)
result_correct = add_gaussian_to_array(flux, fwhm=3.0, amplitude=0.10, 
                                       wavelength=0.532, wavelengths=wavelengths)

# Create visualization
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 10))

# Plot 1: WRONG (center=0.532)
ax1.plot(wavelengths, result_wrong, 'r-', linewidth=2, label='With Gaussian')
ax1.plot(wavelengths, flux, 'gray', alpha=0.5, label='Baseline', linestyle='--')
ax1.axvline(x=0.532, color='green', linestyle=':', linewidth=2, alpha=0.7, 
            label='Target: 0.532 μm')
peak_idx_wrong = np.argmax(result_wrong)
ax1.axvline(x=wavelengths[peak_idx_wrong], color='red', linestyle='--', linewidth=2,
            label=f'Actual peak: {wavelengths[peak_idx_wrong]:.3f} μm')
ax1.set_xlabel('Wavelength (microns)', fontsize=12)
ax1.set_ylabel('Flux', fontsize=12)
ax1.set_title('❌ WRONG: center=0.532 (treats 0.532 as ARRAY INDEX)', 
              fontsize=14, fontweight='bold', color='red')
ax1.legend(fontsize=10)
ax1.grid(True, alpha=0.3)
ax1.set_xlim(0.15, 2.05)

# Add annotation
ax1.annotate('Laser appears\nhere at EDGE!', 
            xy=(wavelengths[peak_idx_wrong], result_wrong[peak_idx_wrong]),
            xytext=(0.4, 0.16),
            fontsize=12, fontweight='bold', color='red',
            arrowprops=dict(arrowstyle='->', color='red', lw=2))
ax1.annotate('You wanted it\nhere!', 
            xy=(0.532, 0.08),
            xytext=(0.7, 0.09),
            fontsize=11, color='green',
            arrowprops=dict(arrowstyle='->', color='green', lw=1.5))

# Plot 2: CORRECT (wavelength=0.532)
ax2.plot(wavelengths, result_correct, 'g-', linewidth=2, label='With Gaussian')
ax2.plot(wavelengths, flux, 'gray', alpha=0.5, label='Baseline', linestyle='--')
ax2.axvline(x=0.532, color='green', linestyle=':', linewidth=2, alpha=0.7, 
            label='Target: 0.532 μm')
peak_idx_correct = np.argmax(result_correct)
ax2.axvline(x=wavelengths[peak_idx_correct], color='darkgreen', linestyle='--', linewidth=2,
            label=f'Actual peak: {wavelengths[peak_idx_correct]:.3f} μm')
ax2.set_xlabel('Wavelength (microns)', fontsize=12)
ax2.set_ylabel('Flux', fontsize=12)
ax2.set_title('✓ CORRECT: wavelength=0.532, wavelengths=array (treats 0.532 as WAVELENGTH)', 
              fontsize=14, fontweight='bold', color='green')
ax2.legend(fontsize=10)
ax2.grid(True, alpha=0.3)
ax2.set_xlim(0.15, 2.05)

# Add annotation
ax2.annotate('Laser appears\nexactly where\nyou want it!', 
            xy=(wavelengths[peak_idx_correct], result_correct[peak_idx_correct]),
            xytext=(0.7, 0.16),
            fontsize=12, fontweight='bold', color='darkgreen',
            arrowprops=dict(arrowstyle='->', color='darkgreen', lw=2))

plt.tight_layout()
plt.savefig('/tmp/edge_injection_comparison.png', dpi=150, bbox_inches='tight')
print("Visualization saved to /tmp/edge_injection_comparison.png")
print()
print("Summary:")
print("-" * 60)
print(f"WRONG (center=0.532):    Peak at {wavelengths[peak_idx_wrong]:.6f} μm (index {peak_idx_wrong})")
print(f"CORRECT (wavelength=0.532): Peak at {wavelengths[peak_idx_correct]:.6f} μm (index {peak_idx_correct})")
print()
print("The difference is dramatic!")
plt.show()
