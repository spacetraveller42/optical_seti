#!/usr/bin/env python
"""
Demonstration: Setting FWHM and getting Gaussian as an array for direct addition

This script clearly demonstrates that:
1. You CAN set the FWHM directly
2. You GET a Gaussian as a numpy array
3. You CAN add it directly to existing arrays
"""

import sys
from pathlib import Path
cwd = Path(__file__).parent.resolve() 
optical_seti_dir = cwd.parent.resolve() 
sys.path.append(str(optical_seti_dir))

import numpy as np
import matplotlib.pyplot as plt
import optical_seti_functions

print("=" * 70)
print("DEMONSTRATION: FWHM Setting and Array-Based Gaussian Addition")
print("=" * 70)

# Create some existing spectral data
existing_data = np.random.randn(300) * 8 + 500
print(f"\n1. Created existing data array with {len(existing_data)} points")
print(f"   Range: [{np.min(existing_data):.2f}, {np.max(existing_data):.2f}]")

# REQUIREMENT 1: Set the FWHM itself
print("\n2. Setting FWHM directly:")
my_fwhm = 12.5  # You can set any FWHM value you want
my_amplitude = 200.0
my_center = 150.0

# Generate Gaussian with the FWHM you specified
gaussian = optical_seti_functions.generate_gaussian(
    fwhm=my_fwhm,           # ✓ FWHM is set directly
    amplitude=my_amplitude, 
    center=my_center, 
    array_length=len(existing_data)
)

print(f"   ✓ Set FWHM = {my_fwhm}")
print(f"   ✓ Generated Gaussian as numpy array")
print(f"   ✓ Array type: {type(gaussian)}")
print(f"   ✓ Array shape: {gaussian.shape}")

# REQUIREMENT 2: Gaussian given as array that can be added directly
print("\n3. Adding Gaussian array directly to existing data:")
print("   Using standard numpy array addition...")

# Direct array addition (this is what you wanted!)
combined_data = existing_data + gaussian

print(f"   ✓ Successfully added: existing_data + gaussian")
print(f"   ✓ Result type: {type(combined_data)}")
print(f"   ✓ Result shape: {combined_data.shape}")
print(f"   ✓ Peak value increased from {np.max(existing_data):.2f} to {np.max(combined_data):.2f}")

# Show that you can add multiple Gaussians with different FWHM values
print("\n4. Adding multiple Gaussians with different FWHM values:")

# Start fresh
result = existing_data.copy()

# Add first Gaussian with FWHM = 5.0
gaussian1 = optical_seti_functions.generate_gaussian(fwhm=5.0, amplitude=150, center=75, array_length=len(result))
result = result + gaussian1
print(f"   ✓ Added Gaussian 1: FWHM = 5.0")

# Add second Gaussian with FWHM = 20.0
gaussian2 = optical_seti_functions.generate_gaussian(fwhm=20.0, amplitude=100, center=180, array_length=len(result))
result = result + gaussian2
print(f"   ✓ Added Gaussian 2: FWHM = 20.0")

# Add third Gaussian with FWHM = 8.0
gaussian3 = optical_seti_functions.generate_gaussian(fwhm=8.0, amplitude=180, center=240, array_length=len(result))
result = result + gaussian3
print(f"   ✓ Added Gaussian 3: FWHM = 8.0")

print(f"\n   Total: Added 3 Gaussians with different FWHM values to the same array")

# Visualize to confirm
print("\n5. Creating visualization to confirm Gaussians are visible...")

fig, axes = plt.subplots(2, 1, figsize=(14, 10))

# Top plot: Single Gaussian addition
axes[0].plot(existing_data, 'gray', alpha=0.5, linewidth=1, label='Original Data')
axes[0].plot(gaussian, 'r--', linewidth=2, alpha=0.7, label=f'Gaussian (FWHM={my_fwhm})')
axes[0].plot(combined_data, 'b-', linewidth=1.5, label='Original + Gaussian')
axes[0].axvline(x=my_center, color='r', linestyle=':', alpha=0.3)
axes[0].set_title(f'Single Gaussian Addition (FWHM={my_fwhm})', fontsize=14, fontweight='bold')
axes[0].set_xlabel('Array Index')
axes[0].set_ylabel('Value')
axes[0].legend(loc='upper right')
axes[0].grid(True, alpha=0.3)

# Bottom plot: Multiple Gaussians
axes[1].plot(existing_data, 'gray', alpha=0.5, linewidth=1, label='Original Data')
axes[1].plot(result, 'b-', linewidth=1.5, label='With 3 Gaussians (FWHM: 5.0, 20.0, 8.0)')
axes[1].axvline(x=75, color='r', linestyle=':', alpha=0.3, label='Gaussian centers')
axes[1].axvline(x=180, color='r', linestyle=':', alpha=0.3)
axes[1].axvline(x=240, color='r', linestyle=':', alpha=0.3)
axes[1].set_title('Multiple Gaussians with Different FWHM Values', fontsize=14, fontweight='bold')
axes[1].set_xlabel('Array Index')
axes[1].set_ylabel('Value')
axes[1].legend(loc='upper right')
axes[1].grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig('demonstration_fwhm_array_addition.png', dpi=150, bbox_inches='tight')
print("   ✓ Visualization saved as 'demonstration_fwhm_array_addition.png'")

print("\n" + "=" * 70)
print("SUMMARY:")
print("=" * 70)
print("✓ You CAN set the FWHM directly (parameter: fwhm)")
print("✓ You GET a numpy array from generate_gaussian()")
print("✓ You CAN add it directly to existing arrays with +")
print("✓ All Gaussians are clearly visible in the plots")
print("=" * 70)

plt.show()
