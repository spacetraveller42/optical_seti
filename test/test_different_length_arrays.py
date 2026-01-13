# Test script for adding Gaussians to arrays of different lengths
# This tests the new array_length parameter in add_gaussian_to_array

import sys
from pathlib import Path
cwd = Path(__file__).parent.resolve() 
optical_seti_dir = cwd.parent.resolve() 
sys.path.append(str(optical_seti_dir))

import numpy as np
import matplotlib.pyplot as plt
import optical_seti_functions

print("=" * 70)
print("TEST: Adding Gaussians to Arrays of Different Lengths")
print("=" * 70)

# Test 1: Gaussian shorter than data array (will be zero-padded)
print("\nTest 1: Shorter Gaussian (50 elements) added to longer data (200 elements)")
data_long = np.ones(200) * 100.0
result1 = optical_seti_functions.add_gaussian_to_array(
    data_long, 
    fwhm=10.0, 
    amplitude=200.0, 
    center=25.0, 
    array_length=50  # Gaussian is only 50 elements
)

print(f"   Data length: {len(data_long)}")
print(f"   Gaussian length specified: 50")
print(f"   Result length: {len(result1)}")
print(f"   Peak in first 50 elements: {np.max(result1[:50]):.2f}")
print(f"   Values after 50th element (should be ~100): {result1[60]:.2f}")
assert len(result1) == 200, "Result should match data length"
assert result1[60] == 100.0, "Zero-padded region should have only baseline"
print("   ✓ PASSED: Shorter Gaussian correctly zero-padded")

# Test 2: Gaussian longer than data array (will be truncated)
print("\nTest 2: Longer Gaussian (200 elements) added to shorter data (50 elements)")
data_short = np.ones(50) * 100.0
result2 = optical_seti_functions.add_gaussian_to_array(
    data_short, 
    fwhm=10.0, 
    amplitude=200.0, 
    center=25.0, 
    array_length=200  # Gaussian is 200 elements
)

print(f"   Data length: {len(data_short)}")
print(f"   Gaussian length specified: 200")
print(f"   Result length: {len(result2)}")
print(f"   Peak value: {np.max(result2):.2f}")
assert len(result2) == 50, "Result should match data length"
print("   ✓ PASSED: Longer Gaussian correctly truncated")

# Test 3: Same length (default behavior - should work as before)
print("\nTest 3: Same length (default behavior, no array_length specified)")
data_same = np.ones(100) * 100.0
result3 = optical_seti_functions.add_gaussian_to_array(
    data_same, 
    fwhm=10.0, 
    amplitude=200.0, 
    center=50.0
)

print(f"   Data length: {len(data_same)}")
print(f"   Result length: {len(result3)}")
print(f"   Peak value: {np.max(result3):.2f}")
assert len(result3) == 100, "Result should match data length"
print("   ✓ PASSED: Default behavior (same length) works")

# Test 4: Visual demonstration
print("\nTest 4: Creating visualization...")

fig, axes = plt.subplots(3, 1, figsize=(12, 10))

# Plot 1: Shorter Gaussian (zero-padded)
axes[0].plot(data_long, 'gray', alpha=0.5, linewidth=1, label='Original data (200 elements)')
axes[0].plot(result1, 'b-', linewidth=2, label='With short Gaussian (50 elem, zero-padded)')
axes[0].axvline(x=50, color='r', linestyle='--', alpha=0.5, label='Gaussian boundary')
axes[0].set_title('Shorter Gaussian Added to Longer Array (Zero-Padded)', fontsize=12, fontweight='bold')
axes[0].set_xlabel('Array Index')
axes[0].set_ylabel('Value')
axes[0].legend()
axes[0].grid(True, alpha=0.3)

# Plot 2: Longer Gaussian (truncated)
axes[1].plot(data_short, 'gray', alpha=0.5, linewidth=1, label='Original data (50 elements)')
axes[1].plot(result2, 'b-', linewidth=2, label='With long Gaussian (200 elem, truncated)')
axes[1].set_title('Longer Gaussian Added to Shorter Array (Truncated)', fontsize=12, fontweight='bold')
axes[1].set_xlabel('Array Index')
axes[1].set_ylabel('Value')
axes[1].legend()
axes[1].grid(True, alpha=0.3)

# Plot 3: Same length (default)
axes[2].plot(data_same, 'gray', alpha=0.5, linewidth=1, label='Original data (100 elements)')
axes[2].plot(result3, 'b-', linewidth=2, label='With same-length Gaussian')
axes[2].set_title('Same Length Gaussian (Default Behavior)', fontsize=12, fontweight='bold')
axes[2].set_xlabel('Array Index')
axes[2].set_ylabel('Value')
axes[2].legend()
axes[2].grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig('test_different_length_gaussians.png', dpi=150, bbox_inches='tight')
print("   ✓ Visualization saved as 'test_different_length_gaussians.png'")

print("\n" + "=" * 70)
print("ALL TESTS PASSED ✓")
print("=" * 70)
print("\nGaussians can now be added to arrays of different lengths!")
print("- Shorter Gaussians are zero-padded")
print("- Longer Gaussians are truncated")
print("- Default behavior (same length) unchanged")
print("=" * 70)

plt.show()
