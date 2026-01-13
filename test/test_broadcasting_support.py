# Test script for broadcasting support in add_gaussian_to_array
# This tests adding Gaussians to multi-dimensional arrays

import sys
from pathlib import Path
cwd = Path(__file__).parent.resolve() 
optical_seti_dir = cwd.parent.resolve() 
sys.path.append(str(optical_seti_dir))

import numpy as np
import matplotlib.pyplot as plt
import optical_seti_functions

print("=" * 70)
print("TEST: Broadcasting Support for Multi-Dimensional Arrays")
print("=" * 70)

# Test 1: Add Gaussian to 2D array along axis 1 (columns)
print("\nTest 1: 2D array - Add Gaussian along axis 1 (broadcasts across rows)")
data_2d = np.ones((50, 100)) * 100.0
print(f"   Input shape: {data_2d.shape}")

result_2d = optical_seti_functions.add_gaussian_to_array(
    data_2d, 
    fwhm=10.0, 
    amplitude=200.0, 
    center=50.0, 
    axis=1
)

print(f"   Output shape: {result_2d.shape}")
print(f"   Peak value: {np.max(result_2d):.2f}")
print(f"   All rows have same pattern: {np.allclose(result_2d[0], result_2d[25])}")
assert result_2d.shape == (50, 100), "Shape should be preserved"
assert np.allclose(result_2d[0], result_2d[25]), "All rows should have same Gaussian"
print("   ✓ PASSED: Gaussian correctly added along axis 1")

# Test 2: Add Gaussian to 2D array along axis 0 (rows)
print("\nTest 2: 2D array - Add Gaussian along axis 0 (broadcasts across columns)")
data_2d = np.ones((100, 50)) * 100.0
print(f"   Input shape: {data_2d.shape}")

result_2d = optical_seti_functions.add_gaussian_to_array(
    data_2d, 
    fwhm=15.0, 
    amplitude=150.0, 
    center=50.0, 
    axis=0
)

print(f"   Output shape: {result_2d.shape}")
print(f"   Peak value: {np.max(result_2d):.2f}")
print(f"   All columns have same pattern: {np.allclose(result_2d[:, 0], result_2d[:, 25])}")
assert result_2d.shape == (100, 50), "Shape should be preserved"
assert np.allclose(result_2d[:, 0], result_2d[:, 25]), "All columns should have same Gaussian"
print("   ✓ PASSED: Gaussian correctly added along axis 0")

# Test 3: Add Gaussian to 3D array
print("\nTest 3: 3D array - Add Gaussian along axis 2 (broadcasts across first two dims)")
data_3d = np.ones((10, 20, 80)) * 100.0
print(f"   Input shape: {data_3d.shape}")

result_3d = optical_seti_functions.add_gaussian_to_array(
    data_3d, 
    fwhm=8.0, 
    amplitude=100.0, 
    center=40.0, 
    axis=2
)

print(f"   Output shape: {result_3d.shape}")
print(f"   Peak value: {np.max(result_3d):.2f}")
assert result_3d.shape == (10, 20, 80), "Shape should be preserved"
print("   ✓ PASSED: Gaussian correctly added to 3D array")

# Test 4: Backward compatibility - 1D array still works
print("\nTest 4: Backward compatibility - 1D array")
data_1d = np.ones(100) * 100.0
result_1d = optical_seti_functions.add_gaussian_to_array(
    data_1d, 
    fwhm=10.0, 
    amplitude=200.0, 
    center=50.0
)

print(f"   Input shape: {data_1d.shape}")
print(f"   Output shape: {result_1d.shape}")
print(f"   Peak value: {np.max(result_1d):.2f}")
assert result_1d.shape == (100,), "Shape should be preserved"
assert np.max(result_1d) == 300.0, "Peak should be baseline + amplitude"
print("   ✓ PASSED: 1D arrays still work (backward compatible)")

# Test 5: Different array_length with broadcasting
print("\nTest 5: Different array_length with 2D array")
data_2d = np.ones((30, 200)) * 100.0
result_2d = optical_seti_functions.add_gaussian_to_array(
    data_2d, 
    fwhm=10.0, 
    amplitude=150.0, 
    center=25.0,
    array_length=50,  # Shorter Gaussian, will be zero-padded
    axis=1
)

print(f"   Input shape: {data_2d.shape}")
print(f"   Output shape: {result_2d.shape}")
print(f"   Gaussian length specified: 50 (zero-padded to 200)")
print(f"   Peak in first 50 elements: {np.max(result_2d[:, :50]):.2f}")
print(f"   Values after 50th element: {result_2d[0, 60]:.2f} (should be ~100)")
assert result_2d.shape == (30, 200), "Shape should be preserved"
assert result_2d[0, 60] == 100.0, "Zero-padded region should have only baseline"
print("   ✓ PASSED: Different array_length works with broadcasting")

# Test 6: Visualization
print("\nTest 6: Creating visualizations...")

fig, axes = plt.subplots(2, 2, figsize=(14, 10))

# Plot 1: 2D array with Gaussian along axis 1
data_2d = np.random.randn(40, 120) * 10 + 500
result_2d = optical_seti_functions.add_gaussian_to_array(
    data_2d, fwhm=12.0, amplitude=300.0, center=60.0, axis=1
)
im1 = axes[0, 0].imshow(result_2d, aspect='auto', cmap='viridis')
axes[0, 0].set_title('2D Array: Gaussian Along Axis 1 (Columns)', fontweight='bold')
axes[0, 0].set_xlabel('Column Index')
axes[0, 0].set_ylabel('Row Index')
plt.colorbar(im1, ax=axes[0, 0])

# Plot 2: Cross-section showing broadcasting
axes[0, 1].plot(result_2d[0, :], 'b-', linewidth=2, label='Row 0')
axes[0, 1].plot(result_2d[20, :], 'r--', linewidth=2, label='Row 20')
axes[0, 1].plot(result_2d[39, :], 'g:', linewidth=2, label='Row 39')
axes[0, 1].set_title('Cross-Sections Showing Broadcasting', fontweight='bold')
axes[0, 1].set_xlabel('Column Index')
axes[0, 1].set_ylabel('Value')
axes[0, 1].legend()
axes[0, 1].grid(True, alpha=0.3)

# Plot 3: 2D array with Gaussian along axis 0
data_2d_ax0 = np.random.randn(120, 40) * 10 + 500
result_2d_ax0 = optical_seti_functions.add_gaussian_to_array(
    data_2d_ax0, fwhm=15.0, amplitude=250.0, center=60.0, axis=0
)
im2 = axes[1, 0].imshow(result_2d_ax0, aspect='auto', cmap='plasma')
axes[1, 0].set_title('2D Array: Gaussian Along Axis 0 (Rows)', fontweight='bold')
axes[1, 0].set_xlabel('Column Index')
axes[1, 0].set_ylabel('Row Index')
plt.colorbar(im2, ax=axes[1, 0])

# Plot 4: Cross-section for axis 0
axes[1, 1].plot(result_2d_ax0[:, 0], 'b-', linewidth=2, label='Column 0')
axes[1, 1].plot(result_2d_ax0[:, 20], 'r--', linewidth=2, label='Column 20')
axes[1, 1].plot(result_2d_ax0[:, 39], 'g:', linewidth=2, label='Column 39')
axes[1, 1].set_title('Cross-Sections for Axis 0 Broadcasting', fontweight='bold')
axes[1, 1].set_xlabel('Row Index')
axes[1, 1].set_ylabel('Value')
axes[1, 1].legend()
axes[1, 1].grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig('test_broadcasting_support.png', dpi=150, bbox_inches='tight')
print("   ✓ Visualization saved as 'test_broadcasting_support.png'")

print("\n" + "=" * 70)
print("ALL TESTS PASSED ✓")
print("=" * 70)
print("\nBroadcasting feature summary:")
print("- ✓ Add Gaussians to 2D arrays along any axis")
print("- ✓ Add Gaussians to 3D and higher dimensional arrays")
print("- ✓ Broadcasting works with different array_length")
print("- ✓ Backward compatible with 1D arrays")
print("- ✓ Gaussian automatically broadcasts across other dimensions")
print("=" * 70)

plt.show()
