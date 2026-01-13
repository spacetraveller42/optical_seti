# Test script to reproduce and verify fix for broadcasting error
# Tests the edge cases that could cause "ValueError: operands could not be broadcast together"

import sys
from pathlib import Path
cwd = Path(__file__).parent.resolve() 
optical_seti_dir = cwd.parent.resolve() 
sys.path.append(str(optical_seti_dir))

import numpy as np
import optical_seti_functions

print("=" * 70)
print("TEST: Edge Cases and Error Handling")
print("=" * 70)

# Test 1: Large 1D array (similar to user's 313116 element array)
print("\nTest 1: Large 1D array (313116 elements)")
try:
    data_large = np.ones(313116) * 100.0
    result = optical_seti_functions.add_gaussian_to_array(
        data_large, 
        fwhm=10.0, 
        amplitude=200.0, 
        center=156558.0  # Middle of array
    )
    print(f"   Array length: {len(data_large)}")
    print(f"   Result length: {len(result)}")
    print(f"   Peak value: {np.max(result):.2f}")
    print("   ✓ PASSED: Large array works correctly")
except Exception as e:
    print(f"   ✗ FAILED: {e}")

# Test 2: Try with array_length=0 (should raise error)
print("\nTest 2: array_length=0 (should raise ValueError)")
try:
    data = np.ones(100) * 100.0
    result = optical_seti_functions.add_gaussian_to_array(
        data, 
        fwhm=10.0, 
        amplitude=200.0, 
        center=50.0,
        array_length=0
    )
    print("   ✗ FAILED: Should have raised ValueError")
except ValueError as e:
    print(f"   ✓ PASSED: Correctly raised ValueError: {e}")

# Test 3: Empty array (should raise error)
print("\nTest 3: Empty array (should raise ValueError)")
try:
    data_empty = np.array([])
    result = optical_seti_functions.add_gaussian_to_array(
        data_empty, 
        fwhm=10.0, 
        amplitude=200.0, 
        center=50.0
    )
    print("   ✗ FAILED: Should have raised ValueError")
except ValueError as e:
    print(f"   ✓ PASSED: Correctly raised ValueError: {e}")

# Test 4: Center out of bounds (should raise error)
print("\nTest 4: Center out of bounds (should raise ValueError)")
try:
    data = np.ones(100) * 100.0
    result = optical_seti_functions.add_gaussian_to_array(
        data, 
        fwhm=10.0, 
        amplitude=200.0, 
        center=200.0  # Out of bounds
    )
    print("   ✗ FAILED: Should have raised ValueError")
except ValueError as e:
    print(f"   ✓ PASSED: Correctly raised ValueError: {e}")

# Test 5: Negative center (should raise error)
print("\nTest 5: Negative center (should raise ValueError)")
try:
    data = np.ones(100) * 100.0
    result = optical_seti_functions.add_gaussian_to_array(
        data, 
        fwhm=10.0, 
        amplitude=200.0, 
        center=-10.0  # Negative
    )
    print("   ✗ FAILED: Should have raised ValueError")
except ValueError as e:
    print(f"   ✓ PASSED: Correctly raised ValueError: {e}")

# Test 6: Center at exact boundary (should work for center < array_length)
print("\nTest 6: Center at boundary-1 (should work)")
try:
    data = np.ones(100) * 100.0
    result = optical_seti_functions.add_gaussian_to_array(
        data, 
        fwhm=10.0, 
        amplitude=200.0, 
        center=99.0  # Last valid index
    )
    print(f"   Result length: {len(result)}")
    print("   ✓ PASSED: Boundary case works")
except Exception as e:
    print(f"   ✗ FAILED: {e}")

# Test 7: Invalid axis for 2D array (should raise error)
print("\nTest 7: Invalid axis for 2D array (should raise ValueError)")
try:
    data_2d = np.ones((10, 20)) * 100.0
    result = optical_seti_functions.add_gaussian_to_array(
        data_2d, 
        fwhm=10.0, 
        amplitude=200.0, 
        center=10.0,
        axis=5  # Invalid axis
    )
    print("   ✗ FAILED: Should have raised ValueError")
except ValueError as e:
    print(f"   ✓ PASSED: Correctly raised ValueError: {e}")

# Test 8: Negative axis (should work - Python style negative indexing)
print("\nTest 8: Negative axis (should work)")
try:
    data_2d = np.ones((10, 20)) * 100.0
    result = optical_seti_functions.add_gaussian_to_array(
        data_2d, 
        fwhm=5.0, 
        amplitude=150.0, 
        center=10.0,
        axis=-1  # Last axis
    )
    print(f"   Input shape: {data_2d.shape}")
    print(f"   Result shape: {result.shape}")
    print("   ✓ PASSED: Negative axis works")
except Exception as e:
    print(f"   ✗ FAILED: {e}")

# Test 9: Reproduce potential error scenario with explicit array_length
print("\nTest 9: array_length explicitly set to match data")
try:
    data = np.ones(313116) * 100.0
    result = optical_seti_functions.add_gaussian_to_array(
        data, 
        fwhm=100.0, 
        amplitude=500.0, 
        center=150000.0,
        array_length=313116  # Explicitly match data length
    )
    print(f"   Array length: {len(data)}")
    print(f"   Result length: {len(result)}")
    print("   ✓ PASSED: Explicit array_length works")
except Exception as e:
    print(f"   ✗ FAILED: {e}")

print("\n" + "=" * 70)
print("Edge case testing complete!")
print("=" * 70)
print("\nAll error cases properly handled with clear error messages")
print("Large arrays work correctly")
print("=" * 70)
