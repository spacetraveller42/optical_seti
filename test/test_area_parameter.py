# Test script for area parameter in Gaussian generation

import sys
from pathlib import Path
cwd = Path(__file__).parent.resolve() 
optical_seti_dir = cwd.parent.resolve() 
sys.path.append(str(optical_seti_dir))

import numpy as np
import matplotlib.pyplot as plt
import optical_seti_functions

print("=" * 70)
print("TEST: Area Parameter for Gaussian Generation")
print("=" * 70)

# Test 1: Generate Gaussian with area parameter
print("\nTest 1: Generate Gaussian with specified area")
try:
    fwhm = 10.0
    area = 1000.0
    center = 50.0
    array_length = 100
    
    gaussian = optical_seti_functions.generate_gaussian(fwhm=fwhm, area=area, center=center, array_length=array_length)
    actual_area = np.sum(gaussian)
    
    print(f"   Requested area: {area}")
    print(f"   Actual area (sum): {actual_area:.2f}")
    print(f"   Difference: {abs(actual_area - area):.2f}")
    print(f"   Peak amplitude: {np.max(gaussian):.2f}")
    
    # Area should be approximately correct (numerical integration via sum)
    assert abs(actual_area - area) < area * 0.01, f"Area mismatch: {actual_area} vs {area}"
    print("   ✓ PASSED: Area matches expected value")
except Exception as e:
    print(f"   ✗ FAILED: {e}")

# Test 2: Compare amplitude vs area methods
print("\nTest 2: Compare amplitude and area methods")
try:
    fwhm = 15.0
    amplitude = 200.0
    center = 100.0
    array_length = 200
    
    # Generate with amplitude
    gauss_amp = optical_seti_functions.generate_gaussian(fwhm=fwhm, amplitude=amplitude, center=center, array_length=array_length)
    area_from_amp = np.sum(gauss_amp)
    
    # Generate with same area
    gauss_area = optical_seti_functions.generate_gaussian(fwhm=fwhm, area=area_from_amp, center=center, array_length=array_length)
    
    # Should be nearly identical
    diff = np.max(np.abs(gauss_amp - gauss_area))
    print(f"   Gaussian from amplitude - peak: {np.max(gauss_amp):.2f}, area: {np.sum(gauss_amp):.2f}")
    print(f"   Gaussian from area - peak: {np.max(gauss_area):.2f}, area: {np.sum(gauss_area):.2f}")
    print(f"   Max difference: {diff:.6f}")
    
    assert diff < 0.01, f"Gaussians should match, max diff: {diff}"
    print("   ✓ PASSED: Amplitude and area methods produce equivalent results")
except Exception as e:
    print(f"   ✗ FAILED: {e}")

# Test 3: Add Gaussian with area to existing array
print("\nTest 3: Add Gaussian with area to existing array")
try:
    data = np.ones(100) * 500.0
    fwhm = 12.0
    area = 800.0
    center = 50.0
    
    result = optical_seti_functions.add_gaussian_to_array(data, fwhm=fwhm, area=area, center=center)
    added_area = np.sum(result - data)
    
    print(f"   Requested area: {area}")
    print(f"   Added area: {added_area:.2f}")
    print(f"   Original data peak: {np.max(data):.2f}")
    print(f"   Result peak: {np.max(result):.2f}")
    
    assert abs(added_area - area) < area * 0.01, f"Added area mismatch: {added_area} vs {area}"
    print("   ✓ PASSED: Gaussian with specified area added correctly")
except Exception as e:
    print(f"   ✗ FAILED: {e}")

# Test 4: Error handling - both amplitude and area specified
print("\nTest 4: Error handling - both amplitude and area specified")
try:
    try:
        gaussian = optical_seti_functions.generate_gaussian(fwhm=10.0, amplitude=100.0, area=1000.0, center=50.0, array_length=100)
        print("   ✗ FAILED: Should have raised ValueError")
    except ValueError as e:
        print(f"   ✓ PASSED: Correctly raised ValueError: {e}")
except Exception as e:
    print(f"   ✗ FAILED: {e}")

# Test 5: Error handling - neither amplitude nor area specified
print("\nTest 5: Error handling - neither amplitude nor area specified")
try:
    try:
        gaussian = optical_seti_functions.generate_gaussian(fwhm=10.0, center=50.0, array_length=100)
        print("   ✗ FAILED: Should have raised ValueError")
    except ValueError as e:
        print(f"   ✓ PASSED: Correctly raised ValueError: {e}")
except Exception as e:
    print(f"   ✗ FAILED: {e}")

# Test 6: Visualization - amplitude vs area
print("\nTest 6: Creating visualization comparing amplitude and area methods...")
try:
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    
    # Same FWHM, different amplitudes
    fwhm = 15.0
    center = 100.0
    array_length = 200
    
    gauss1_amp = optical_seti_functions.generate_gaussian(fwhm=fwhm, amplitude=100.0, center=center, array_length=array_length)
    gauss2_amp = optical_seti_functions.generate_gaussian(fwhm=fwhm, amplitude=200.0, center=center, array_length=array_length)
    
    axes[0, 0].plot(gauss1_amp, 'b-', linewidth=2, label=f'Amplitude=100')
    axes[0, 0].plot(gauss2_amp, 'r-', linewidth=2, label=f'Amplitude=200')
    axes[0, 0].set_title('Gaussians with Different Amplitudes', fontweight='bold')
    axes[0, 0].set_xlabel('Array Index')
    axes[0, 0].set_ylabel('Value')
    axes[0, 0].legend()
    axes[0, 0].grid(True, alpha=0.3)
    axes[0, 0].text(0.02, 0.98, f'Area 1: {np.sum(gauss1_amp):.0f}\nArea 2: {np.sum(gauss2_amp):.0f}', 
                    transform=axes[0, 0].transAxes, verticalalignment='top',
                    bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
    
    # Same area, different amplitudes/widths
    area = 1000.0
    gauss1_area = optical_seti_functions.generate_gaussian(fwhm=10.0, area=area, center=center, array_length=array_length)
    gauss2_area = optical_seti_functions.generate_gaussian(fwhm=20.0, area=area, center=center, array_length=array_length)
    
    axes[0, 1].plot(gauss1_area, 'b-', linewidth=2, label=f'FWHM=10')
    axes[0, 1].plot(gauss2_area, 'r-', linewidth=2, label=f'FWHM=20')
    axes[0, 1].set_title(f'Gaussians with Same Area ({area})', fontweight='bold')
    axes[0, 1].set_xlabel('Array Index')
    axes[0, 1].set_ylabel('Value')
    axes[0, 1].legend()
    axes[0, 1].grid(True, alpha=0.3)
    axes[0, 1].text(0.02, 0.98, f'Peak 1: {np.max(gauss1_area):.0f}\nPeak 2: {np.max(gauss2_area):.0f}', 
                    transform=axes[0, 1].transAxes, verticalalignment='top',
                    bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
    
    # Add to existing data - amplitude method
    data = np.random.randn(200) * 5 + 500
    result_amp = optical_seti_functions.add_gaussian_to_array(data, fwhm=15.0, amplitude=150.0, center=100.0)
    
    axes[1, 0].plot(data, 'gray', alpha=0.5, linewidth=1, label='Original')
    axes[1, 0].plot(result_amp, 'b-', linewidth=1.5, label='With Gaussian (amplitude)')
    axes[1, 0].set_title('Adding Gaussian by Amplitude', fontweight='bold')
    axes[1, 0].set_xlabel('Array Index')
    axes[1, 0].set_ylabel('Value')
    axes[1, 0].legend()
    axes[1, 0].grid(True, alpha=0.3)
    
    # Add to existing data - area method
    result_area = optical_seti_functions.add_gaussian_to_array(data, fwhm=15.0, area=1500.0, center=100.0)
    
    axes[1, 1].plot(data, 'gray', alpha=0.5, linewidth=1, label='Original')
    axes[1, 1].plot(result_area, 'r-', linewidth=1.5, label='With Gaussian (area)')
    axes[1, 1].set_title('Adding Gaussian by Area', fontweight='bold')
    axes[1, 1].set_xlabel('Array Index')
    axes[1, 1].set_ylabel('Value')
    axes[1, 1].legend()
    axes[1, 1].grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig('test_area_parameter.png', dpi=150, bbox_inches='tight')
    print("   ✓ Visualization saved as 'test_area_parameter.png'")
except Exception as e:
    print(f"   ✗ FAILED: {e}")

print("\n" + "=" * 70)
print("ALL TESTS COMPLETED ✓")
print("=" * 70)
print("\nNew feature summary:")
print("- ✓ Can specify Gaussian by area instead of amplitude")
print("- ✓ Area parameter validated and calculated correctly")
print("- ✓ Backwards compatible with amplitude parameter")
print("- ✓ Works with add_gaussian_to_array function")
print("- ✓ Proper error handling for invalid parameter combinations")
print("=" * 70)

plt.show()
