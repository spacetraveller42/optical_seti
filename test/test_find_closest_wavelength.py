#!/usr/bin/env python
"""
Test suite for find_closest_wavelength_index function.

Tests the new utility function that finds the closest wavelength in an array
to a specified target wavelength.
"""

import numpy as np
import sys
import os

# Add parent directory to path
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from gaussian_generator import find_closest_wavelength_index


def test_basic_functionality():
    """Test basic functionality with simple array."""
    print("\n" + "="*70)
    print("TEST 1: Basic Functionality")
    print("="*70)
    
    # Simple wavelength array
    wavelengths = np.array([0.5, 0.6, 0.7, 0.8, 0.9, 1.0])
    target = 0.75
    
    idx, wl, diff = find_closest_wavelength_index(target, wavelengths)
    
    print(f"Wavelengths: {wavelengths}")
    print(f"Target: {target}")
    print(f"Result: index={idx}, wavelength={wl}, difference={diff}")
    
    # Should find index 2 (0.7) or 3 (0.8), depending on which is closer
    assert idx in [2, 3], f"Expected index 2 or 3, got {idx}"
    assert abs(wl - target) == diff, "Difference calculation incorrect"
    
    print("✓ PASSED\n")
    return True


def test_exact_match():
    """Test when target exactly matches an array element."""
    print("\n" + "="*70)
    print("TEST 2: Exact Match")
    print("="*70)
    
    wavelengths = np.array([0.5, 0.6, 0.7, 0.8, 0.9, 1.0])
    target = 0.7  # Exact match
    
    idx, wl, diff = find_closest_wavelength_index(target, wavelengths)
    
    print(f"Wavelengths: {wavelengths}")
    print(f"Target: {target} (exact match)")
    print(f"Result: index={idx}, wavelength={wl}, difference={diff}")
    
    assert idx == 2, f"Expected index 2, got {idx}"
    assert wl == 0.7, f"Expected wavelength 0.7, got {wl}"
    assert diff == 0.0, f"Expected difference 0.0, got {diff}"
    
    print("✓ PASSED - Found exact match\n")
    return True


def test_user_wavelength_array():
    """Test with user's actual 151-point wavelength array."""
    print("\n" + "="*70)
    print("TEST 3: User's Wavelength Array (151 points)")
    print("="*70)
    
    # User's actual wavelength array
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
    
    target = 0.532  # User's green laser wavelength
    
    print(f"Array length: {len(wavelengths)}")
    print(f"Wavelength range: {wavelengths[0]:.6f} to {wavelengths[-1]:.6f} microns")
    print(f"Target wavelength: {target} microns (green laser)")
    
    idx, wl, diff = find_closest_wavelength_index(target, wavelengths)
    
    print(f"\nResult:")
    print(f"  Index: {idx}")
    print(f"  Closest wavelength: {wl:.6f} microns")
    print(f"  Difference: {diff:.6f} microns ({diff*1000:.2f} nm)")
    
    # Verify result is reasonable
    assert 0 <= idx < len(wavelengths), f"Index {idx} out of range"
    assert wavelengths.min() <= wl <= wavelengths.max(), "Wavelength out of range"
    assert diff >= 0, "Difference should be non-negative"
    
    # Should be close to index 13 (0.533657) based on previous analysis
    print(f"\nNearby wavelengths:")
    for i in range(max(0, idx-2), min(len(wavelengths), idx+3)):
        marker = " ← CLOSEST" if i == idx else ""
        print(f"  Index {i}: {wavelengths[i]:.6f} microns{marker}")
    
    print("✓ PASSED\n")
    return True


def test_edge_cases():
    """Test edge cases and error handling."""
    print("\n" + "="*70)
    print("TEST 4: Edge Cases and Error Handling")
    print("="*70)
    
    wavelengths = np.array([0.5, 0.6, 0.7, 0.8, 0.9, 1.0])
    
    # Test 4a: Target at beginning of range
    print("\nTest 4a: Target at beginning of range")
    target = 0.45
    idx, wl, diff = find_closest_wavelength_index(target, wavelengths)
    print(f"Target: {target}, Closest: index {idx}, wavelength {wl}")
    assert idx == 0, f"Expected index 0 for target at start, got {idx}"
    print("✓ Passed")
    
    # Test 4b: Target at end of range
    print("\nTest 4b: Target at end of range")
    target = 1.05
    idx, wl, diff = find_closest_wavelength_index(target, wavelengths)
    print(f"Target: {target}, Closest: index {idx}, wavelength {wl}")
    assert idx == len(wavelengths) - 1, f"Expected last index for target at end, got {idx}"
    print("✓ Passed")
    
    # Test 4c: Empty array should raise error
    print("\nTest 4c: Empty array should raise error")
    try:
        empty_array = np.array([])
        find_closest_wavelength_index(0.5, empty_array)
        print("✗ FAILED - Should have raised ValueError for empty array")
        return False
    except ValueError as e:
        print(f"✓ Correctly raised ValueError: {e}")
    
    # Test 4d: Non-monotonic array should raise error
    print("\nTest 4d: Non-monotonic array should raise error")
    try:
        non_monotonic = np.array([0.5, 0.7, 0.6, 0.8])  # Not monotonic
        find_closest_wavelength_index(0.65, non_monotonic)
        print("✗ FAILED - Should have raised ValueError for non-monotonic array")
        return False
    except ValueError as e:
        print(f"✓ Correctly raised ValueError: {e}")
    
    print("\n✓ ALL EDGE CASES PASSED\n")
    return True


def test_decreasing_wavelengths():
    """Test with decreasing wavelength array."""
    print("\n" + "="*70)
    print("TEST 5: Decreasing Wavelength Array")
    print("="*70)
    
    # Wavelengths in decreasing order (some instruments output this way)
    wavelengths = np.array([1.0, 0.9, 0.8, 0.7, 0.6, 0.5])
    target = 0.75
    
    print(f"Wavelengths (decreasing): {wavelengths}")
    print(f"Target: {target}")
    
    idx, wl, diff = find_closest_wavelength_index(target, wavelengths)
    
    print(f"Result: index={idx}, wavelength={wl}, difference={diff}")
    
    # Should work correctly with decreasing array
    assert idx in [2, 3], f"Expected index 2 or 3, got {idx}"
    
    print("✓ PASSED - Works with decreasing arrays\n")
    return True


def test_different_units():
    """Test with different wavelength units."""
    print("\n" + "="*70)
    print("TEST 6: Different Wavelength Units")
    print("="*70)
    
    # Test in nanometers instead of microns
    wavelengths_nm = np.array([400, 450, 500, 550, 600, 650, 700])  # nm
    target_nm = 532  # nm (green laser)
    
    print(f"Wavelengths (nm): {wavelengths_nm}")
    print(f"Target: {target_nm} nm")
    
    idx, wl, diff = find_closest_wavelength_index(target_nm, wavelengths_nm)
    
    print(f"Result: index={idx}, wavelength={wl} nm, difference={diff} nm")
    
    # Should find 550 nm or 500 nm as closest
    assert wl in [500, 550], f"Expected 500 or 550 nm, got {wl}"
    
    print("✓ PASSED - Works with different units\n")
    return True


def run_all_tests():
    """Run all tests."""
    print("\n" + "="*70)
    print("TESTING find_closest_wavelength_index() FUNCTION")
    print("="*70)
    
    tests = [
        ("Basic Functionality", test_basic_functionality),
        ("Exact Match", test_exact_match),
        ("User's Wavelength Array", test_user_wavelength_array),
        ("Edge Cases", test_edge_cases),
        ("Decreasing Array", test_decreasing_wavelengths),
        ("Different Units", test_different_units),
    ]
    
    results = []
    for name, test_func in tests:
        try:
            result = test_func()
            results.append((name, result))
        except Exception as e:
            print(f"\n✗ TEST FAILED: {name}")
            print(f"Error: {e}")
            import traceback
            traceback.print_exc()
            results.append((name, False))
    
    # Summary
    print("\n" + "="*70)
    print("TEST SUMMARY")
    print("="*70)
    
    passed = sum(1 for _, result in results if result)
    total = len(results)
    
    for name, result in results:
        status = "✓ PASSED" if result else "✗ FAILED"
        print(f"{status}: {name}")
    
    print(f"\nTotal: {passed}/{total} tests passed")
    
    if passed == total:
        print("\n🎉 ALL TESTS PASSED! 🎉")
        return True
    else:
        print(f"\n⚠️  {total - passed} test(s) failed")
        return False


if __name__ == "__main__":
    success = run_all_tests()
    sys.exit(0 if success else 1)
