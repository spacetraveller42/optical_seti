# Summary: New find_closest_wavelength_index() Function

## What You Asked For

"I want another separate function added to the module, which can tell me the index of the element in the wavelength array with a wavelength value that's the closest to my specified value for the wavelength-defined center."

## What You Got

A new function `find_closest_wavelength_index()` that does exactly that!

## The Function

```python
def find_closest_wavelength_index(target_wavelength, wavelengths):
    """Find the index of the wavelength closest to a target wavelength."""
    # Returns: (index, closest_wavelength, difference)
```

## How To Use It

```python
from gaussian_generator import find_closest_wavelength_index
import numpy as np

# Your wavelength array (151 points, 0.2-1.995 microns)
wavelengths = np.array([0.2, 0.2285714286, ..., 1.995336705])

# Find closest wavelength to 0.532 microns
idx, wl, diff = find_closest_wavelength_index(0.532, wavelengths)

print(f"Closest wavelength to 0.532 μm:")
print(f"  Index: {idx}")              # 13
print(f"  Wavelength: {wl:.6f} μm")  # 0.533657 μm
print(f"  Difference: {diff*1000:.2f} nm")  # 1.66 nm
```

## What It Returns

Three values:
1. **index** (int) - The array index of the closest wavelength
2. **closest_wavelength** (float) - The actual wavelength value at that index
3. **difference** (float) - How far off the closest wavelength is from your target

## Real Example With Your Data

```python
>>> wavelengths = np.array([0.2, 0.2285714286, ..., 1.995336705])  # 151 points
>>> idx, wl, diff = find_closest_wavelength_index(0.532, wavelengths)
>>> print(f"Index: {idx}, Wavelength: {wl:.6f}, Diff: {diff*1000:.2f} nm")
Index: 13, Wavelength: 0.533657, Diff: 1.66 nm
```

## Where It's Located

The function is now in:
- ✅ `gaussian_generator.py` (main module)
- ✅ `gaussian_generator_COPYABLE.py` (copy-pasteable version)
- ✅ `GAUSSIAN_GENERATOR_COMPLETE_MODULE.py` (complete version)

## Test Results

All 6 tests pass ✅

Tested with your exact 151-point wavelength array:
- ✓ Finds index 13 for 0.532 microns
- ✓ Returns 0.533657 as closest wavelength
- ✓ Calculates 1.66 nm difference
- ✓ Works perfectly!

## Documentation

- **Quick Reference:** `FIND_CLOSEST_WAVELENGTH_QUICK_REF.md`
- **Test Suite:** `test/test_find_closest_wavelength.py` (6 tests)
- **Examples:** `test/example_find_closest_wavelength.py` (5 examples)
- **Function Docstring:** Complete API docs in code

## Common Uses

### 1. Find where Gaussian will peak
```python
idx, wl, diff = find_closest_wavelength_index(0.532, wavelengths)
print(f"Gaussian will peak at index {idx}")
```

### 2. Extract flux at specific wavelength
```python
idx, wl, diff = find_closest_wavelength_index(0.589, wavelengths)
flux_at_sodium = flux_array[idx]
```

### 3. Find multiple calibration lines
```python
for target in [0.532, 0.633, 0.780]:  # Green, Red, NIR
    idx, wl, diff = find_closest_wavelength_index(target, wavelengths)
    print(f"{target:.3f} μm → index {idx}")
```

## Features

✅ Simple 2-parameter interface  
✅ Returns 3 useful values  
✅ Works with any units (μm, nm, Å)  
✅ Fast (uses numpy)  
✅ Input validation  
✅ Clear error messages  
✅ Handles edge cases  
✅ Fully tested  
✅ Well documented  

## Status

✅ **Implemented**  
✅ **Tested** (6/6 tests passing)  
✅ **Documented** (Quick ref + examples + docstring)  
✅ **Ready to use**  

## How To Get It

**Option 1: Copy the updated module**
- Copy `gaussian_generator_COPYABLE.py`
- Paste into your IDE as `gaussian_generator.py`

**Option 2: Import from this repo**
```python
from gaussian_generator import find_closest_wavelength_index
```

That's it! The function is ready to use. 🎯
