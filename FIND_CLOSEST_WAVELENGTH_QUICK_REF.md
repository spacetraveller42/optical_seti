# Quick Reference: find_closest_wavelength_index()

## Purpose
Find the index of the wavelength in an array that's closest to a target wavelength.

## Function Signature
```python
find_closest_wavelength_index(target_wavelength, wavelengths)
```

## Parameters
- **target_wavelength** (float): The wavelength you're searching for
- **wavelengths** (array-like): Your wavelength array

## Returns
A tuple of three values:
1. **index** (int): Index of closest wavelength
2. **closest_wavelength** (float): Actual wavelength value at that index  
3. **difference** (float): Absolute difference (target - closest)

## Quick Example
```python
from gaussian_generator import find_closest_wavelength_index
import numpy as np

# Your wavelength array
wavelengths = np.array([0.2, 0.3, 0.4, 0.5, 0.6, 0.7])

# Find closest to 0.532 microns
idx, wl, diff = find_closest_wavelength_index(0.532, wavelengths)

print(f"Index: {idx}")              # Index: 4
print(f"Wavelength: {wl}")          # Wavelength: 0.6
print(f"Difference: {diff}")        # Difference: 0.068
```

## Your Use Case
```python
# Your 151-point wavelength array (0.2-1.995 microns)
wavelengths = np.array([0.2, 0.2285714286, ..., 1.995336705])

# Find where 0.532 micron laser will be
idx, wl, diff = find_closest_wavelength_index(0.532, wavelengths)

print(f"Green laser (0.532 μm) closest to:")
print(f"  Index: {idx}")              # 13
print(f"  Wavelength: {wl:.6f} μm")  # 0.533657 μm
print(f"  Difference: {diff*1000:.2f} nm")  # 1.66 nm

# Extract flux at that wavelength
flux_at_laser = flux_array[idx]
```

## Common Use Cases

### 1. Find injection position
```python
idx, wl, diff = find_closest_wavelength_index(0.532, wavelengths)
print(f"Gaussian will peak at index {idx}")
```

### 2. Extract flux at wavelength
```python
idx, wl, diff = find_closest_wavelength_index(0.589, wavelengths)
flux_at_sodium = flux_array[idx]
```

### 3. Find multiple calibration lines
```python
lines = [0.532, 0.633, 0.780]  # Green, Red, NIR lasers
for target in lines:
    idx, wl, diff = find_closest_wavelength_index(target, wavelengths)
    print(f"{target:.3f} μm → index {idx}, actual {wl:.6f} μm")
```

### 4. Check if wavelength is in range
```python
target = 0.532
if target < wavelengths.min() or target > wavelengths.max():
    print(f"Out of range!")
else:
    idx, wl, diff = find_closest_wavelength_index(target, wavelengths)
    print(f"In range: index {idx}")
```

### 5. Use with Gaussian injection
```python
# Find where it will be
idx, wl, diff = find_closest_wavelength_index(0.532, wavelengths)

# Inject Gaussian
result = add_gaussian_to_array(flux, fwhm=3.0, amplitude=0.05,
                               wavelength=0.532, wavelengths=wavelengths)

# Verify peak is at expected location
assert np.argmax(result) == idx
```

## Features

✅ Works with any wavelength units (microns, nm, Angstroms)  
✅ Handles both increasing and decreasing arrays  
✅ Fast (O(n) with numpy)  
✅ Returns comprehensive information  
✅ Validates input (checks for empty, non-monotonic, NaN)  
✅ Clear error messages  

## Error Handling

```python
# Empty array
try:
    idx, wl, diff = find_closest_wavelength_index(0.5, [])
except ValueError as e:
    print(e)  # "wavelengths array is empty"

# Non-monotonic array
try:
    idx, wl, diff = find_closest_wavelength_index(0.5, [0.5, 0.7, 0.6])
except ValueError as e:
    print(e)  # "wavelengths array must be monotonically increasing or decreasing"
```

## Important Notes

- Returns **nearest neighbor** index (does not interpolate)
- For interpolation between points, use `wavelength` parameter in `generate_gaussian()` or `add_gaussian_to_array()`
- Both target and array must use same units
- Works with decreasing wavelength arrays (some instruments output this way)

## See Also

- `test/test_find_closest_wavelength.py` - Comprehensive test suite
- `test/example_find_closest_wavelength.py` - 5 usage examples
- Function docstring - Full API documentation
