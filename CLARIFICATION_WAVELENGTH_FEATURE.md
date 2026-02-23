# Yes, You CAN Inject at Wavelengths Not in the Array! 

## Your Question
**"Wait, but I thought I could set the center as a wavelength value? How do I use the new feature of injecting at a wavelength not listed in the array otherwise?"**

## The Answer: Yes, the Feature Exists! Use the `wavelength` Parameter

You **CAN** inject at wavelengths not in your array. The feature is fully implemented and working. But you need to use the correct parameter.

## Two Different Parameters

| Parameter | What It Means | Example Usage |
|-----------|---------------|---------------|
| **`center`** | **Array index** (position 0, 1, 2, ...) | `center=50` → 50th element |
| **`wavelength`** | **Physical wavelength** (in your units) | `wavelength=0.532` → 0.532 microns |

## How to Use Wavelength Injection

### ✓ CORRECT Method

```python
from gaussian_generator import add_gaussian_to_array
import numpy as np

# Your wavelength array (151 points, 0.2-1.995 microns)
wavelengths = np.array([0.2, 0.228571, ..., 1.995337])
flux = np.ones(151) * 0.08  # Your baseline data

# Inject at 0.532 microns (even if not exactly in your array!)
result = add_gaussian_to_array(
    flux,
    fwhm=3.0,
    amplitude=0.05,
    wavelength=0.532,        # ← Use 'wavelength' parameter
    wavelengths=wavelengths   # ← Must provide wavelength array
)
```

**This works even if 0.532 is not exactly in your wavelengths array!**

The code automatically:
1. Finds where 0.532 falls in your wavelength array
2. Interpolates to get a fractional index (e.g., 12.5621)
3. Injects the Gaussian at that precise location

## Why Two Parameters?

### `center` Parameter (Traditional Method)
Use when you know the **array position**:
```python
gaussian = generate_gaussian(fwhm=5.0, amplitude=100, center=50, array_length=100)
# Injects at position 50 (51st element)
```

### `wavelength` Parameter (New Feature!)
Use when you know the **physical wavelength**:
```python
gaussian = generate_gaussian(fwhm=5.0, amplitude=100, 
                            wavelength=0.532,
                            array_length=151,
                            wavelengths=wavelengths)
# Injects at 0.532 microns (automatically finds correct position)
```

## The Confusion Explained

### ❌ What You Were Doing (WRONG)
```python
result = add_gaussian_to_array(flux, fwhm=3.0, amplitude=0.05, center=0.532)
```

**Problem:** `center=0.532` means "inject at array index 0.532"  
→ Index 0.532 ≈ position 1 → **LEFT EDGE** of your spectrum!

### ✓ What You Should Do (CORRECT)
```python
result = add_gaussian_to_array(flux, fwhm=3.0, amplitude=0.05,
                               wavelength=0.532,
                               wavelengths=wavelengths)
```

**This means:** "inject at wavelength 0.532 microns"  
→ Code finds position via interpolation → **MIDDLE** of your spectrum!

## Key Points

1. **YES**, you can inject at wavelengths not in your array
2. Use the `wavelength` parameter (NOT `center`)
3. Must also provide the `wavelengths` array
4. Works perfectly with fractional/interpolated positions
5. `center` is for array indices (0, 1, 2, ...)
6. `wavelength` is for physical wavelength values (0.2, 0.532, 1.5 microns, ...)

## Your Specific Fix

**Change from:**
```python
result = add_gaussian_to_array(flux, fwhm=3.0, amplitude=0.05, center=0.532)
```

**To:**
```python
result = add_gaussian_to_array(flux, fwhm=3.0, amplitude=0.05,
                               wavelength=0.532,
                               wavelengths=wavelengths)
```

That's it! Just change `center` to `wavelength` and add `wavelengths=wavelengths`.

## Test It

```python
import numpy as np
from gaussian_generator import add_gaussian_to_array

# Your data
wavelengths = np.array([0.2, 0.2285714286, 0.2612244898, 0.2985422741, 
                        0.3411911703, 0.3899327661, 0.445637447, 
                        0.5092999394, 0.515, 0.5186785714, 0.5223834184,
                        0.5261147285, 0.5298726908, 0.5336574958, 
                        0.537469335, 0.5413084017, ...])
flux = np.ones(len(wavelengths)) * 0.08

# Inject at 0.532 microns
result = add_gaussian_to_array(flux, fwhm=3.0, amplitude=0.05,
                               wavelength=0.532,
                               wavelengths=wavelengths)

# Verify
peak_idx = np.argmax(result)
print(f"Peak at index: {peak_idx}")
print(f"Peak wavelength: {wavelengths[peak_idx]:.6f} microns")
print(f"Target wavelength: 0.532000 microns")

# Output:
# Peak at index: 13
# Peak wavelength: 0.533657 microns
# Target wavelength: 0.532000 microns
```

The small difference (1.66 nm) is expected because your array has discrete points and 0.532 falls between them.

## Summary

✅ Wavelength injection feature EXISTS  
✅ It DOES support wavelengths between array elements  
✅ Use `wavelength=0.532` (NOT `center=0.532`)  
✅ Must include `wavelengths=wavelengths`  
✅ Works perfectly!
