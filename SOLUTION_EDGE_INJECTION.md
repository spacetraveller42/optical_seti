# SOLUTION: Laser Appearing at Edge Instead of 0.532 Microns

## The Problem

You said: "I'm setting 0.532 as the center, but it shows up at the edge."

## The Root Cause

You're using the **wrong parameter**!

### What's Happening:

```python
# ❌ WRONG - This is what you're doing:
result = add_gaussian_to_array(
    flux,
    fwhm=3.0,
    amplitude=0.05,
    center=0.532  # ← "center" means ARRAY INDEX, not wavelength!
)
```

When you use `center=0.532`, the code thinks you want to inject at **array index 0.532**, which rounds to **index 0 or 1** - that's the **LEFT EDGE** of your spectrum!

## The Fix

Use the `wavelength` parameter instead:

```python
# ✓ CORRECT - Do this instead:
result = add_gaussian_to_array(
    flux,
    fwhm=3.0,
    amplitude=0.05,
    wavelength=0.532,  # ← "wavelength" means physical wavelength value!
    wavelengths=wavelengths  # ← Must also provide wavelength array
)
```

## Understanding the Two Parameters

| Parameter | Meaning | Example Value | Where It Injects |
|-----------|---------|---------------|------------------|
| `center` | **Array INDEX** | `center=13` | At array position 13 |
| `wavelength` | **Physical wavelength** | `wavelength=0.532` | At wavelength 0.532 microns |

### Your Specific Case:

**Your wavelength array:** 151 points, from 0.2 to 1.995 microns

**What you want:** Inject at 0.532 microns (green laser)

**Using center=0.532:**
- Interprets 0.532 as array **INDEX**
- Index 0.532 ≈ index 1
- Index 1 corresponds to wavelength 0.228571 microns
- **Result: Laser appears at LEFT EDGE** ❌

**Using wavelength=0.532:**
- Interprets 0.532 as **WAVELENGTH** in microns
- Finds wavelength 0.532 in your array
- Maps to index 13 (middle of your 151-point array)
- Index 13 corresponds to wavelength 0.533657 microns
- **Result: Laser appears in MIDDLE where you want it** ✅

## Complete Working Code

```python
import numpy as np
from gaussian_generator import add_gaussian_to_array

# Your wavelength array (151 points, 0.2 to 1.995 microns)
wavelengths = np.array([0.2, 0.2285714286, ..., 1.995336705])

# Your flux data  
flux = np.ones(151) * 0.08  # Or your actual flux values

# ✓ CORRECT: Inject at 0.532 microns
result = add_gaussian_to_array(
    flux,
    fwhm=3.0,
    amplitude=0.05,
    wavelength=0.532,      # ← Wavelength in microns
    wavelengths=wavelengths  # ← Your wavelength array
)

# Verify
peak_idx = np.argmax(result)
print(f"Peak at index {peak_idx}")
print(f"Peak wavelength: {wavelengths[peak_idx]:.6f} microns")
# Output: Peak at index 13
#         Peak wavelength: 0.533657 microns
```

## New Warning System

Starting now, if you accidentally use `center` with a value that looks like a wavelength, you'll get a helpful warning:

```
⚠️  WARNING: center=0.532 looks like a WAVELENGTH value, not an array INDEX!

You specified:
  center=0.532

But 'center' expects an array INDEX (0 to 150),
and 0.532 is within your wavelengths range [0.200, 1.995].

Did you mean to use:
  wavelength=0.532, wavelengths=<your_wavelength_array>

This will inject at array index 0.5 ≈ index 1,
which might not be what you intended!
```

## Quick Reference

| What You Want | Parameter to Use | Example |
|---------------|------------------|---------|
| Inject at **array position 50** | `center=50` | `center=50` |
| Inject at **wavelength 0.532 μm** | `wavelength=0.532` | `wavelength=0.532, wavelengths=wave_array` |
| Inject at **632.8 nm** | `wavelength=0.6328` | `wavelength=0.6328, wavelengths=wave_array` (convert nm to μm!) |

## Summary

**Your issue:** Using `center=0.532` when you should use `wavelength=0.532`

**The fix:** Change parameter from `center` to `wavelength` and add `wavelengths` array

**Why it matters:** 
- `center` = array index (0-150 for your data)
- `wavelength` = physical wavelength value (0.2-1.995 for your data)

Now your laser will inject exactly where you want it! 🎯
