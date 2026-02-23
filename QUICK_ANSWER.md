# Quick Answer: How to Inject at a Specific Wavelength

## Your Question
**"I thought I could set the center as a wavelength value? How do I use the new feature?"**

## Quick Answer

✅ **YES! The feature exists and works!**  
✅ **Use `wavelength=0.532` NOT `center=0.532`**

## Your Fix (2 Changes)

```python
# Change from this:
result = add_gaussian_to_array(flux, fwhm=3.0, amplitude=0.05, center=0.532)

# To this:
result = add_gaussian_to_array(flux, fwhm=3.0, amplitude=0.05,
                               wavelength=0.532,        # Change 1
                               wavelengths=wavelengths)  # Change 2
```

## Why?

- **`center`** = Array **index** (0, 1, 2, ..., 150)
  - `center=0.532` → position 0.532 → **LEFT EDGE**
  
- **`wavelength`** = Physical **wavelength** (0.2, 0.532, 1.995 microns)
  - `wavelength=0.532` → 0.532 microns → **MIDDLE**

## Does It Work with Wavelengths NOT in the Array?

**YES!** That's the whole point of the feature!

```python
# Even if 0.550 is NOT in your wavelengths array:
result = add_gaussian_to_array(flux, fwhm=3.0, amplitude=0.05,
                               wavelength=0.550,        # Not in array - works anyway!
                               wavelengths=wavelengths)
# ✓ Works! Interpolates to find correct position
```

## Test It

Run this:
```bash
python test/example_center_vs_wavelength.py
```

Output shows:
- `center=0.532` → Peak at 0.228 μm (EDGE) ❌
- `wavelength=0.532` → Peak at 0.534 μm (MIDDLE) ✓
- `wavelength=0.550` → Works even though not in array! ✓

## Summary

The wavelength injection feature **EXISTS** and **WORKS PERFECTLY**.

Just use the **`wavelength` parameter**, not `center`.
