# YOUR FIX: Laser Appearing at Edge Instead of 0.532 Microns

## The Problem You Had

You said: **"I'm setting 0.532 as the center, but it shows up at the edge."**

## What Was Wrong

You were using:
```python
center=0.532
```

But you should have been using:
```python
wavelength=0.532, wavelengths=your_wavelength_array
```

## Why This Matters

These are **TWO COMPLETELY DIFFERENT THINGS**:

### `center` = Array Position (Index)
- **Means:** "Inject at array position X"
- **Example:** `center=50` means "inject at the 50th element"
- **Your case:** `center=0.532` means "inject at position 0.532" → rounds to position 0 or 1 → **LEFT EDGE!**

### `wavelength` = Physical Wavelength Value  
- **Means:** "Inject at wavelength X microns"
- **Example:** `wavelength=0.532` means "inject at 0.532 microns"
- **Your case:** `wavelength=0.532` means "inject at 0.532 μm" → middle of your spectrum → **WHERE YOU WANT IT!**

## The Fix (Just Change Your Code)

### ❌ WRONG (What you were doing):
```python
result = add_gaussian_to_array(
    flux,
    fwhm=3.0,
    amplitude=0.05,
    center=0.532  # ← This means "array index 0.532" = EDGE!
)
```

### ✅ CORRECT (What you should do):
```python
result = add_gaussian_to_array(
    flux,
    fwhm=3.0,
    amplitude=0.05,
    wavelength=0.532,      # ← This means "wavelength 0.532 μm" = MIDDLE!
    wavelengths=wavelengths  # ← Must include your wavelength array!
)
```

## That's It!

Just change `center` to `wavelength` and add your `wavelengths` array. Your laser will now appear exactly where you want it (at 0.532 microns, in the middle of your spectrum).

## Test Results

**Before (center=0.532):**
- Injected at: wavelength 0.228571 μm (index 1)
- Location: **LEFT EDGE** ❌

**After (wavelength=0.532):**
- Injected at: wavelength 0.533657 μm (index 13)
- Location: **MIDDLE** ✅

## Warning System

Now if you accidentally use `center` with a wavelength-like value, you'll get a helpful warning:

```
⚠️  WARNING: center=0.532 looks like a WAVELENGTH value, not an array INDEX!

Did you mean to use:
  wavelength=0.532, wavelengths=<your_wavelength_array>
```

## Summary

- **Problem:** Used `center=0.532` (array index)
- **Solution:** Use `wavelength=0.532, wavelengths=array` (physical wavelength)
- **Result:** Laser now injects at 0.532 microns where you want it!

Done! 🎯
