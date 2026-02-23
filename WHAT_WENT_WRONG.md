# Your 0.532 Micron Laser Injection - What Went Wrong?

## TL;DR: **The Code Works! Here's How to Use It**

Your wavelength injection **is working correctly**. The 0.532 micron (532 nm) laser was successfully injected at the right position.

---

## Quick Start - Copy This Code

```python
import numpy as np
from gaussian_generator import add_gaussian_to_array

# Your wavelength array (in microns!)
wavelengths = np.array([0.2, 0.2285714286, ..., 1.995336705])

# Your flux data
flux = your_flux_array  # Replace with your actual data

# Inject 0.532 micron laser (532 nm green laser)
result = add_gaussian_to_array(
    flux,
    fwhm=3.0,           # Width in pixels
    amplitude=0.05,     # Peak height (adjust to be visible)
    wavelength=0.532,   # ← TARGET IN MICRONS
    wavelengths=wavelengths  # ← YOUR WAVELENGTH ARRAY
)

# Done! The laser is now in your spectrum.
```

---

## What Actually Happened

When you inject at **0.532 microns**:
- Target: 0.532000 microns (532.00 nm)
- Peak appears at: 0.533657 microns (533.66 nm)
- **Difference: 1.66 nm** ← This is CORRECT!

### Why the 1.66 nm difference?

Your wavelength array has **discrete points**:
```
Index 12: 0.5298726908 microns (529.87 nm)
Index 13: 0.5336574958 microns (533.66 nm) ← Peak ends up here
Target:   0.5320000000 microns (532.00 nm) ← Your request
```

The target (0.532) falls **between** indices 12 and 13. The Gaussian peak appears at the nearest grid point (index 13). This is **expected and correct** for discrete data!

---

## Most Common Mistakes

### 1. **Unit Mismatch** ⚠️ (Most Likely Issue!)

```python
# ❌ WRONG - Mixing units!
wavelengths_nm = np.array([500, 510, 520, 530, 540])  # nanometers
result = add_gaussian_to_array(..., wavelength=0.532,  # microns!
                               wavelengths=wavelengths_nm)

# ✓ CORRECT - Both in microns
wavelengths_microns = np.array([0.50, 0.51, 0.52, 0.53, 0.54])  
result = add_gaussian_to_array(..., wavelength=0.532,  # microns
                               wavelengths=wavelengths_microns)

# ✓ ALSO CORRECT - Both in nanometers
wavelengths_nm = np.array([500, 510, 520, 530, 540])
result = add_gaussian_to_array(..., wavelength=532,  # nanometers!
                               wavelengths=wavelengths_nm)
```

**YOUR DATA: Your wavelengths are in MICRONS, so use 0.532 (not 532)!**

### 2. **Amplitude Too Small**

```python
# ❌ WRONG - Too small to see
flux = np.ones(100) * 0.08  # Baseline at 0.08
result = add_gaussian_to_array(..., amplitude=0.001)  # Only adds 0.001!

# ✓ CORRECT - Visible amplitude
result = add_gaussian_to_array(..., amplitude=0.05)  # Adds 0.05 (62% increase)
```

### 3. **Missing wavelengths Parameter**

```python
# ❌ WRONG - Forgot wavelengths array!
result = add_gaussian_to_array(flux, fwhm=3.0, amplitude=0.05,
                               wavelength=0.532)  # Error: missing wavelengths!

# ✓ CORRECT
result = add_gaussian_to_array(flux, fwhm=3.0, amplitude=0.05,
                               wavelength=0.532,
                               wavelengths=wavelengths)  # ← Don't forget this!
```

---

## How to Check If It Worked

```python
# After injection
peak_idx = np.argmax(result)
peak_wavelength = wavelengths[peak_idx]

print(f"Peak at: {peak_wavelength:.6f} microns ({peak_wavelength*1000:.2f} nm)")
print(f"Target: 0.532 microns (532.00 nm)")
print(f"Difference: {abs(peak_wavelength - 0.532)*1000:.2f} nm")

# If difference is < 5 nm, it worked!
if abs(peak_wavelength - 0.532) * 1000 < 5.0:
    print("✓ SUCCESS - Laser injected at correct wavelength!")
else:
    print("⚠ Warning - Check your parameters")
```

---

## Visualization

```python
import matplotlib.pyplot as plt

# Inject laser
result = add_gaussian_to_array(flux, fwhm=3.0, amplitude=0.05,
                               wavelength=0.532, wavelengths=wavelengths)

# Plot
plt.figure(figsize=(12, 5))
plt.plot(wavelengths, flux, 'gray', alpha=0.5, label='Original')
plt.plot(wavelengths, result, 'b-', linewidth=2, label='With 532 nm laser')
plt.axvline(0.532, color='r', linestyle='--', label='Target (532 nm)')
plt.xlabel('Wavelength (microns)')
plt.ylabel('Flux')
plt.title('Laser Injection at 532 nm')
plt.legend()
plt.grid(True, alpha=0.3)
plt.savefig('my_laser_injection.png')
plt.show()
```

---

## Still Not Working?

### Run the Diagnostic:

```bash
cd /path/to/optical_seti
python test/test_user_wavelength_issue.py
```

This will show you **exactly** what's happening with your data.

### Check These:

1. **Units:**
   ```python
   print(f"Min wavelength: {wavelengths.min()}")
   print(f"Max wavelength: {wavelengths.max()}")
   # If these are 400-700 → nanometers
   # If these are 0.4-0.7 → microns
   ```

2. **Array Lengths:**
   ```python
   print(f"Flux length: {len(flux)}")
   print(f"Wavelength length: {len(wavelengths)}")
   # These MUST match!
   ```

3. **Target in Range:**
   ```python
   if 0.532 < wavelengths.min() or 0.532 > wavelengths.max():
       print("ERROR: 0.532 is out of range!")
   ```

---

## Your Exact Working Code

Based on your wavelength array (0.2 to 1.995 microns, 151 points):

```python
import numpy as np
from gaussian_generator import add_gaussian_to_array

# YOUR wavelength array (in microns)
wavelengths = np.array([
    0.2, 0.2285714286, 0.2612244898, 0.2985422741, 0.3411911703, 
    # ... (all 151 points) ...
    1.939525859, 1.967233371, 1.995336705
])

# YOUR flux data (use your actual flux array here)
flux = np.ones(151) * 0.08  # Replace with your actual data

# INJECT THE LASER
result = add_gaussian_to_array(
    flux,
    fwhm=3.0,        # Adjust width as needed
    amplitude=0.05,  # Adjust brightness as needed
    wavelength=0.532,  # 532 nm = 0.532 microns
    wavelengths=wavelengths
)

# VERIFY
peak_idx = np.argmax(result)
print(f"✓ Laser injected at index {peak_idx}")
print(f"  Wavelength: {wavelengths[peak_idx]:.6f} microns")
print(f"  Expected: 0.532 microns")
print(f"  Difference: {abs(wavelengths[peak_idx] - 0.532)*1000:.2f} nm")
```

**This WILL work with your data!** The diagnostic test confirms it.

---

## Example Output You Should See

```
✓ Laser injected at index 13
  Wavelength: 0.533657 microns
  Expected: 0.532 microns
  Difference: 1.66 nm
```

The 1.66 nm difference is **normal and expected**!

---

## Need More Help?

Run these diagnostic scripts:
- `python test/test_user_wavelength_issue.py` - Full diagnostic
- `python test/example_user_injection.py` - Simple working example

Check the troubleshooting guide:
- `TROUBLESHOOTING_WAVELENGTH_INJECTION.md`

The code is working correctly. If you're still having issues, the problem is likely:
1. Unit mismatch (microns vs nanometers)
2. Amplitude too small to see
3. Missing or incorrect parameters

**Double-check your units!** Your wavelengths are in MICRONS, so use `wavelength=0.532` (not 532).
