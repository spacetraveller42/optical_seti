# How to Use the Complete Gaussian Generator Module

## 📋 Quick Reference

### Getting the Module

**File:** `GAUSSIAN_GENERATOR_COMPLETE_MODULE.py`

**To use:**
1. Open `GAUSSIAN_GENERATOR_COMPLETE_MODULE.py` in GitHub
2. Click "Raw" to see plain text
3. Select all (Ctrl+A / Cmd+A)
4. Copy (Ctrl+C / Cmd+C)
5. Paste into your IDE as `gaussian_generator.py`

OR just download the file directly from this PR.

---

## 🚀 Your Use Case: Inject Laser at 0.532 Microns

```python
import numpy as np
from gaussian_generator import add_gaussian_to_array

# Your wavelength array (151 points, 0.2-1.995 microns)
wavelengths = np.array([0.2, 0.2285714286, 0.2612244898, ..., 1.995336705])

# Your flux data
flux = np.array([...])  # Your flux values

# Inject laser at 0.532 microns
result = add_gaussian_to_array(
    flux,
    fwhm=3.0,              # Width in array indices
    amplitude=0.05,         # Peak height
    wavelength=0.532,       # ← Use wavelength parameter!
    wavelengths=wavelengths # ← Must provide wavelengths array
)
```

**KEY POINT:** Use `wavelength=0.532` (NOT `center=0.532`)!

---

## 📖 Complete Usage Examples

### Example 1: Index-Based Center (Traditional)

```python
import numpy as np
from gaussian_generator import generate_gaussian

# Generate Gaussian at array position 50
gaussian = generate_gaussian(
    fwhm=10.0,
    amplitude=100.0,
    center=50.0,      # Array index (0-99)
    array_length=100
)
```

### Example 2: Wavelength-Based Center (NEW!)

```python
import numpy as np
from gaussian_generator import generate_gaussian

# Your wavelength array
wavelengths = np.linspace(400, 700, 100)  # 400-700 nm

# Generate Gaussian at 550 nm
gaussian = generate_gaussian(
    fwhm=10.0,
    amplitude=100.0,
    wavelength=550.0,       # Wavelength value
    array_length=100,
    wavelengths=wavelengths # Required for wavelength mode
)
```

### Example 3: Wavelength Between Array Elements

```python
# Inject at 525.5 nm (even if it's not in your array!)
gaussian = generate_gaussian(
    fwhm=5.0,
    amplitude=100.0,
    wavelength=525.5,       # Between array elements - will interpolate!
    array_length=100,
    wavelengths=wavelengths
)
```

### Example 4: Using Area Instead of Amplitude

```python
# Specify area under curve instead of peak height
gaussian = generate_gaussian(
    fwhm=10.0,
    area=1000.0,            # Total area under curve
    wavelength=550.0,
    array_length=100,
    wavelengths=wavelengths
)
```

### Example 5: Add to Existing Data

```python
import numpy as np
from gaussian_generator import add_gaussian_to_array

# Your data
data = np.ones(100) * 500

# Add Gaussian
result = add_gaussian_to_array(
    data,
    fwhm=10.0,
    amplitude=200.0,
    wavelength=550.0,
    wavelengths=wavelengths
)
```

### Example 6: Multi-Dimensional Arrays

```python
# 2D array (50 rows x 100 columns)
data_2d = np.ones((50, 100)) * 500

# Add Gaussian along axis 1 (broadcasts across all rows)
result = add_gaussian_to_array(
    data_2d,
    fwhm=10.0,
    amplitude=200.0,
    center=50.0,
    axis=1  # Apply to columns
)
```

---

## 🎯 Parameter Reference

### Required Parameters

| Parameter | Type | Description |
|-----------|------|-------------|
| `fwhm` | float | Full Width Half Maximum (in array indices) |
| `array_length` | int | Length of output array |

### Center Specification (Pick ONE)

| Parameter | Type | Description | When to Use |
|-----------|------|-------------|-------------|
| `center` | float | Array index (0 to N-1) | You know the position |
| `wavelength` | float | Physical wavelength | You know the wavelength |

**IMPORTANT:** If you use `wavelength`, you MUST also provide `wavelengths` array!

### Amplitude Specification (Pick ONE)

| Parameter | Type | Description |
|-----------|------|-------------|
| `amplitude` | float | Peak height of Gaussian |
| `area` | float | Total area under curve |

### Optional Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `wavelengths` | array | None | Required if using `wavelength` |
| `axis` | int | -1 | Axis for multi-dimensional arrays |
| `array_length` | int | None | For different Gaussian length |

---

## ⚠️ Common Mistakes

### Mistake 1: Using center with wavelength value

```python
# ❌ WRONG
result = add_gaussian_to_array(flux, fwhm=3.0, amplitude=0.05, center=0.532)
# This treats 0.532 as array INDEX → injects at edge!

# ✅ CORRECT
result = add_gaussian_to_array(flux, fwhm=3.0, amplitude=0.05, 
                               wavelength=0.532, wavelengths=wavelengths)
# This treats 0.532 as WAVELENGTH → injects at correct position!
```

### Mistake 2: Forgetting wavelengths array

```python
# ❌ WRONG
gaussian = generate_gaussian(fwhm=10.0, amplitude=100, wavelength=550.0, array_length=100)
# Error: wavelengths array required!

# ✅ CORRECT
gaussian = generate_gaussian(fwhm=10.0, amplitude=100, wavelength=550.0, 
                             array_length=100, wavelengths=wavelengths)
```

### Mistake 3: Unit mismatch

```python
# ❌ WRONG - mixing units
wavelengths = [0.4, 0.5, 0.6, 0.7]  # microns
gaussian = generate_gaussian(..., wavelength=550, ...)  # nanometers!

# ✅ CORRECT - consistent units
wavelengths = [0.4, 0.5, 0.6, 0.7]  # microns
gaussian = generate_gaussian(..., wavelength=0.550, ...)  # microns
```

---

## 🔍 Two Parameters, Two Meanings

| You Know | Use This | Example |
|----------|----------|---------|
| Array **position** | `center=50` | 50th element |
| Physical **wavelength** | `wavelength=0.532` | 0.532 microns |

**The module will warn you if you accidentally mix them up!**

---

## 📊 Module Features

✅ **FWHM-based generation** - Specify width at half maximum  
✅ **Amplitude or area** - Choose peak height or total area  
✅ **Index or wavelength center** - Flexible positioning  
✅ **Wavelength interpolation** - Works between array elements  
✅ **Broadcasting** - Multi-dimensional array support  
✅ **Length flexibility** - Different Gaussian/data lengths  
✅ **Input validation** - Catches errors with helpful messages  
✅ **Warning system** - Helps prevent parameter confusion  

---

## 🧪 Running Tests

```python
# If you saved as gaussian_generator.py
python gaussian_generator.py

# Will run:
# 1. Built-in test suite
# 2. Comprehensive demonstration
# 3. Generate example plots
```

---

## 📞 Quick Help

**Q: My laser appears at the edge, not where I want!**  
A: You're probably using `center=0.532` instead of `wavelength=0.532`. See "Common Mistakes" above.

**Q: Can I inject at wavelengths not in my array?**  
A: YES! The module interpolates - works perfectly for wavelengths between array elements.

**Q: What units should I use?**  
A: Any units, as long as they're CONSISTENT. If wavelengths are in microns, use wavelength in microns. If wavelengths are in nm, use wavelength in nm.

**Q: How do I specify total signal strength?**  
A: Use `area` parameter instead of `amplitude`.

---

## 📦 Complete Module Contents

**Main Functions:**
- `generate_gaussian()` - Core generation function
- `add_gaussian_to_array()` - Add to existing data

**Helper Functions:**
- `plot_gaussian_comparison()` - Visualization
- `run_tests()` - Test suite
- `run_demo()` - Demonstration

**Constants:**
- `FWHM_TO_SIGMA_FACTOR` - Conversion factor

**Total:** ~640 lines of clean, documented code

---

## ✅ Ready to Use!

Just copy `GAUSSIAN_GENERATOR_COMPLETE_MODULE.py` into your IDE and you're good to go! 🎯
