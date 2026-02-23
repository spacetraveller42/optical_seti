# Gaussian Generator - Finalized Standalone Module

## 📄 File Information

**Filename:** `gaussian_generator.py`  
**Size:** 14 KB  
**Location:** `/home/runner/work/optical_seti/optical_seti/gaussian_generator.py`

## ✅ What This File Does

This is a **complete, standalone Python module** that provides:

1. ✅ **Set FWHM directly** - Specify Full Width Half Maximum as a parameter
2. ✅ **Get Gaussian as array** - Returns numpy array for direct manipulation
3. ✅ **Add to existing arrays** - Simple array addition with `+` operator
4. ✅ **Visible in plots** - Gaussians are clearly visible when plotted

## 🎯 Core Functions

### 1. `generate_gaussian(fwhm, amplitude, center, array_length)`
Generates a Gaussian curve with specified FWHM and returns it as a numpy array.
- **NEW:** Now supports `wavelength` parameter for wavelength-based injection
- Supports both `amplitude` and `area` parameters

### 2. `add_gaussian_to_array(data, fwhm, amplitude, center)`
Convenience function that adds a Gaussian to existing data in one step.
- **NEW:** Supports wavelength-based injection with `wavelength` and `wavelengths` parameters
- Works with wavelengths between grid points (interpolates automatically)

### 3. `plot_gaussian_comparison(data, fwhm, amplitude, center, title)`
Utility function to visualize original vs. modified data.

## 💡 Quick Start Examples

### Example 1: Basic Usage
```python
import numpy as np
from gaussian_generator import generate_gaussian

# Generate a Gaussian with FWHM=10
gaussian = generate_gaussian(fwhm=10.0, amplitude=100, center=50, array_length=100)

# Add to your data
your_data = np.random.randn(100) * 5 + 500
result = your_data + gaussian
```

### Example 2: Multiple Gaussians
```python
from gaussian_generator import add_gaussian_to_array

# Start with baseline
data = np.ones(500) * 1000

# Add narrow laser signal
data = add_gaussian_to_array(data, fwhm=5.0, amplitude=300, center=200)

# Add broader emission line
data = add_gaussian_to_array(data, fwhm=15.0, amplitude=150, center=350)
```

### Example 3: Wavelength-Based Injection (NEW!)
```python
from gaussian_generator import add_gaussian_to_array
import numpy as np

# Load or create wavelength array (e.g., from HARPS spectral data)
wavelengths = np.linspace(400, 700, 1000)  # 400-700 nm
flux = np.ones(1000) * 1000  # Baseline spectrum

# Inject Gaussian at specific wavelength (not limited to grid points!)
flux_with_spike = add_gaussian_to_array(
    flux,
    fwhm=5.0,
    amplitude=500.0,
    wavelength=532.0,      # Green laser at 532.0 nm
    wavelengths=wavelengths
)

# Works with wavelengths BETWEEN grid points too!
flux_with_spike = add_gaussian_to_array(
    flux,
    fwhm=3.0,
    area=1000.0,           # Can use area instead of amplitude
    wavelength=589.5,      # Between grid points - will interpolate!
    wavelengths=wavelengths
)
```

### Example 4: With Visualization
```python
from gaussian_generator import plot_gaussian_comparison
import matplotlib.pyplot as plt

# Your spectral data
spectral_data = np.random.randn(200) * 10 + 500

# Plot comparison
fig, axes = plot_gaussian_comparison(
    spectral_data, 
    fwhm=12.0, 
    amplitude=200.0, 
    center=100.0
)
plt.show()
```

## 🚀 How to Use on Your Laptop

### Option 1: Copy the File
1. Download `gaussian_generator.py` from the repository
2. Save it to your project directory
3. Import and use:
```python
from gaussian_generator import generate_gaussian, add_gaussian_to_array
```

### Option 2: Run Standalone
The file includes built-in tests and demonstration:
```bash
python gaussian_generator.py
```

This will:
- ✅ Run automated tests
- ✅ Show demonstration examples
- ✅ Create visualization plots
- ✅ Confirm all functionality works

## 📦 Dependencies

Only requires standard scientific Python packages:
- `numpy` - Array operations
- `matplotlib` - Plotting (optional, only for visualization functions)

Install with:
```bash
pip install numpy matplotlib
```

## 🔬 Technical Details

- **FWHM Conversion:** Uses exact formula: FWHM = 2 × √(2 × ln(2)) × σ ≈ 2.35482 × σ
- **Input Validation:** All functions validate inputs and provide clear error messages
- **Array-Based:** Returns standard numpy arrays for maximum compatibility
- **Well-Documented:** Complete docstrings with examples for all functions

## 📊 What You Get

When you run the standalone file, it generates:
- `gaussian_demo_single.png` - Single Gaussian addition visualization
- `gaussian_demo_multiple.png` - Multiple Gaussians visualization

## ✨ Features

- ✅ **Standalone** - No dependencies on other project files
- ✅ **Wavelength-Based Injection** - Inject at specific wavelengths (not just array indices)
- ✅ **Interpolation Support** - Works with wavelengths between grid points
- ✅ **Tested** - Built-in test suite confirms functionality
- ✅ **Documented** - Comprehensive docstrings and examples
- ✅ **Production-Ready** - Input validation and error handling
- ✅ **Easy to Use** - Simple, intuitive API
- ✅ **Broadcasting Support** - Works with multi-dimensional arrays

## 📝 Complete Function Signatures

```python
# Index-based (traditional)
generate_gaussian(fwhm, amplitude, center, array_length) -> numpy.ndarray

# Wavelength-based (NEW!)
generate_gaussian(fwhm, amplitude, wavelength, array_length, wavelengths) -> numpy.ndarray

# Also supports area parameter instead of amplitude
generate_gaussian(fwhm, area, center, array_length) -> numpy.ndarray

# Add to existing array (index-based)
add_gaussian_to_array(data, fwhm, amplitude, center) -> numpy.ndarray

# Add to existing array (wavelength-based)
add_gaussian_to_array(data, fwhm, amplitude, wavelength, wavelengths) -> numpy.ndarray

plot_gaussian_comparison(data, fwhm, amplitude, center, title) -> (figure, axes)

run_demo() -> None  # Built-in demonstration

run_tests() -> bool  # Built-in test suite
```

## 🎓 Use Cases

Perfect for:
- Testing spike detection algorithms
- Creating synthetic spectral data
- Simulating laser signals
- Adding emission lines to spectra
- Educational demonstrations
- Algorithm validation

## ✅ Final Status

**COMPLETE AND READY TO USE**

All requirements met:
- ✅ Set FWHM directly
- ✅ Get Gaussian as array
- ✅ Add to existing arrays
- ✅ Visible when plotted
- ✅ Fully tested
- ✅ Standalone module
- ✅ Ready to save to laptop

---

**File Location:** `/home/runner/work/optical_seti/optical_seti/gaussian_generator.py`  
**Status:** ✅ Finalized and ready for export
