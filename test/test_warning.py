#!/usr/bin/env python
"""
Test the new warning system
"""

import sys
import os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import warnings
import numpy as np
from gaussian_generator import add_gaussian_to_array

# Enable all warnings
warnings.simplefilter("always")

wavelengths = np.linspace(0.2, 2.0, 151)
flux = np.ones(151) * 0.08

print("Testing with center=0.532 (should trigger warning if wavelengths are passed):")
print()

try:
    # This should trigger the warning
    result = add_gaussian_to_array(
        flux,
        fwhm=3.0,
        amplitude=0.05,
        center=0.532,
        wavelengths=wavelengths  # Passing wavelengths even though using center!
    )
    print("Injection succeeded (with warning above)")
except Exception as e:
    print(f"Error: {e}")
